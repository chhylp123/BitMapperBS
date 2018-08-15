/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_update.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section LICENCE
 *
 * This file is part of pSAscan v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2015
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_UPDATE_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_UPDATE_H_INCLUDED

#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "../gap_buffer.h"
#include "../utils.h"
#include "inmem_gap_array.h"


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// This object creates a given number of threads that will perform gap array
// updates. Most of the time all threads are sleeping on a conditional variable.
// Once the gap buffer is available for processing, they are all woken up and
// perform the update in parallel. The caller then waits until all threads are
// finished and then puts the gap buffer in the poll of empty buffers.
//
// Only one object of this class should exist.
//==============================================================================
template<typename block_offset_type>
struct gap_parallel_updater {

  template<typename T>
  static void parallel_update(gap_parallel_updater<T> *updater, int id) {
    while (true) {
      // Wait until there is a gap buffer available or the
      // message 'no more buffers' arrives.
      std::unique_lock<std::mutex> lk(updater->m_avail_mutex);
      while (!(updater->m_avail[id]) && !(updater->m_avail_no_more))
        updater->m_avail_cv.wait(lk);

      if (!(updater->m_avail[id]) && updater->m_avail_no_more) {
        // No more buffers -- exit.
        lk.unlock();
        return;
      }

      updater->m_avail[id] = false;
      lk.unlock();

      // Safely perform the update.
      const gap_buffer<T> *buf = updater->m_buffer;
      inmem_gap_array *gap = updater->m_gap_array;
      int beg = buf->sblock_beg[id];
      int end = beg + buf->sblock_size[id];

      for (int i = beg; i < end; ++i) {
        T x = buf->m_content[i];
        gap->m_count[x]++;

        // Check if values wrapped-around.
        if (gap->m_count[x] == 0) {
          gap->m_excess_mutex.lock();
          gap->m_excess.push_back(x);
          gap->m_excess_mutex.unlock();
        }
      }

      // Update the number of finished threads.
      bool finished_last = false;
      std::unique_lock<std::mutex> lk2(updater->m_finished_mutex);
      updater->m_finished++;
      if (updater->m_finished == updater->m_threads_cnt)
        finished_last = true;
      lk2.unlock();

      // If this was the last thread finishing, let the caller know.
      if (finished_last)
        updater->m_finished_cv.notify_one();
    }
  }

  gap_parallel_updater(inmem_gap_array *gap_array, int threads_cnt)
      : m_gap_array(gap_array),
        m_threads_cnt(threads_cnt),
        m_avail_no_more(false) {
    m_avail = new bool[m_threads_cnt];
    std::fill(m_avail, m_avail + m_threads_cnt, false);
    m_threads = new std::thread*[m_threads_cnt];

    // After this, threads immediately hang up on m_avail_cv.
    for (int i = 0; i < m_threads_cnt; ++i)
      m_threads[i] = new std::thread(parallel_update<block_offset_type>, this, i);
  }

  ~gap_parallel_updater() {
    // Signal all threads to finish.
    std::unique_lock<std::mutex> lk(m_avail_mutex);
    m_avail_no_more = true;
    lk.unlock();
    m_avail_cv.notify_all();

    // Wait until all threads finish and release memory.
    for (int i = 0; i < m_threads_cnt; ++i) {
      m_threads[i]->join();
      delete m_threads[i];
    }
    delete[] m_threads;
    delete[] m_avail;
  }

  void update(const gap_buffer<block_offset_type> *buffer) {
    // Prepare a message for each thread that new buffer is available.
    std::unique_lock<std::mutex> lk(m_avail_mutex);
    m_finished = 0;
    m_buffer = buffer;
    for (int i = 0; i < m_threads_cnt; ++i)
      m_avail[i] = true;
    lk.unlock();

    // Wake up all threads to perform the update.
    m_avail_cv.notify_all();

    // Wait until all threads report that they are done.
    std::unique_lock<std::mutex> lk2(m_finished_mutex);
    while (m_finished != m_threads_cnt)
      m_finished_cv.wait(lk2);
    lk2.unlock();

    // We are done processing the buffer. The caller of this method
    // can now place the buffer into the poll of empty buffers.
  }

private:
  inmem_gap_array *m_gap_array;

  std::thread **m_threads;
  int m_threads_cnt;

  const gap_buffer<block_offset_type> *m_buffer;

  // For notifying threads about available buffer.
  std::mutex m_avail_mutex;
  std::condition_variable m_avail_cv;
  bool *m_avail;
  bool m_avail_no_more;

  // The mutex below is to protect m_finished. The condition
  // variable allows the caller to wait (and to be notified when done)
  // until threads complete processing their section of the buffer.
  int m_finished;
  std::mutex m_finished_mutex;
  std::condition_variable m_finished_cv;
};

template<typename block_offset_type>
void inmem_gap_updater(gap_buffer_poll<block_offset_type> *full_gap_buffers,
    gap_buffer_poll<block_offset_type> *empty_gap_buffers,
    inmem_gap_array *gap, long n_increasers) {

  gap_parallel_updater<block_offset_type> *updater =
    new gap_parallel_updater<block_offset_type>(gap, n_increasers);

  while (true) {
    // Get a buffer from the poll of full buffers.
    std::unique_lock<std::mutex> lk(full_gap_buffers->m_mutex);
    while (!full_gap_buffers->available() && !full_gap_buffers->finished())
      full_gap_buffers->m_cv.wait(lk);

    if (!full_gap_buffers->available() && full_gap_buffers->finished()) {
      // There will be no more full buffers -- exit.
      lk.unlock();
      break;
    }

    gap_buffer<block_offset_type> *b = full_gap_buffers->get();
    lk.unlock();

    // Process buffer.
    updater->update(b);

    // Add the buffer to the poll of empty buffers and notify
    // the waiting thread.
    std::unique_lock<std::mutex> lk2(empty_gap_buffers->m_mutex);
    empty_gap_buffers->add(b);
    lk2.unlock();
    empty_gap_buffers->m_cv.notify_one();
  }

  delete updater;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_UPDATE_H_INCLUDED
