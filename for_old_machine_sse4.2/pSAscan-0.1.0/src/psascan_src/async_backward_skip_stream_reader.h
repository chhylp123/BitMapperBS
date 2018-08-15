/**
 * @file    src/psascan_src/async_backward_skip_stream_reader.h
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

#ifndef __PSASCAN_SRC_ASYNC_BACKWARD_SKIP_STREAM_READER_H_INCLUDED
#define __PSASCAN_SRC_ASYNC_BACKWARD_SKIP_STREAM_READER_H_INCLUDED

#include <cstdio>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "utils.h"


namespace psascan_private {

template<typename value_type>
struct async_backward_skip_stream_reader {
  template<typename T>
  static void io_thread_code(async_backward_skip_stream_reader<T> *reader) {
    while (true) {
      // Wait until the passive buffer is available.
      std::unique_lock<std::mutex> lk(reader->m_mutex);
      while (!(reader->m_avail) && !(reader->m_finished))
        reader->m_cv.wait(lk);

      if (!(reader->m_avail) && (reader->m_finished)) {
        // We're done, terminate the thread.
        lk.unlock();
        return;
      }
      lk.unlock();

      // Safely read the data from disk.
      long filepos = std::ftell(reader->m_file) / sizeof(T);
      long toread = std::min(reader->m_buf_size, filepos - reader->m_active_buf_filled);
      if (toread > 0) {
        std::fseek(reader->m_file, -((reader->m_active_buf_filled + toread) * sizeof(T)), SEEK_CUR);
        reader->m_passive_buf_filled = std::fread(reader->m_passive_buf, sizeof(T), toread, reader->m_file);
      }

      // Let the caller know that the I/O thread finished reading.
      lk.lock();
      reader->m_avail = false;
      lk.unlock();
      reader->m_cv.notify_one();
    }
  }

  async_backward_skip_stream_reader(std::string filename, long skip_elems, long bufsize = (4 << 20)) {
    m_file = utils::open_file(filename.c_str(), "r");
    std::fseek(m_file, -(skip_elems * sizeof(value_type)), SEEK_END);


    // Initialize buffers.
    long elems = std::max(2UL, (bufsize + sizeof(value_type) - 1) / sizeof(value_type));
    m_buf_size = elems / 2;

    m_active_buf_filled = 0L;
    m_passive_buf_filled = 0L;
    m_active_buf_pos = -1L;
    m_active_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));
    m_passive_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));

    m_finished = false;
    
    // Start the I/O thread and immediately start reading.
    m_avail = true;
    m_thread = new std::thread(io_thread_code<value_type>, this);
  }
  
  ~async_backward_skip_stream_reader() {
    // Let the I/O thread know that we're done.
    std::unique_lock<std::mutex> lk(m_mutex);
    m_finished = true;
    lk.unlock();
    m_cv.notify_one();

    // Wait for the thread to finish.
    m_thread->join();
    
    // Clean up.
    delete m_thread;
    free(m_active_buf);
    free(m_passive_buf);
    std::fclose(m_file);
  }

  // This function checks if the reading thread has already
  // prefetched the next buffer (the request should have been
  // issued before), and waits in case the prefetching was not
  // completed yet.
  void receive_new_buffer() {
    // Wait until the I/O thread finishes reading the previous
    // buffer. In most cases this step is instantaneous.
    std::unique_lock<std::mutex> lk(m_mutex);
    while (m_avail == true)
      m_cv.wait(lk);

    // Set the new active buffer.
    std::swap(m_active_buf, m_passive_buf);
    m_active_buf_filled = m_passive_buf_filled;
    m_active_buf_pos = m_active_buf_filled - 1L;

    // Let the I/O thread know that it can now prefetch
    // another buffer.
    m_avail = true;
    lk.unlock();
    m_cv.notify_one();
  }

  inline value_type read() {
    if (m_active_buf_pos < 0L) {
      // The active buffer run out of data.
      // At this point we need to swap it with the passive
      // buffer. The request to read that passive buffer should
      // have been scheduled long time ago, so hopefully the
      // buffer is now available. We check for that, but we
      // also might wait, if the reading has not yet been
      // finished. At this point we also already schedule
      // the next read.
      receive_new_buffer();
    }

    return m_active_buf[m_active_buf_pos--];
  }

private:
  value_type *m_active_buf;
  value_type *m_passive_buf;

  long m_buf_size;
  long m_active_buf_pos;
  long m_active_buf_filled;
  long m_passive_buf_filled;

  // Used for synchronization with the I/O thread.
  std::mutex m_mutex;
  std::condition_variable m_cv;
  bool m_avail;
  bool m_finished;

  std::FILE *m_file;
  std::thread *m_thread;
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_ASYNC_BACKWARD_SKIP_STREAM_READER_H_INCLUDED
