/**
 * @file    src/psascan_src/async_stream_writer.h
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

#ifndef __PSASCAN_SRC_ASYNC_STREAM_WRITER_H_INCLUDED
#define __PSASCAN_SRC_ASYNC_STREAM_WRITER_H_INCLUDED

#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "utils.h"


namespace psascan_private {

template<typename value_type>
struct async_stream_writer {
  template<typename T>
  static void io_thread_code(async_stream_writer<T> *writer) {
    while (true) {
      // Wait until the passive buffer is available.
      std::unique_lock<std::mutex> lk(writer->m_mutex);
      while (!(writer->m_avail) && !(writer->m_finished))
        writer->m_cv.wait(lk);

      if (!(writer->m_avail) && (writer->m_finished)) {
        // We're done, terminate the thread.
        lk.unlock();
        return;
      }
      lk.unlock();

      // Safely write the data to disk.
      utils::add_objects_to_file(writer->m_passive_buf,
          writer->m_passive_buf_filled, writer->m_file);

      // Let the caller know that the I/O thread finished writing.
      lk.lock();
      writer->m_avail = false;
      lk.unlock();
      writer->m_cv.notify_one();
    }
  }

  async_stream_writer(std::string filename, long bufsize = (4 << 20)) {
    m_file = utils::open_file(filename.c_str(), "w");

    // Initialize buffers.    
    long elems = std::max(2UL, (bufsize + sizeof(value_type) - 1) / sizeof(value_type));
    m_buf_size = elems / 2;  // both buffers are of the same size

    m_active_buf_filled = 0L;
    m_passive_buf_filled = 0L;

    m_active_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));
    m_passive_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));

    m_avail = false;
    m_finished = false;
    
    // Start the I/O thread.
    m_thread = new std::thread(io_thread_code<value_type>, this);
  }
  
  ~async_stream_writer() {
    // Write the partially filled active buffer to disk.
    if (m_active_buf_filled > 0L)
      send_active_buf_to_write();

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

  // Passes on the active buffer (full, unless it's the last one,
  // partially filled, buffer passed from destructor) to the I/O thread.
  void send_active_buf_to_write() {
    // Wait until the I/O thread finishes writing the previous buffer.
    std::unique_lock<std::mutex> lk(m_mutex);
    while (m_avail == true)
      m_cv.wait(lk);
      
    // Set the new passive buffer.
    std::swap(m_active_buf, m_passive_buf);
    m_passive_buf_filled = m_active_buf_filled;
    m_active_buf_filled = 0L;

    // Let the I/O thread know that the buffer is waiting.
    m_avail = true;
    lk.unlock();
    m_cv.notify_one();
  }
  
  inline void write(value_type x) {
    m_active_buf[m_active_buf_filled++] = x;
    
    // If the active buffer was full, send it to I/O thread.
    // This function may wait a bit until the I/O thread
    // finishes writing the previous passive buffer.
    if (m_active_buf_filled == m_buf_size)
      send_active_buf_to_write();
  }

private:
  value_type *m_active_buf;
  value_type *m_passive_buf;
  
  long m_buf_size;  // size of each of the buffers
  long m_active_buf_filled;
  long m_passive_buf_filled;

  // Used for synchronization with the I/O thread.
  bool m_avail;     // signals availability of buffer for I/O thread
  bool m_finished;  // signals the end of writing
  std::mutex m_mutex;
  std::condition_variable m_cv;

  std::FILE *m_file;
  std::thread *m_thread;
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_ASYNC_STREAM_WRITER_H_INCLUDED
