/**
 * @file    src/psascan_src/async_vbyte_stream_reader.h
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

#ifndef __PSASCAN_SRC_ASYNC_VBYTE_STREAM_READER_H_INCLUDED
#define __PSASCAN_SRC_ASYNC_VBYTE_STREAM_READER_H_INCLUDED

#include <cstdio>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "utils.h"


namespace psascan_private {

template<typename value_type>
struct async_vbyte_stream_reader {
  static void io_thread_code(async_vbyte_stream_reader *reader) {
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
      long count = std::fread(reader->m_passive_buf, 1, reader->m_buf_size + 128, reader->m_file);
      if (count > reader->m_buf_size) {
        reader->m_passive_buf_filled = reader->m_buf_size;
        std::fseek(reader->m_file, reader->m_buf_size - count, SEEK_CUR);
      } else reader->m_passive_buf_filled = count;
 
      // Let the caller know that the I/O thread finished reading.
      lk.lock();
      reader->m_avail = false;
      lk.unlock();
      reader->m_cv.notify_one();
    }
  }

  async_vbyte_stream_reader(std::string filename, long bufsize = (4L << 20)) {
    m_file = utils::open_file(filename.c_str(), "r");

    // Initialize buffers.
    long elems = std::max(4096L, bufsize);
    m_buf_size = elems / 2;

    m_active_buf_filled = 0L;
    m_passive_buf_filled = 0L;
    m_active_buf_pos = 0L;
    m_active_buf = (unsigned char *)malloc(m_buf_size + 128);
    m_passive_buf = (unsigned char *)malloc(m_buf_size + 128);

    m_finished = false;
    
    // Start the I/O thread and immediately start reading.
    m_avail = true;
    m_thread = new std::thread(io_thread_code, this);
  }

  ~async_vbyte_stream_reader() {
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
  void receive_new_buffer(long skipped_bytes) {
    // Wait until the I/O thread finishes reading the previous
    // buffer. In most cases, this step is instantaneous.
    std::unique_lock<std::mutex> lk(m_mutex);
    while (m_avail == true)
      m_cv.wait(lk);

    // Set the new active buffer.
    std::swap(m_active_buf, m_passive_buf);
    m_active_buf_filled = m_passive_buf_filled;
    m_active_buf_pos = skipped_bytes;

    // Let the I/O thread know that it can now prefetch
    // another buffer.
    m_avail = true;
    lk.unlock();
    m_cv.notify_one();
  }

  inline value_type read() {
    if (m_active_buf_pos >= m_active_buf_filled) {
      // The active buffer run out of data.
      // At this point we need to swap it with the passive
      // buffer. The request to read that passive buffer should
      // have been scheduled long time ago, so hopefully the
      // buffer is now available. We check for that, but we
      // also might wait, if the reading has not yet been finished.
      // At this point we also already schedule the next read.
      receive_new_buffer(m_active_buf_pos - m_active_buf_filled);
    }

    value_type result = 0L;
    long offset = 0L;
    while (m_active_buf[m_active_buf_pos] & 0x80) {
      result |= (((value_type)m_active_buf[m_active_buf_pos++] & 0x7F) << offset);
      offset += 7;
    }
    result |= ((value_type)m_active_buf[m_active_buf_pos++] << offset);

    return result;
  }

private:
  unsigned char *m_active_buf;
  unsigned char *m_passive_buf;

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

#endif  // __PSASCAN_SRC_ASYNC_VBYTE_STREAM_READER_H_INCLUDED
