/**
 * @file    src/psascan_src/async_multifile_bit_stream_reader.h
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

#ifndef __PSASCAN_SRC_ASYNC_MULTIFILE_BIT_STREAM_READER_H_INCLUDED
#define __PSASCAN_SRC_ASYNC_MULTIFILE_BIT_STREAM_READER_H_INCLUDED

#include <thread>
#include <mutex>
#include <vector>
#include <algorithm>
#include <condition_variable>

#include "utils.h"
#include "multifile.h"


namespace psascan_private {

struct async_multifile_bit_stream_reader {
  async_multifile_bit_stream_reader(const multifile *m, long start_pos = 0L,
      long bufsize = (4L << 20)) {
    m_files_info = m->files_info;
  
    long items = std::max(2L, bufsize);
    m_buf_size = items / 2L;

    // Reset counters.
    m_active_buf_filled = 0;
    m_passive_buf_filled = 0;
    m_active_buf_pos = 0;

    // Initialize buffers.
    m_active_buf = (unsigned char *)malloc(m_buf_size);
    m_passive_buf = (unsigned char *)malloc(m_buf_size);

    // Initialize the reading.
    init(start_pos);
  }

  void init(long start_pos) {
    m_total_read_buf = start_pos;

    m_file = NULL;
    for (size_t j = 0; j < m_files_info.size(); ++j) {
      if (m_files_info[j].m_beg <= start_pos && start_pos < m_files_info[j].m_end) {
        m_file_id = j;
        m_file = utils::open_file(m_files_info[j].m_filename, "r");
        break;
      }
    }

    if (m_file != NULL) {
      long offset = start_pos - m_files_info[m_file_id].m_beg;
      std::fseek(m_file, offset >> 3, SEEK_SET);

      m_cur_byte = 0;
      m_cur_bit = (offset & 7L);
      m_active_buf_pos = m_cur_bit;
      m_total_read_buf -= m_cur_bit;

      long file_left = m_files_info[m_file_id].m_end - m_total_read_buf;
      m_active_buf_filled = std::min(file_left, 8L * m_buf_size);
      long toread_bytes = (m_active_buf_filled + 7L) / 8L;
      utils::read_n_objects_from_file(m_active_buf, toread_bytes, m_file);
      m_total_read_buf += m_active_buf_filled;
      if (m_total_read_buf == m_files_info[m_file_id].m_end) {
        std::fclose(m_file);
        m_file = NULL;
      }
    }

    m_avail = true;
    m_finished = false;
    m_thread = new std::thread(async_io_code, this);
  }

  inline bool read() {
    if (m_active_buf_pos == m_active_buf_filled)
      receive_new_buffer();

    bool result = (m_active_buf[m_cur_byte] & (1 << m_cur_bit));
    ++m_cur_bit;
    ++m_active_buf_pos;
    if (m_cur_bit == 8) {
      m_cur_bit = 0;
      ++m_cur_byte;
    }

    return result;
  }

  ~async_multifile_bit_stream_reader() {
    // Let the I/O thread know that we are done.
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
    if (m_file)
      std::fclose(m_file);
  }

  static void async_io_code(async_multifile_bit_stream_reader *file) {
    while (true) {
      // Wait until the passive buffer is available.
      std::unique_lock<std::mutex> lk(file->m_mutex);
      while (!(file->m_avail) && !(file->m_finished))
        file->m_cv.wait(lk);

      if (!(file->m_avail) && (file->m_finished)) {
        // We're done, terminate the thread.
        lk.unlock();
        return;
      }
      lk.unlock();
      
      if (file->m_file == NULL) {
        // Find the next file to open.
        for (size_t j = 0; j < file->m_files_info.size(); ++j)
          if (file->m_files_info[j].m_beg == file->m_total_read_buf) {
            file->m_file_id = j;
            file->m_file = utils::open_file(file->m_files_info[j].m_filename, "r");
            break;
          }
      }

      // If file ID was found, we perform the read.
      // Otherwise there is no more data to prefetch.
      if (file->m_file != NULL) {
        long file_left = file->m_files_info[file->m_file_id].m_end - file->m_total_read_buf;
        file->m_passive_buf_filled = std::min(file_left, 8L * (file->m_buf_size));
        long toread_bytes = (file->m_passive_buf_filled + 7L) / 8L;
        utils::read_n_objects_from_file(file->m_passive_buf, toread_bytes, file->m_file);
        file->m_total_read_buf += file->m_passive_buf_filled;
        if (file->m_total_read_buf == file->m_files_info[file->m_file_id].m_end) {
          std::fclose(file->m_file);
          file->m_file = NULL;
        }
      }

      // Let the caller know that the I/O thread finished reading.
      lk.lock();
      file->m_avail = false;
      lk.unlock();
      file->m_cv.notify_one();
    }
  }

  void receive_new_buffer() {
    // Wait until the I/O thread finishes reading the previous
    // buffer. Most of the time this step is instantaneous.
    std::unique_lock<std::mutex> lk(m_mutex);
    while (m_avail == true)
      m_cv.wait(lk);

    // Set the new active buffer.
    std::swap(m_active_buf, m_passive_buf);
    m_active_buf_filled = m_passive_buf_filled;
    m_active_buf_pos = 0;
    m_cur_byte = 0;
    m_cur_bit = 0;

    // Let the I/O thread know that it can now
    // prefetch another buffer.
    m_avail = true;
    lk.unlock();
    m_cv.notify_one();
  }

private:
  std::FILE *m_file;       // file handler  
  long m_total_read_buf;   // total number of items read from files into buffers
  long m_file_id;
  std::vector<single_file_info> m_files_info;

  // Buffers used for asynchronous reading.
  unsigned char *m_active_buf;
  unsigned char *m_passive_buf;
  long m_buf_size;
  long m_active_buf_pos;
  long m_active_buf_filled;
  long m_passive_buf_filled;

  long m_cur_byte;
  long m_cur_bit;

  // For synchronization with thread doing asynchronous reading.
  std::thread *m_thread;
  std::mutex m_mutex;
  std::condition_variable m_cv;
  bool m_finished;
  bool m_avail;
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_ASYNC_MULTIFILE_BIT_STREAM_READER_H_INCLUDED
