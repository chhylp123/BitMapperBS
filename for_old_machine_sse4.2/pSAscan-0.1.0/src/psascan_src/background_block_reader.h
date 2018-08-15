/**
 * @file    src/psascan_src/background_block_reader.h
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

#ifndef __PSASCAN_SRC_BACKGROUND_BLOCK_READER_H_INCLUDED
#define __PSASCAN_SRC_BACKGROUND_BLOCK_READER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "utils.h"


namespace psascan_private {

struct background_block_reader {
  public:
    unsigned char *m_data;
    long m_start;
    long m_size;

  private:
    static const long k_chunk_size;

    // These variables are protected by m_mutex.
    long m_fetched;
    bool m_signal_stop;
    bool m_joined;

    std::mutex m_mutex;

    // This condition variable is used by the I/O thread to notify
    // the waiting threads when the next chunk is read.
    std::condition_variable m_cv;

    std::thread *m_thread;
    std::FILE *m_file;

  private:
    static void io_thread_main(background_block_reader &reader) {
      while (true) {
        std::unique_lock<std::mutex> lk(reader.m_mutex);
        long fetched = reader.m_fetched;
        bool signal_stop = reader.m_signal_stop;
        lk.unlock();

        if (fetched == reader.m_size || signal_stop) break;

        long toread = std::min(reader.m_size - fetched, reader.k_chunk_size);
        unsigned char *dest = reader.m_data + fetched;
        utils::read_n_objects_from_file(dest, toread, reader.m_file);

        lk.lock();
        reader.m_fetched += toread;
        lk.unlock();
        reader.m_cv.notify_all();
      }
      
      // Close the file and exit.
      std::fclose(reader.m_file);
    }

  public:
    background_block_reader(std::string filename, long start, long size) {
      m_start = start;
      m_size = size;
         
      // Initialize file and buffer.
      m_data = (unsigned char *)malloc(m_size);
      m_file = utils::open_file(filename, "r");
      std::fseek(m_file, m_start, SEEK_SET);
      m_fetched = 0;

      // Start the I/O thread.
      m_signal_stop = false;
      m_joined = false;
      m_thread = new std::thread(io_thread_main, std::ref(*this));
    }

    ~background_block_reader() {
      if (!m_joined) {
        fprintf(stderr, "\nError: the I/O thread is still not joined when "
          "destroying an object of backgroud_block_reader.\n");
        std::exit(EXIT_FAILURE);
      }
      
      // Note: m_file is already closed.
      delete m_thread;
      free(m_data);
    }

    inline void stop() {
      // Set the flag for the thread to stop.
      std::unique_lock<std::mutex> lk(m_mutex);
      m_signal_stop = true;
      lk.unlock();

      // Wait until the thread notices the flag and exits. Possibly the thread
      // is already not running, but in this case this call will do nothing.
      m_thread->join();
      
      // To detect (in the destructor) if stop() was called.
      lk.lock();
      m_joined = true;
      lk.unlock();
    }

    inline void wait(long target_fetched) {
      std::unique_lock<std::mutex> lk(m_mutex);
      while (m_fetched < target_fetched)
        m_cv.wait(lk);
      lk.unlock();
    }
};

const long background_block_reader::k_chunk_size = (1L << 20);

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_BACKGROUND_BLOCK_READER_H_INCLUDED
