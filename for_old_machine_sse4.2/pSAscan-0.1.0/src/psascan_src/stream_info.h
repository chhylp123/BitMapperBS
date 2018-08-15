/**
 * @file    src/psascan_src/stream_info.h
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

#ifndef __PSASCAN_SRC_STREAM_INFO_H_INCLUDED
#define __PSASCAN_SRC_STREAM_INFO_H_INCLUDED

#include <mutex>
#include <algorithm>

#include "utils.h"


namespace psascan_private {

//=============================================================================
// Used to store progress information for different threads during streaming.
//=============================================================================
struct stream_info {
  stream_info(long thread_count, long tostream)
    : m_update_count(0L),
      m_thread_count(thread_count),
      m_tostream(tostream) {
    m_streamed = new long[thread_count];
    std::fill(m_streamed, m_streamed + thread_count, 0L);

    m_idle_update = new long double[thread_count];
    m_idle_work  = new long double[thread_count];
    std::fill(m_idle_update, m_idle_update + thread_count, 0.L);
    std::fill(m_idle_work, m_idle_work + thread_count, 0.L);

    m_timestamp = utils::wclock();
  }

  ~stream_info() {
    delete[] m_streamed;
    delete[] m_idle_work;
    delete[] m_idle_update;
  }

  long m_update_count;     // number of updates
  long m_thread_count;     // number of threads
  long m_tostream;         // total text length to stream
  long double m_timestamp; // when the streaming started
  long *m_streamed;        // how many bytes streamed by each thread
  long double *m_idle_update;
  long double *m_idle_work;

  std::mutex m_mutex;
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_STREAM_INFO_H_INCLUDED
