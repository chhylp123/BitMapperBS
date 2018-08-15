/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_bwt_from_sa.h
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_H_INCLUDED

#include <algorithm>
#include <thread>

#include "../utils.h"
#include "bwtsa.h"


namespace psascan_private {
namespace inmem_psascan_private {

template<typename saidx_t>
void compute_bwt_in_bwtsa_aux(const unsigned char *text, long beg,
    long end, bwtsa_t<saidx_t> *dest, long *i0) {
  *i0 = -1;
  for (long j = beg; j < end; ++j) {
    if (dest[j].sa) dest[j].bwt = text[dest[j].sa - 1];
    else { dest[j].bwt = 0; *i0 = j; }
  }
}

template<typename saidx_t>
void compute_bwt_in_bwtsa(const unsigned char *text, long length,
  bwtsa_t<saidx_t> *dest, long max_threads, long &result) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;
  long *index_0 = new long[n_blocks];

  // Compute bwt and find i0, where sa[i0] == 0.
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);

    threads[i] = new std::thread(compute_bwt_in_bwtsa_aux<saidx_t>,
        text, block_beg, block_end, dest, index_0 + i);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;

  // Find and return i0.
  result = -1;
  for (long i = 0; i < n_blocks; ++i)
    if (index_0[i] != -1) result = index_0[i];
  delete[] index_0;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_INCLUDED
