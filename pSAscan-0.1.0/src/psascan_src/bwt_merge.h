/**
 * @file    src/psascan_src/bwt_merge.h
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

#ifndef __PSASCAN_SRC_BWT_MERGE_H_INCLUDED
#define __PSASCAN_SRC_BWT_MERGE_H_INCLUDED

#include <thread>
#include <algorithm>

#include "bitvector.h"
#include "ranksel_support.h"


namespace psascan_private {

//==============================================================================
// Compute bwt[beg..end).
//==============================================================================
void merge_bwt_aux(long beg, long end, long left_ptr, long right_ptr,
    const unsigned char *left_bwt, const unsigned char *right_bwt, unsigned char *bwt,
    const bitvector *bv) {
  for (long i = beg; i < end; ++i) {
    if (bv->get(i)) bwt[i] = right_bwt[right_ptr++];
    else bwt[i] = left_bwt[left_ptr++];
  }
}

void compute_initial_rank(long i, const ranksel_support *ranksel, long &result) {
  result = ranksel->rank(i);
}

//==============================================================================
// Merge partial bwt of half-blocks (of size left_size and right_size) into
// partial bwt of the whole block.
//==============================================================================
long merge_bwt(const unsigned char *left_bwt, const unsigned char *right_bwt,
    long left_size, long right_size, long left_block_i0, long right_block_i0,
    unsigned char left_block_last, unsigned char *bwt, const bitvector *bv,
    long max_threads) {
  long block_size = left_size + right_size;

  // 1
  //
  // Initialize rank/select queries support for bv.
  ranksel_support *bv_ranksel = new ranksel_support(bv, block_size, max_threads);

  // 2
  //
  // Compute range size.
  long max_range_size = (block_size + max_threads - 1) / max_threads;
  long n_ranges = (block_size + max_range_size - 1) / max_range_size;

  // 3
  //
  // Compute starting parameters for each thread.
  long *left_ptr = new long[n_ranges];
  long *right_ptr = new long[n_ranges];
  long *rank_at_range_beg = new long[n_ranges];

  std::thread **threads = new std::thread*[n_ranges];
  for (long t = 0; t < n_ranges; ++t) {
    long range_beg = t * max_range_size;
    threads[t] = new std::thread(compute_initial_rank,
        range_beg, bv_ranksel, std::ref(rank_at_range_beg[t]));
  }

  for (long t = 0; t < n_ranges; ++t) threads[t]->join();
  for (long t = 0; t < n_ranges; ++t) delete threads[t];

  for (long t = 0; t < n_ranges; ++t) {
    long range_beg = t * max_range_size;
    left_ptr[t] = range_beg - rank_at_range_beg[t];
    right_ptr[t] = rank_at_range_beg[t];
  }
  delete[] rank_at_range_beg;

  // 4
  //
  // Merge BWTs in parallel.
  for (long t = 0; t < n_ranges; ++t) {
    long range_beg = max_range_size * t;
    long range_end = std::min(range_beg + max_range_size, block_size);

    threads[t] = new std::thread(merge_bwt_aux, range_beg, range_end,
        left_ptr[t], right_ptr[t], left_bwt, right_bwt, bwt, bv);
  }

  for (long t = 0; t < n_ranges; ++t) threads[t]->join();
  for (long t = 0; t < n_ranges; ++t) delete threads[t];
  delete[] threads;
  delete[] left_ptr;
  delete[] right_ptr;

  // 5
  //
  // Find position j = select_1(bv, right_block_i0) and replace bwt[j] with
  // left_block_last. To speed up the search for j, we use sparse_rank.
  bwt[bv_ranksel->select1(right_block_i0)] = left_block_last;

  // 6
  //
  // Compute the returned value.
  long block_i0 = bv_ranksel->select0(left_block_i0);

  // 7
  //
  // Clean up and exit.
  delete bv_ranksel;
  return block_i0;
}

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_BWT_MERGE_H_INCLUDED
