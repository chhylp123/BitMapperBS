/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_bwtsa_merge.h
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWTSA_MERGE_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWTSA_MERGE_H_INCLUDED

#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

#include "../bitvector.h"
#include "../multifile.h"
#include "inmem_gap_array.h"
#include "inmem_compute_gap.h"
#include "parallel_merge.h"
#include "pagearray.h"
#include "bwtsa.h"
#include "merge_schedule.h"


namespace psascan_private {
namespace inmem_psascan_private {

template<typename saidx_t, unsigned pagesize_log>
pagearray<bwtsa_t<saidx_t>, pagesize_log> *inmem_bwtsa_merge(
    const unsigned char *text,
    long text_length,
    bwtsa_t<saidx_t> *bwtsa,
    bitvector *gt,
    long max_block_size,
    long range_beg,
    long range_end,
    long max_threads,
    bool need_gt,
    bool need_bwt,
    long &result_i0,
    MergeSchedule &schedule,
    long text_beg,
    long text_end,
    long supertext_length,
    std::string supertext_filename,
    const multifile *tail_gt_begin_reversed,
    long *i0_array,
    long **block_rank_matrix) {
  typedef pagearray<bwtsa_t<saidx_t>, pagesize_log> pagearray_type;

  long shift = (max_block_size - text_length % max_block_size) % max_block_size;
  long range_size = range_end - range_beg;

  if (range_size == 1) {
    long block_beg = range_beg * max_block_size;
    long block_end = block_beg + max_block_size;
    block_beg = std::max(0L, block_beg - shift);
    block_end -= shift;

    result_i0 = i0_array[range_beg];
    pagearray_type *bwtsa_pagearray =
      new pagearray_type(bwtsa + block_beg, bwtsa + block_end);
    return bwtsa_pagearray;
  }

  //----------------------------------------------------------------------------
  // STEP 1: Split the blocks in the left and right group.
  //----------------------------------------------------------------------------
  long lrange_size = schedule.left_size(range_size);
  long rrange_size = range_size - lrange_size;

  long lrange_beg = range_beg;
  long lrange_end = range_beg + lrange_size;
  long rrange_beg = lrange_end;
  long rrange_end = rrange_beg + rrange_size;

  long lbeg = lrange_beg * max_block_size;
  long rbeg = rrange_beg * max_block_size;
  long lend = rbeg;
  long rend = rbeg + rrange_size * max_block_size;
  lbeg = std::max(0L, lbeg - shift);
  rbeg -= shift;
  lend -= shift;
  rend -= shift;

  long lsize = lend - lbeg;
  long rsize = rend - rbeg;

  //----------------------------------------------------------------------------
  // STEP 2: Compute partial SAs and BWTs for left and right block.
  //----------------------------------------------------------------------------

  // 2.a
  //
  // Left block
  long left_i0;
  pagearray_type *l_bwtsa = inmem_bwtsa_merge<saidx_t, pagesize_log>(text,
      text_length, bwtsa, gt, max_block_size, lrange_beg, lrange_end,
      max_threads, need_gt, true, left_i0, schedule, text_beg, text_end,
      supertext_length, supertext_filename, tail_gt_begin_reversed, i0_array,
      block_rank_matrix);

  // 2.b
  // 
  // Right block
  long right_i0;
  pagearray_type *r_bwtsa = inmem_bwtsa_merge<saidx_t, pagesize_log>(text,
      text_length, bwtsa, gt, max_block_size, rrange_beg, rrange_end,
      max_threads, true, need_bwt, right_i0, schedule, text_beg, text_end,
      supertext_length, supertext_filename, tail_gt_begin_reversed, i0_array,
      block_rank_matrix);

  //----------------------------------------------------------------------------
  // STEP 3: Merge partial SAs and BWTs.
  //----------------------------------------------------------------------------
  fprintf(stderr, "Merging blocks %ld-%ld with %ld-%ld\n",
      lrange_beg + 1, lrange_end, rrange_beg + 1, rrange_end);
  long double start = utils::wclock();

  // 3.a
  //
  // Compute gap
  fprintf(stderr, "  Computing gap:\n");
  inmem_gap_array *gap;
  long double rank_init_time;
  long double streaming_time;
  long double start1 = utils::wclock();
  inmem_compute_gap<saidx_t, pagesize_log>(text, text_length, lbeg, lsize,
      rsize, *l_bwtsa, gt, gap, max_threads, need_gt, left_i0, (1L << 21),
      rank_init_time, streaming_time, block_rank_matrix, lrange_beg,
      lrange_size, rrange_size);
  fprintf(stderr, "  Time: %.2Lf\n", utils::wclock() - start1);

  // 3.b
  //
  // Merge partial SAs and BWTs
  fprintf(stderr, "  Merging SA/BWT:  ");
  start1 = utils::wclock();
  long delta_i0;
  if (need_bwt)
    (*r_bwtsa)[right_i0].bwt = text[rbeg - 1];
  pagearray_type *result = parallel_merge(l_bwtsa, r_bwtsa, gap,
      max_threads, left_i0, delta_i0, lsize);
  result_i0 = left_i0 + delta_i0;
  long double merging_time = utils::wclock() - start1;
  fprintf(stderr, "total: %.2Lf\n", merging_time);

  // 3.c
  //
  // Clean up.
  start1 = utils::wclock();
  delete l_bwtsa;
  delete r_bwtsa;
  delete gap;
  long double cleaning_time = utils::wclock() - start1;
  if (cleaning_time > 0.2L)
    fprintf(stderr, "Cleaning: %.2Lf\n", cleaning_time);

  long double time_per_elem_left = merging_time / (lsize + rsize) + rank_init_time / lsize;
  long double time_per_elem_right = merging_time / (lsize + rsize) + streaming_time / rsize;
  long double ratio = time_per_elem_right / time_per_elem_left;
  fprintf(stderr, "Time: %.2Lf (rl_ratio = %.3Lf)\n",
      utils::wclock() - start, ratio);

  return result;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWTSA_MERGE_H_INCLUDED
