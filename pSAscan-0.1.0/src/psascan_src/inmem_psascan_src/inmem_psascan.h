/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_psascan.h
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_PSASCAN_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_PSASCAN_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

#include "../bitvector.h"
#include "../multifile.h"
#include "../background_block_reader.h"
#include "inmem_gap_array.h"
#include "compute_initial_gt_bitvectors.h"
#include "initial_partial_sufsort.h"
#include "change_gt_reference_point.h"
#include "inmem_bwt_from_sa.h"
#include "inmem_compute_initial_ranks.h"
#include "parallel_merge.h"
#include "inmem_bwtsa_merge.h"
#include "pagearray.h"
#include "bwtsa.h"
#include "parallel_shrink.h"
#include "merge_schedule.h"


namespace psascan_private {
namespace inmem_psascan_private {

template<typename saidx_t, unsigned pagesize_log = 12>
void inmem_psascan(
    unsigned char *text,
    long text_length,
    unsigned char *sa_bwt,
    long max_threads = 1,
    bool compute_bwt = false,
    bool compute_gt_begin = false,
    bitvector *gt_begin = NULL,
    long max_blocks = -1,
    long text_beg = 0,
    long text_end = 0,
    long supertext_length = 0,
    std::string supertext_filename = "",
    const multifile *tail_gt_begin_reversed = NULL,
    long *i0 = NULL,
    unsigned char *tail_prefix_preread = NULL) {
  static const unsigned pagesize = (1U << pagesize_log);
  long double absolute_start = utils::wclock();
  long double start;

  if ((long)std::numeric_limits<saidx_t>::max() < text_length) {
    fprintf(stderr, "Error: text is too long (%ld bytes),\n", text_length);
    fprintf(stderr, "       std::numeric_limits<saidx_t>::max() = %ld\n",
        (long)std::numeric_limits<saidx_t>::max());
    std::exit(EXIT_FAILURE);
  }

  if (max_blocks == -1)
    max_blocks = max_threads;

  if (text_end == 0) {
    supertext_length = text_length;
    text_end = text_length;
    text_beg = 0;
    supertext_filename = "";
    tail_gt_begin_reversed = NULL;
  }

  bool has_tail = (text_end != supertext_length);

  if (!has_tail && tail_prefix_preread != NULL) {
    fprintf(stderr, "Error: has_tail == false but tail_prefix_preread != NULL\n");
    std::exit(EXIT_FAILURE);
  }

  long alignment_unit = (long)std::max(pagesize, 8U);
  long max_block_size = (text_length + max_blocks - 1) / max_blocks;
  while ((max_block_size & (alignment_unit - 1)) && max_block_size < text_length)
    ++max_block_size;

  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  if (!compute_gt_begin) {
    if (gt_begin) {
      fprintf(stderr, "Error: check gt_begin == NULL failed\n");
      std::exit(EXIT_FAILURE);
    }
    if (n_blocks > 1 || has_tail)
      gt_begin = new bitvector(text_length);
  } else {
    if (!gt_begin) {
      fprintf(stderr, "inmem_sascan: gt_begin was requested but is not allocated!\n");
      std::exit(EXIT_FAILURE);
    }
  }

  fprintf(stderr, "Text length = %ld (%.2LfMiB)\n", text_length, text_length / (1024.L * 1024));
  fprintf(stderr, "Max block size = %ld (%.2LfMiB)\n", max_block_size, max_block_size / (1024.L * 1024));
  fprintf(stderr, "Max blocks = %ld\n", max_blocks);
  fprintf(stderr, "Number of blocks = %ld\n", n_blocks);
  fprintf(stderr, "Max threads = %ld\n", max_threads);
  fprintf(stderr, "sizeof(saidx_t) = %lu\n", sizeof(saidx_t));
  fprintf(stderr, "Pagesize = %u\n", (1U << pagesize_log));
  fprintf(stderr, "Compute bwt = %s\n", compute_bwt ? "true" : "false");
  fprintf(stderr, "Compute gt begin = %s\n", compute_gt_begin ? "true" : "false");
  fprintf(stderr, "Text beg = %ld\n", text_beg);
  fprintf(stderr, "Text end = %ld\n", text_end);
  fprintf(stderr, "Supertext length = %ld (%.2LfMiB)\n", supertext_length, supertext_length / (1024.L * 1024));
  fprintf(stderr, "Supertext filename = %s\n", supertext_filename.c_str());
  fprintf(stderr, "Has tail = %s\n", has_tail ? "true" : "false");
  fprintf(stderr, "\n");

  bwtsa_t<saidx_t> *bwtsa = (bwtsa_t<saidx_t> *)sa_bwt;

  // Initialize reading of the tail prefix in the background.
  long tail_length = supertext_length - text_end;
  long tail_prefix_length = std::min(text_length, tail_length);

  background_block_reader *tail_prefix_background_reader = NULL;
  if (has_tail && tail_prefix_preread == NULL)
    tail_prefix_background_reader =
      new background_block_reader(supertext_filename, text_end, tail_prefix_length);

  //----------------------------------------------------------------------------
  // STEP 1: compute initial bitvectors, and partial suffix arrays.
  //----------------------------------------------------------------------------
  if (n_blocks > 1 || compute_gt_begin || has_tail) {
    fprintf(stderr, "Compute initial bitvectors:\n");
    start = utils::wclock();
    compute_initial_gt_bitvectors(text, text_length, gt_begin, max_block_size,
        max_threads, text_end, supertext_length, tail_gt_begin_reversed,
        tail_prefix_background_reader, tail_prefix_preread);
    fprintf(stderr, "Time: %.2Lf\n\n", utils::wclock() - start);
  }

  fprintf(stderr, "Initial sufsort:\n");
  start = utils::wclock();
  initial_partial_sufsort(text, text_length, gt_begin, bwtsa, max_block_size, max_threads, has_tail);
  fprintf(stderr, "Time: %.2Lf\n", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 2: compute matrix of block ranks.
  //----------------------------------------------------------------------------
  fprintf(stderr, "Compute matrix of initial ranks: ");
  start = utils::wclock();
  long **block_rank_matrix = new long*[n_blocks];
  for (long j = 0; j < n_blocks; ++j)
    block_rank_matrix[j] = new long[n_blocks];
  compute_block_rank_matrix<saidx_t>(text, text_length, bwtsa,
      max_block_size, text_beg, supertext_length, supertext_filename,
      tail_gt_begin_reversed, tail_prefix_background_reader,
      tail_prefix_preread, block_rank_matrix);

  // Stop reading next block in the background or free memory taken by next block.
  if (has_tail) {
    if (tail_prefix_background_reader != NULL) {
      tail_prefix_background_reader->stop();
      delete tail_prefix_background_reader;
    } else free(tail_prefix_preread);
  }

  fprintf(stderr, "%.2Lf\n\n", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 3: compute the gt bitvectors for blocks that will be on the right
  //         side during the merging.
  //----------------------------------------------------------------------------
  if (n_blocks > 1 || compute_gt_begin) {
    fprintf(stderr, "Overwriting gt_end with gt_begin: ");
    start = utils::wclock();
    gt_end_to_gt_begin(text, text_length, gt_begin, max_block_size);
    fprintf(stderr, "%.2Lf\n\n", utils::wclock() - start);
  }

  float rl_ratio = 10.L;  // estimated empirically
  long max_ram_usage_per_input_byte = 10L;  // peak ram usage = 10n
  int max_left_size = std::max(1, (int)floor(n_blocks * (((long double)max_ram_usage_per_input_byte - (2.125L + sizeof(saidx_t))) / 5.L)));
  fprintf(stderr, "Assumed rl_ratio: %.2f\n", rl_ratio);
  fprintf(stderr, "Max left size = %d\n", max_left_size);
  fprintf(stderr, "Peak memory usage during last merging = %.3Lfn\n",
      (2.125L + sizeof(saidx_t)) + (5.L * max_left_size) / n_blocks);
  MergeSchedule schedule(n_blocks, rl_ratio, max_left_size);

  fprintf(stderr, "Skewed merge schedule:\n");
  print_schedule(schedule, n_blocks);
  fprintf(stderr, "\n");

  long *i0_array = new long[n_blocks];
  if (n_blocks > 1 || compute_bwt) {
    for (long block_id = 0; block_id < n_blocks; ++block_id) {
      long block_end = text_length - (n_blocks - 1 - block_id) * max_block_size;
      long block_beg = std::max(0L, block_end - max_block_size);
      long block_size = block_end - block_beg;

      if (block_id + 1 != n_blocks || compute_bwt) {
        fprintf(stderr, "Computing BWT for block %ld: ", block_id + 1);
        long double bwt_start = utils::wclock();
        compute_bwt_in_bwtsa<saidx_t>(text + block_beg, block_size,
            bwtsa + block_beg, max_threads, i0_array[block_id]);
        fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
      }
    }
    fprintf(stderr, "\n");
  }

  if (n_blocks > 1) {
    long i0_result;
    pagearray<bwtsa_t<saidx_t>, pagesize_log> *result =
      inmem_bwtsa_merge<saidx_t, pagesize_log>(text, text_length, bwtsa,
          gt_begin, max_block_size, 0, n_blocks, max_threads, compute_gt_begin,
          compute_bwt, i0_result, schedule, text_beg, text_end,
          supertext_length, supertext_filename, tail_gt_begin_reversed,
          i0_array, block_rank_matrix);
    if (i0) *i0 = i0_result;

    // Permute SA to plain array.
    fprintf(stderr, "\nPermuting the resulting SA to plain array: ");
    start = utils::wclock();
    result->permute_to_plain_array(max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

    delete result;
  } else if (compute_bwt) {
    if (i0) *i0 = i0_array[0];
  }
  delete[] i0_array;
  for (long j = 0; j < n_blocks; ++j)
    delete[] block_rank_matrix[j];
  delete[] block_rank_matrix;

  if (!compute_gt_begin && (n_blocks > 1 || has_tail)) {
    delete gt_begin;
    gt_begin = NULL;
  }

  unsigned char *bwt = NULL;
  if (compute_bwt) {
    // Allocate aux, copy bwt into aux.
    fprintf(stderr, "Copying bwtsa.bwt into aux memory: ");
    start = utils::wclock();
    bwt = (unsigned char *)malloc(text_length);
    parallel_copy<bwtsa_t<saidx_t>, unsigned char>(bwtsa, bwt, text_length, max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  fprintf(stderr, "Shrinking bwtsa.sa into sa: ");
  start = utils::wclock();

  parallel_shrink<bwtsa_t<saidx_t>, saidx_t>(bwtsa, text_length, max_threads);

  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  if (compute_bwt) {
    // Copy from aux into the end of bwtsa.
    fprintf(stderr, "Copying bwt from aux memory to the end of bwtsa: ");
    start = utils::wclock();
    unsigned char *dest = (unsigned char *)(((saidx_t *)bwtsa) + text_length);
    parallel_copy<unsigned char, unsigned char>(bwt, dest, text_length, max_threads);
    free(bwt);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  long double total_sascan_time = utils::wclock() - absolute_start;
  fprintf(stderr, "\nTotal time:\n");
  fprintf(stderr, "\tabsolute: %.2Lf\n", total_sascan_time);
  fprintf(stderr, "\trelative: %.4Lfs/MiB\n", total_sascan_time / ((long double)text_length / (1 << 20)));
  fprintf(stderr, "Speed: %.2LfMiB/s\n", ((long double)text_length / (1 << 20)) / total_sascan_time);
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_PSASCAN_H_INCLUDED
