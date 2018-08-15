/**
 * @file    src/psascan_src/compute_left_gap.h
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

#ifndef __PSASCAN_SRC_COMPUTE_LEFT_GAP_H_INCLUDED
#define __PSASCAN_SRC_COMPUTE_LEFT_GAP_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "bitvector.h"
#include "ranksel_support.h"
#include "gap_array.h"
#include "parallel_utils.h"


namespace psascan_private {

//==============================================================================
// Compute the range_gap values corresponging to bv[part_beg..part_end).
//==============================================================================
void lblock_handle_bv_part(long part_beg, long part_end, long range_beg,
    long *range_gap, const gap_array_2n *block_gap, const bitvector *bv,
    const ranksel_support *bv_ranksel, long &res_sum, long &res_rank) {
  size_t excess_ptr = std::lower_bound(block_gap->m_excess.begin(),
      block_gap->m_excess.end(), part_beg) - block_gap->m_excess.begin();

  // Initialize j.
  long j = part_beg;

  // Compute gap[j].
  long gap_j = block_gap->m_count[j];
  while (excess_ptr < block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == j) {
    ++excess_ptr;
    gap_j += (1L << 16);
  }

  // Initialize sum.
  long sum = gap_j + 1;

  while (j != part_end - 1 && bv->get(j) == 1) {
    // Update j.
    ++j;

    // Compute gap[j].
    gap_j = block_gap->m_count[j];
    while (excess_ptr < block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == j) {
      ++excess_ptr;
      gap_j += (1L << 16);
    }

    // Update sum.
    sum += gap_j + 1;
  }
  if (bv->get(j) == 0) --sum;

  // Store gap[part_beg] + .. + gap[j] and bv.rank0(part_beg) (== bv.rank0(j)).
  res_sum = sum;
  res_rank = bv_ranksel->rank0(part_beg);

  if (j == part_end - 1)
    return;

  sum = 0L;
  long range_gap_ptr = res_rank + 1;
  while (j != part_end - 1) {
    // Update j.
    ++j;

    // Compute gap[j].
    gap_j = block_gap->m_count[j];
    while (excess_ptr < block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == j) {
      ++excess_ptr;
      gap_j += (1L << 16);
    }

    // Update sum.
    sum += gap_j + 1;

    // Update range_gap.
    if (bv->get(j) == 0) {
      range_gap[range_gap_ptr - range_beg] = sum - 1;
      ++range_gap_ptr;
      sum = 0L;
    }
  }

  if (bv->get(j) == 1)
    range_gap[range_gap_ptr - range_beg] = sum;
}


void lblock_async_write_code(unsigned char* &slab, long &length, std::mutex &mtx,
    std::condition_variable &cv, bool &avail, bool &finished, std::string filename) {
  while (true) {
    // Wait until the passive buffer is available.
    std::unique_lock<std::mutex> lk(mtx);
    while (!avail && !finished)
      cv.wait(lk);

    if (!avail && finished) {
      // We're done, terminate the thread.
      lk.unlock();
      return;
    }
    lk.unlock();

    // Safely write the data to disk.
    utils::add_objects_to_file(slab, length, filename);

    // Let the caller know that the I/O thread finished writing.
    lk.lock();
    avail = false;
    lk.unlock();
    cv.notify_one();
  }
}


//==============================================================================
// Given the gap array of the block (representation using 2 bytes per elements)
// and the gap array of the left half-block wrt right half-block (bitvector
// representation), compute the gap array (wrt tail) of the left half-block
// and write to a given file using v-byte encoding.
//
// The whole computation is performed under given ram budget. It is fully
// parallelized and uses asynchronous I/O as much as possible.
//==============================================================================
void compute_left_gap(long left_block_size, long right_block_size,
    const gap_array_2n *block_gap, bitvector *bv, std::string out_filename,
    long max_threads, long ram_budget) {
  long block_size = left_block_size + right_block_size;
  long left_gap_size = left_block_size + 1;

  // NOTE: we require that bv has room for one extra bit at the end
  //       which we use as a sentinel. The actual value of that bit
  //       prior to calling this function does not matter.
  bv->reset(block_size);
  long bv_size = block_size + 1;

  fprintf(stderr, "  Compute gap array for left half-block: ");
  long compute_gap_start = utils::wclock();

  //----------------------------------------------------------------------------
  // STEP 1: Preprocess left_block_gap_bv for rank and select queries,
  //         i.e., compute sparse_gap.
  //----------------------------------------------------------------------------
  ranksel_support *bv_ranksel = new ranksel_support(bv, bv_size, max_threads);


  //----------------------------------------------------------------------------
  // STEP 2: compute the values of the right gap array, one range at a time.
  //----------------------------------------------------------------------------
  long max_range_size = std::max(1L, ram_budget / (3L * (long)sizeof(long)));
  long n_ranges = (left_gap_size + max_range_size - 1) / max_range_size;

  // To ensure that asynchronous I/O is really taking
  // place, we try to make 8 parts.
  if (n_ranges < 8L) {
    max_range_size = (left_gap_size + 7L) / 8L;
    n_ranges = (left_gap_size + max_range_size - 1) / max_range_size;
  }

  long *range_gap = (long *)malloc(max_range_size * sizeof(long));
  unsigned char *active_vbyte_slab = (unsigned char *)malloc(max_range_size * sizeof(long));
  unsigned char *passive_vbyte_slab = (unsigned char *)malloc(max_range_size * sizeof(long));
  long active_vbyte_slab_length;
  long passive_vbyte_slab_length;

  // Used for communication with thread doing asynchronous writes.  
  std::mutex mtx;
  std::condition_variable cv;
  bool avail = false;
  bool finished = false;
  
  // Start the thread doing asynchronous writes.
  std::thread *async_writer = new std::thread(lblock_async_write_code,
      std::ref(passive_vbyte_slab), std::ref(passive_vbyte_slab_length),
      std::ref(mtx), std::ref(cv), std::ref(avail), std::ref(finished),
      out_filename);

  for (long range_id = 0L; range_id < n_ranges; ++range_id) {
    // Compute the range [range_beg..range_end) of values in the left gap
    // array (which is indexed [0..left_gap_size)).
    long range_beg = range_id * max_range_size;
    long range_end = std::min(range_beg + max_range_size, left_gap_size);
    long range_size = range_end - range_beg;

    // 2.a
    //
    // Find the section in the bitvector that contains
    // the bits necessary to compute the answer.
    long bv_section_beg = 0L;
    long bv_section_end = 0L;
    if (range_beg > 0)
      bv_section_beg = bv_ranksel->select0(range_beg - 1) + 1;
    bv_section_end = bv_ranksel->select0(range_end - 1) + 1;
    long bv_section_size = bv_section_end - bv_section_beg;

    // Split the current bitvector section into
    // equal parts. Each thread handles one part.
    long max_part_size = (bv_section_size + max_threads - 1) / max_threads;
    long n_parts = (bv_section_size + max_part_size - 1) / max_part_size;

    parallel_utils::parallel_fill<long>(range_gap, range_size, 0L, max_threads);

    // Allocate arrays used to store the answers for part boundaries.
    long *res_sum = new long[n_parts];
    long *res_rank = new long[n_parts];

    std::thread **threads = new std::thread*[n_parts];
    for (long t = 0; t < n_parts; ++t) {
      long part_beg = bv_section_beg + t * max_part_size;
      long part_end = std::min(part_beg + max_part_size, bv_section_end);

      threads[t] = new std::thread(lblock_handle_bv_part, part_beg, part_end, range_beg,
          range_gap, block_gap, bv, bv_ranksel, std::ref(res_sum[t]), std::ref(res_rank[t]));
    }

    for (long t = 0; t < n_parts; ++t) threads[t]->join();
    for (long t = 0; t < n_parts; ++t) delete threads[t];
    delete[] threads;

    // Update range_gap with values computed at part boundaries.
    for (long t = 0; t < n_parts; ++t)
      range_gap[res_rank[t] - range_beg] += res_sum[t];
    delete[] res_sum;
    delete[] res_rank;

    // 2.c
    //
    // Convert the range_gap to the slab of vbyte encoding.
    active_vbyte_slab_length = parallel_utils::convert_array_to_vbyte_slab(
        range_gap, range_size, active_vbyte_slab, max_threads);

    // 2.d
    //
    // Schedule asynchronous write of the slab.
    // First, wait for the I/O thread to finish writing.
    std::unique_lock<std::mutex> lk(mtx);
    while (avail == true)
      cv.wait(lk);

    // Set the new passive slab.
    std::swap(active_vbyte_slab, passive_vbyte_slab);
    passive_vbyte_slab_length = active_vbyte_slab_length;

    // Let the I/O thread know that the slab is waiting.
    avail = true;
    lk.unlock();
    cv.notify_one();
  }

  // Let the I/O thread know that we're done.
  std::unique_lock<std::mutex> lk(mtx);
  finished = true;
  lk.unlock();
  cv.notify_one();
  
  // Wait for the thread to finish.
  async_writer->join();

  // Clean up.
  delete async_writer;
  delete bv_ranksel;
  free(range_gap);
  free(active_vbyte_slab);
  free(passive_vbyte_slab);

  long double compute_gap_time = utils::wclock() - compute_gap_start;
  long double compute_gap_speed = (block_size / (1024.L * 1024)) / compute_gap_time;
  fprintf(stderr, "%.2Lfs (%.2LfMiB/s)\n", compute_gap_time, compute_gap_speed);
}

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_COMPUTE_LEFT_GAP_H_INCLUDED
