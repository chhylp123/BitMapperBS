/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_gap_array.h
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_SASCAN_INMEM_GAP_ARRAY_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_SASCAN_INMEM_GAP_ARRAY_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <mutex>
#include <stack>
#include <thread>


namespace psascan_private {
namespace inmem_psascan_private {

struct inmem_gap_array {
  unsigned char *m_count;
  long m_length;

  std::vector<long> m_excess;
  std::mutex m_excess_mutex;

  inmem_gap_array(long length)
    : m_length(length) {
    m_count = (unsigned char *)calloc(m_length, sizeof(unsigned char));
  }

  ~inmem_gap_array() {
    free(m_count);
  }
  
  //==============================================================================
  // Find and smallest j such that j + gap[0] + .. + gap[j] >= a. Store
  // the value of j into b and gap[0] + .. + gap[j] into c. To speed up the
  // algorithm, we have array gapsum defined as
  //
  //    gapsum[i] = gap[0] + .. + gap[i * block_size - 1].
  //
  //==============================================================================
  static void answer_single_gap_query(const inmem_gap_array *gap, long block_size,
      const long *gapsum, long a, long &b, long &c) {
    long n_blocks = (gap->m_length + block_size - 1) / block_size;

    // Find the block containing the correct index. To do that  find the largest
    // j such that gapsum[j] + block_size * j - 1 < a and start searching from
    // j * block_size.
    long j = 0;
    while (j + 1 < n_blocks && gapsum[j + 1] + block_size * (j + 1) - 1 < a) ++j;
    // Invariant: the j we are searching for is > j * block_size - 1.

    long sum = gapsum[j];
    j = block_size * j;
    size_t excess_ptr = std::lower_bound(gap->m_excess.begin(),
        gap->m_excess.end(), j) - gap->m_excess.begin();
    while (true) {
      // Invariant: sum = gap[0] + .. + gap[j - 1].
      // Compute gap[j] using small gap array representation.
      long gap_j = gap->m_count[j];
      while (excess_ptr < gap->m_excess.size() && gap->m_excess[excess_ptr] == j) {
        gap_j += 256L;
        ++excess_ptr;
      }

      if (j + sum + gap_j >= a) { b = j; c = sum + gap_j; return; }
      else { sum += gap_j; ++j; }
    }
  }

  //==============================================================================
  // Compute gap[0] + gap[1] + .. + gap[j - 1] with the help of gapsum array.
  //==============================================================================
  static long compute_sum3(const inmem_gap_array *gap, long j,
      long max_block_size, long *gapsum) {
    long block_id = j / max_block_size;
    long result = gapsum[block_id];

    long scan_beg = block_id * max_block_size;
    long scan_end = j;
    long occ = std::upper_bound(gap->m_excess.begin(), gap->m_excess.end(), scan_end - 1)
             - std::lower_bound(gap->m_excess.begin(), gap->m_excess.end(), scan_beg);
    result += 256L * std::max(0L, occ);
    for (long i = block_id * max_block_size; i < j; ++i)
      result += gap->m_count[i];

    return result;
  }

  //==============================================================================
  // Compute sum of gap values for blocks in range [range_beg..range_end).
  // The sum for each block is stored in gapsum array.
  //==============================================================================
  static void compute_sum2(const inmem_gap_array *gap, long range_beg,
      long range_end, long max_block_size, long *gapsum) {
    for (long block_id = range_beg; block_id < range_end; ++block_id) {
      long block_beg = block_id * max_block_size;
      long block_end = std::min(block_beg + max_block_size, gap->m_length);

      // Process block.
      long occ = std::upper_bound(gap->m_excess.begin(), gap->m_excess.end(), block_end - 1)
               - std::lower_bound(gap->m_excess.begin(), gap->m_excess.end(), block_beg);
      long block_gap_sum = 256L * std::max(0L, occ);
      for (long j = block_beg; j < block_end; ++j)
        block_gap_sum += gap->m_count[j];

      gapsum[block_id] = block_gap_sum;
    }
  }

  //==============================================================================
  // Parallel computaton of answers to n_queries queries of the form:
  // What is the smallest j such that j + gap[0] + .. + gap[j] >= a[i]"
  //   - the answer to i-th query is stored in b[i]
  //   - in addition we also return gap[0] + .. + gap[j] in c[i]
  //
  // To do that we first split the gap array into blocks of size of about
  // length / max_threads and (in parallel) compute sums of gap values inside
  // these blocks. We the accumulate these sums into array of prefix sums.
  //
  // To answer each of the queries we start a separate thread. Each thread uses
  // the partial sums of gap array at block boundaries to find a good starting
  // point for search and then scans the gap array from there.
  //==============================================================================
  long answer_queries(long n_queries, const long *a, long *b, long *c,
      long max_threads, long i0) const {
    //----------------------------------------------------------------------------
    // STEP 1: split gap array into at most max_threads blocks
    // and in parallel compute sum of values inside each block.
    //----------------------------------------------------------------------------
    long max_block_size = std::min(4L << 20, (m_length + max_threads - 1) / max_threads);
    long n_blocks = (m_length + max_block_size - 1) / max_block_size;
    long *gapsum = new long[n_blocks];
  
    // Each thread handles range of blocks.
    long range_size = (n_blocks + max_threads - 1) / max_threads;
    long n_ranges = (n_blocks + range_size - 1) / range_size;
    std::thread **threads = new std::thread*[max_threads];
    for (long range_id = 0; range_id < n_ranges; ++range_id) {
      long range_beg = range_id * range_size;
      long range_end = std::min(range_beg + range_size, n_blocks);

      threads[range_id] = new std::thread(compute_sum2, this,
          range_beg, range_end, max_block_size, gapsum);
    }
    for (long i = 0; i < n_ranges; ++i) threads[i]->join();
    for (long i = 0; i < n_ranges; ++i) delete threads[i];
    delete[] threads;

    //----------------------------------------------------------------------------
    // STEP 2: compute partial sum from block counts.
    //----------------------------------------------------------------------------
    for (long i = 0, s = 0, t; i < n_blocks; ++i)
      { t = gapsum[i]; gapsum[i] = s; s += t; }

    //----------------------------------------------------------------------------
    // STEP 3: Answer the queries in parallel.
    //----------------------------------------------------------------------------
    threads = new std::thread*[n_queries];
    for (long i = 0; i < n_queries; ++i)
      threads[i] = new std::thread(answer_single_gap_query, this,
        max_block_size, gapsum, a[i], std::ref(b[i]), std::ref(c[i]));
    for (long i = 0; i < n_queries; ++i) threads[i]->join();
    for (long i = 0; i < n_queries; ++i) delete threads[i];
    delete[] threads;

    long result = -1;
    if (i0 != -1) 
      result = compute_sum3(this, i0 + 1, max_block_size, gapsum);
  
    delete[] gapsum;
  
    return result;
  }
};

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_GAP_ARRAY_H_INCLUDED
