/**
 * @file    src/psascan_src/ranksel_support.h
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

#ifndef __PSASCAN_SRC_RANKSEL_SUPPORT_H_INCLUDED
#define __PSASCAN_SRC_RANKSEL_SUPPORT_H_INCLUDED

#include <thread>
#include <algorithm>

#include "bitvector.h"


namespace psascan_private {

struct ranksel_support {
  //============================================================================
  // Compute sparse_rank[group_beg..group_end).
  //============================================================================
  static void process_group_of_chunks(long group_beg, long group_end,
      long chunk_size, long *sparse_rank, const bitvector *bv) {
    for (long chunk_id = group_beg; chunk_id < group_end; ++chunk_id) {
      long chunk_beg = chunk_id * chunk_size;
      long chunk_end = chunk_beg + chunk_size;

      sparse_rank[chunk_id] = bv->range_sum(chunk_beg, chunk_end);
    }
  }


  //============================================================================
  // Constructor.
  //============================================================================
  ranksel_support(const bitvector *bv, long length, long max_threads) {
    m_bv = bv;
    m_length = length;
    
    // 1
    //
    // Compute chunk size and allocate m_sparse_rank.
    m_chunk_size = std::min((1L << 20), (m_length + max_threads - 1) / max_threads);
    n_chunks = m_length / m_chunk_size;  // we exclude the last partial chunk
    m_sparse_rank = (long *)malloc((n_chunks + 1) * sizeof(long));

    // 2
    //
    // Compute the sum of 1-bits inside each chunk and write to m_sparse_rank.
    // Since there can be more chunks than threads, we split chunks
    // into groups and let each thread handle the group of chunks.
    long chunk_max_group_size = (n_chunks + max_threads - 1) / max_threads;
    long n_chunk_groups = (n_chunks + chunk_max_group_size - 1) / chunk_max_group_size;

    std::thread **threads = new std::thread*[n_chunk_groups];
    for (long t = 0; t < n_chunk_groups; ++t) {
      long chunk_group_beg = t * chunk_max_group_size;
      long chunk_group_end = std::min(chunk_group_beg + chunk_max_group_size, n_chunks);
      threads[t] = new std::thread(process_group_of_chunks, chunk_group_beg,
          chunk_group_end, m_chunk_size, m_sparse_rank, m_bv);
    }

    for (long t = 0; t < n_chunk_groups; ++t) threads[t]->join();
    for (long t = 0; t < n_chunk_groups; ++t) delete threads[t];
    delete[] threads;
    
    // 3
    //
    // Compute partial (exclusive) sum on m_sparse_rank.
    long ones = 0L;
    for (long i = 0; i < n_chunks; ++i) {
      long temp = m_sparse_rank[i];
      m_sparse_rank[i] = ones;
      ones += temp;
    }
    m_sparse_rank[n_chunks] = ones;
  }


  //============================================================================
  // Find the largest position j such that the number of 0s in bv[0..j) is <= i.
  // In other words, find the position of i-th 0-bit in bv (i = 0, 1, ..).
  // 0 <= i < number of 0-bits in bv.
  //============================================================================
  inline long select0(long i) const {
    // Fast-forward through chunks preceding the chunk with the answer.
    long j = 0L;
    while (j < n_chunks && ((j + 1) * m_chunk_size) - m_sparse_rank[j + 1] <= i)
      ++j;

    long zero_cnt_j = (j * m_chunk_size) - m_sparse_rank[j];
    j *= m_chunk_size;

    // Find the final position in a single chunk.
    while (zero_cnt_j + (1 - m_bv->get(j)) <= i)
      zero_cnt_j += (1 - m_bv->get(j++));

    return j;
  }


  //============================================================================
  // Find the largest position j such that the number of 1s in bv[0..j) is <= i.
  // In other words, find the position of i-th 1-bit in bv (i = 0, 1, ..).
  // 0 <= i < number of 1-bits in bv.
  //============================================================================
  inline long select1(long i) const {
    // Fast-forward through chunks preceding the chunk with the answer.
    long j = 0L;
    while (j < n_chunks && m_sparse_rank[j + 1] <= i)
      ++j;

    long rank_j = m_sparse_rank[j];
    j *= m_chunk_size;

    // Find the final position in a single chunk.
    while (rank_j + m_bv->get(j) <= i)
      rank_j += m_bv->get(j++);

    return j;
  }
  
  //============================================================================
  // Compute the number of 1-bits in bv[0..i) with the help of sparse_rank.
  // Note:
  // - i is an integer in the range from 0 to length of bv (inclusive),
  // - sparse_rank[k] = number of 1-bits in bv[0..k * chunk_size),
  //============================================================================
  inline long rank(long i) const {
    long j = i / m_chunk_size;
    long result = m_sparse_rank[j];
    j *= m_chunk_size;

    while (j < i)
      result += m_bv->get(j++);

    return result;
  }


  //============================================================================
  // Compute the number of 0-bits in bv[0..i).
  // 0 <= i <= m_length.
  //============================================================================
  inline long rank0(long i) const {
    return i - rank(i);
  }


  ~ranksel_support() {
    free(m_sparse_rank);
  }
  
  long m_length;      // length of bitvector
  long m_chunk_size;  // chunk size
  long n_chunks;      // number of chunks
  long *m_sparse_rank;

  const bitvector *m_bv;
};

}  // psascan_private

#endif  // __PSASCAN_SRC_RANKSEL_SUPPORT_H_INCLUDED
