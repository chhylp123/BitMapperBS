/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_stream.h
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_STREAM_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_STREAM_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <queue>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "../bitvector.h"
#include "../gap_buffer.h"
#include "../utils.h"
#include "rank.h"
#include "inmem_update.h"


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// The main streaming function.
//
// Note:
//    * it reads and writes bits in range
//      [stream_block_beg..stream_block_end) from gt bitvector right to left.
//==============================================================================
template<typename rank_type, typename block_offset_type>
void inmem_parallel_stream(
    const unsigned char *text,
    long text_length,
    long stream_block_beg,
    long stream_block_end,
    unsigned char last,
    const long *count,
    gap_buffer_poll<block_offset_type> *full_gap_buffers,
    gap_buffer_poll<block_offset_type> *empty_gap_buffers,
    block_offset_type i,
    block_offset_type i0,
    const rank_type *rank,
    long gap_range_size,
    long n_increasers,
    bitvector *gt,
    block_offset_type *temp,
    int *oracle,
    bool need_gt) {

  //----------------------------------------------------------------------------
  // STEP 1: initialize structures necessary to do the buffer partitions.
  //----------------------------------------------------------------------------
  static const int max_buckets = 4096;
  int *block_id_to_sblock_id = new int[max_buckets];

  long bucket_size = 1;
  long bucket_size_bits = 0;
  while ((gap_range_size + bucket_size - 1) / bucket_size > max_buckets)
    bucket_size <<= 1, ++bucket_size_bits;
  long n_buckets = (gap_range_size + bucket_size - 1) / bucket_size;
  int *block_count = new int[n_buckets];

  static const long buffer_sample_size = 512;
  std::vector<block_offset_type> samples(buffer_sample_size);
  long *ptr = new long[n_increasers];
  block_offset_type *bucket_lbound = new block_offset_type[n_increasers + 1];

  //----------------------------------------------------------------------------
  // STEP 2: perform the actual streaming.
  //----------------------------------------------------------------------------
  long j = stream_block_end;
  bool gt_bit = gt->get(text_length - j);
  while (j > stream_block_beg) {
    // 2.a
    //
    // Get a buffer from the poll of empty buffers.
    std::unique_lock<std::mutex> lk(empty_gap_buffers->m_mutex);
    while (!empty_gap_buffers->available()) empty_gap_buffers->m_cv.wait(lk);
    gap_buffer<block_offset_type> *b = empty_gap_buffers->get();
    lk.unlock();
    empty_gap_buffers->m_cv.notify_one();

    // 2.b
    //
    // Process buffer, i.e., fill with gap values.
    long left = j - stream_block_beg;
    b->m_filled = std::min(left, b->m_size);
    std::fill(block_count, block_count + n_buckets, 0);

    if (need_gt) {
      for (long t = 0; t < b->m_filled; ++t) {
        bool new_gt_bit = (i > i0);
        if (new_gt_bit) gt->set(text_length - j);
        else gt->reset(text_length - j);

        unsigned char c = text[j - 1];

        // Compute new i.
        int delta = (new_gt_bit && c == 0);
        i = (block_offset_type)(count[c] + rank->rank(i, c) - delta);
        if (c == last && gt_bit) ++i;

        temp[t] = i;
        block_count[i >> bucket_size_bits]++;

        --j;
        gt_bit = gt->get(text_length - j);
      }
    } else {
      for (long t = 0; t < b->m_filled; ++t) {
        bool new_gt_bit = (i > i0);

        unsigned char c = text[j - 1];

        // Compute new i.
        int delta = (new_gt_bit && c == 0);
        i = (block_offset_type)(count[c] + rank->rank(i, c) - delta);
        if (c == last && gt_bit) ++i;

        temp[t] = i;
        block_count[i >> bucket_size_bits]++;

        --j;
        gt_bit = gt->get(text_length - j);
      }

    }

    // 2.c
    //
    // Partition the buffer into equal n_increasers parts.

    // Compute super-buckets.
    long ideal_sblock_size = (b->m_filled + n_increasers - 1) / n_increasers;
    long max_sbucket_size = 0;
    long bucket_id_beg = 0;
    for (long t = 0; t < n_increasers; ++t) {
      long bucket_id_end = bucket_id_beg, size = 0L;
      while (bucket_id_end < n_buckets && size < ideal_sblock_size)
        size += block_count[bucket_id_end++];
      b->sblock_size[t] = size;
      max_sbucket_size = std::min(max_sbucket_size, size);
      for (long id = bucket_id_beg; id < bucket_id_end; ++id)
        block_id_to_sblock_id[id] = t;
      bucket_id_beg = bucket_id_end;
    }

    if (max_sbucket_size < 4L * ideal_sblock_size) {
      // The quick partition was good enough.
      for (long t = 0, curbeg = 0; t < n_increasers; curbeg += b->sblock_size[t++])
        b->sblock_beg[t] = ptr[t] = curbeg;

      // Permute the elements of the buffer.
      for (long t = 0; t < b->m_filled; ++t) {
        long id = (temp[t] >> bucket_size_bits);
        long sblock_id = block_id_to_sblock_id[id];
        oracle[t] = ptr[sblock_id]++;
      }

      for (long t = 0; t < b->m_filled; ++t) {
        long addr = oracle[t];
        b->m_content[addr] = temp[t];
      }
    } else {
      // Repeat the partition into sbuckets, this time using random sample.
      // This is a fallback mechanism in case the quick partition failed,
      // and is expected to happen very rarely.
      
      // Compute random sample of elements in the buffer.
      for (long t = 0; t < buffer_sample_size; ++t)
        samples[t] = temp[utils::random_long(0L, b->m_filled - 1)];
      std::sort(samples.begin(), samples.end());
      samples.erase(std::unique(samples.begin(), samples.end()), samples.end());

      // Compute bucket boundaries (lower bound is enough).
      std::fill(bucket_lbound, bucket_lbound + n_increasers + 1, gap_range_size);

      long step = (samples.size() + n_increasers - 1) / n_increasers;
      for (size_t t = 1, p = step; p < samples.size(); ++t, p += step)
        bucket_lbound[t] = (samples[p - 1] + samples[p] + 1) / 2;
      bucket_lbound[0] = 0;

      // Compute bucket sizes and sblock id into oracle array.
      std::fill(b->sblock_size, b->sblock_size + n_increasers, 0L);
      for (long t = 0; t < b->m_filled; ++t) {
        block_offset_type x = temp[t];
        int id = n_increasers;
        while (bucket_lbound[id] > x) --id;
        oracle[t] = id;
        b->sblock_size[id]++;
      }

      // Permute elements into their own buckets using oracle.
      for (long t = 0, curbeg = 0; t < n_increasers; curbeg += b->sblock_size[t++])
        b->sblock_beg[t] = ptr[t] = curbeg;

      for (long t = 0; t < b->m_filled; ++t) {
        long sblock_id = oracle[t];
        oracle[t] = ptr[sblock_id]++;
      }

      for (long t = 0; t < b->m_filled; ++t) {
        long addr = oracle[t];
        b->m_content[addr] = temp[t];
      }
    }

    // 2.d
    //
    // Add the buffer to the poll of full buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(full_gap_buffers->m_mutex);
    full_gap_buffers->add(b);
    lk2.unlock();
    full_gap_buffers->m_cv.notify_one();
  }

  //---------------------------------------------------------------------------
  // STEP 3: Clean up.
  //---------------------------------------------------------------------------

  // Report that another thread has finished.
  std::unique_lock<std::mutex> lk(full_gap_buffers->m_mutex);
  full_gap_buffers->increment_finished_workers();
  lk.unlock();

  // Notify waiting update threads in case no more buffers
  // are going to be produced by streaming threads.
  full_gap_buffers->m_cv.notify_one();

  delete[] block_count;
  delete[] block_id_to_sblock_id;
  delete[] ptr;
  delete[] bucket_lbound;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_STREAM_H_INCLUDED
