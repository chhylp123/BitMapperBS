/**
 * @file    src/psascan_src/stream.h
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

#ifndef __PSASCAN_SRC_STREAM_H_INCLUDED
#define __PSASCAN_SRC_STREAM_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <mutex>
#include <algorithm>

#include "utils.h"
#include "rank.h"
#include "gap_buffer.h"
#include "update.h"
#include "stream_info.h"
#include "multifile.h"
#include "multifile_bit_stream_reader.h"
#include "async_multifile_bit_stream_reader.h"
#include "async_backward_skip_stream_reader.h"
#include "async_bit_stream_writer.h"


namespace psascan_private {

std::mutex stdout_mutex;

template<typename block_offset_type>
void parallel_stream(
    gap_buffer_poll<block_offset_type> *full_gap_buffers,
    gap_buffer_poll<block_offset_type> *empty_gap_buffers,
    long stream_block_beg,
    long stream_block_end,
    block_offset_type i,
    const long *count,
    block_offset_type whole_suffix_rank,
    const rank4n<> *rank,
    unsigned char last,
    std::string text_filename,
    long length,
    std::string &tail_gt_filename,
    stream_info *info,
    int thread_id,
    long gap_range_size,
    long gap_buf_size,
    const multifile *tail_gt_begin,
    long n_increasers) {

  static const int max_buckets = 4096;
  int *block_id_to_sblock_id = new int[max_buckets];

  long bucket_size = 1;
  long bucket_size_bits = 0;
  while ((gap_range_size + bucket_size - 1) / bucket_size > max_buckets)
    bucket_size <<= 1, ++bucket_size_bits;
  long n_buckets = (gap_range_size + bucket_size - 1) / bucket_size;
  int *block_count = new int[n_buckets];

  long max_buffer_elems = gap_buf_size / sizeof(block_offset_type);
  block_offset_type *temp = new block_offset_type[max_buffer_elems];
  int *oracle = new int[max_buffer_elems];

  static const long buffer_sample_size = 512;
  std::vector<block_offset_type> samples(buffer_sample_size);
  long *ptr = new long[n_increasers];
  block_offset_type *bucket_lbound = new block_offset_type[n_increasers + 1];

  typedef async_multifile_bit_stream_reader bit_stream_reader_type;
  typedef async_backward_skip_stream_reader<unsigned char> text_reader_type;
  typedef async_bit_stream_writer bit_stream_writer_type;

  text_reader_type *text_streamer = new text_reader_type(text_filename, length - stream_block_end, 4L << 20);
  bit_stream_writer_type *gt_out = new bit_stream_writer_type(tail_gt_filename, 1L << 20);
  bit_stream_reader_type gt_in(tail_gt_begin, length - stream_block_end, 1L << 20);

  long j = stream_block_end, dbg = 0L;
  while (j > stream_block_beg) {
    if (dbg > (1 << 26)) {
      info->m_mutex.lock();
      info->m_streamed[thread_id] = stream_block_end - j;
      info->m_update_count += 1;
      if (info->m_update_count == info->m_thread_count) {
        info->m_update_count = 0L;
        long double elapsed = utils::wclock() - info->m_timestamp;
        long total_streamed = 0L;

        for (long t = 0; t < info->m_thread_count; ++t)
          total_streamed += info->m_streamed[t];
        long double speed = (total_streamed / (1024.L * 1024)) / elapsed;

        stdout_mutex.lock();
        fprintf(stderr, "\r    Stream: %.2Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s",
            (total_streamed * 100.L) / info->m_tostream, elapsed, speed);
        stdout_mutex.unlock();
      }
      info->m_mutex.unlock();
      dbg = 0L;
    }

    // Get a gap buffer from the poll of empty buffers.
    std::unique_lock<std::mutex> lk(empty_gap_buffers->m_mutex);
    while (!empty_gap_buffers->available())
      empty_gap_buffers->m_cv.wait(lk);

    gap_buffer<block_offset_type> *b = empty_gap_buffers->get();
    lk.unlock();
    empty_gap_buffers->m_cv.notify_one(); // let others know they should re-check

    // Process buffer -- fill with gap values.
    long left = j - stream_block_beg;
    b->m_filled = std::min(left, b->m_size);
    dbg += b->m_filled;
    std::fill(block_count, block_count + n_buckets, 0);

    for (long t = 0L; t < b->m_filled; ++t, --j) {
      unsigned char c = text_streamer->read();

      gt_out->write(i > whole_suffix_rank);
      bool next_gt = (gt_in.read());

      int delta = (i > whole_suffix_rank && c == 0);
      i = (block_offset_type)(count[c] + rank->rank((long)i, c) - delta);
      if (c == last && next_gt) ++i;
      temp[t] = i;
      block_count[i >> bucket_size_bits]++;
    }

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
      // This is a fallback mechanism in case the quick partition failed.
      // It is not suppose to happen to often.

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

    // Add the buffer to the poll of full buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(full_gap_buffers->m_mutex);
    full_gap_buffers->add(b);
    lk2.unlock();
    full_gap_buffers->m_cv.notify_one();
  }

  delete text_streamer;
  delete gt_out;
  
  // Report that another worker thread has finished.
  std::unique_lock<std::mutex> lk(full_gap_buffers->m_mutex);
  full_gap_buffers->increment_finished_workers();
  lk.unlock();

  // Notify waiting update threads in case no more buffers
  // are going to be produces by worker threads.
  full_gap_buffers->m_cv.notify_one();

  delete[] block_count;
  delete[] block_id_to_sblock_id;
  delete[] temp;
  delete[] oracle;
  delete[] ptr;
  delete[] bucket_lbound;
}

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_STREAM_H_INCLUDED
