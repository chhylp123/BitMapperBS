/**
 * @file    src/psascan_src/compute_gap.h
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

#ifndef __PSASCAN_SRC_COMPUTE_GAP_H_INCLUDED
#define __PSASCAN_SRC_COMPUTE_GAP_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <algorithm>
#include <vector>

#include "utils.h"
#include "rank.h"
#include "gap_array.h"
#include "gap_buffer.h"
#include "stream.h"
#include "update.h"
#include "stream_info.h"
#include "multifile.h"


namespace psascan_private {

//==============================================================================
// Compute the gap for an arbitrary range of suffixes of tail. This version is
// more general, and can be used also when processing half-blocks.
//==============================================================================
template<typename block_offset_type>
void compute_gap(const rank4n<> *rank, buffered_gap_array *gap,
    long tail_begin, long tail_end, long text_length, long max_threads,
    long block_isa0, long gap_buf_size, unsigned char block_last_symbol,
    std::vector<long> initial_ranks, std::string text_filename, std::string output_filename,
    const multifile *tail_gt_begin_rev, multifile *newtail_gt_begin_rev) {
  long tail_length = tail_end - tail_begin;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  fprintf(stderr, "    Stream:");
  long double stream_start = utils::wclock();

  // 1
  //
  // Get symbol counts of a block and turn into exclusive partial sum.
  long *count = new long[256];
  std::copy(rank->m_count, rank->m_count + 256, count);
  ++count[block_last_symbol];
  --count[0];
  for (long j = 0, s = 0, t; j < 256; ++j) {
    t = count[j];
    count[j] = s;
    s += t;
  }

  // 2
  //
  // Allocate gap buffers.
  long n_gap_buffers = 2 * max_threads;
  gap_buffer<block_offset_type> **gap_buffers = new gap_buffer<block_offset_type>*[n_gap_buffers];
  for (long i = 0L; i < n_gap_buffers; ++i)
    gap_buffers[i] = new gap_buffer<block_offset_type>(gap_buf_size, max_threads);

  // 3
  //
  // Create poll of empty and full buffers.
  gap_buffer_poll<block_offset_type> *empty_gap_buffers = new gap_buffer_poll<block_offset_type>();
  gap_buffer_poll<block_offset_type> *full_gap_buffers = new gap_buffer_poll<block_offset_type>(n_threads);

  // 4
  //
  // Add all buffers to the poll of empty buffers.
  for (long i = 0L; i < n_gap_buffers; ++i)
    empty_gap_buffers->add(gap_buffers[i]);

  // 5
  //
  // Start threads doing the backward search.
  stream_info info(n_threads, tail_length);
  std::thread **streamers = new std::thread*[n_threads];
  std::vector<std::string> gt_filenames(n_threads);

  for (long t = 0L; t < n_threads; ++t) {
    long stream_block_beg = tail_begin + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, tail_end);

    gt_filenames[t] = output_filename + ".gt_tail." + utils::random_string_hash();
    newtail_gt_begin_rev->add_file(text_length - stream_block_end, text_length - stream_block_beg, gt_filenames[t]);

    streamers[t] = new std::thread(parallel_stream<block_offset_type>, full_gap_buffers, empty_gap_buffers, stream_block_beg,
        stream_block_end, initial_ranks[t], count, block_isa0, rank, block_last_symbol, text_filename, text_length,
        std::ref(gt_filenames[t]), &info, t, gap->m_length, gap_buf_size, tail_gt_begin_rev, max_threads);
  }

  // 6
  //
  // Start threads doing the gap array updates.
  std::thread *updater = new std::thread(gap_updater<block_offset_type>,
        full_gap_buffers, empty_gap_buffers, gap, max_threads);

  // 7
  //
  // Wait for all threads to finish.
  for (long i = 0L; i < n_threads; ++i) streamers[i]->join();
  updater->join();

  // 8
  //
  // Clean up.
  for (long i = 0L; i < n_threads; ++i) delete streamers[i];
  for (long i = 0L; i < n_gap_buffers; ++i) delete gap_buffers[i];
  delete updater;
  delete[] streamers;
  delete[] gap_buffers;
  delete empty_gap_buffers;
  delete full_gap_buffers;
  delete[] count;

  // 9
  //
  // Print summary and exit.
  long double stream_time = utils::wclock() - stream_start;
  long double speed = (tail_length / (1024.L * 1024)) / stream_time;
  fprintf(stderr,"\r    Stream: 100.0%%. Time: %.2Lfs. Speed: %.2LfMiB/s\n",
      stream_time, speed);
}

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_COMPUTE_GAP_H_INCLUDED
