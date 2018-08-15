/**
 * @file    src/psascan_src/psascan.h
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

#ifndef __PSASCAN_SRC_PSASCAN_H_INCLUDED
#define __PSASCAN_SRC_PSASCAN_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/resource.h>

#include "partial_sufsort.h"
#include "merge.h"
#include "utils.h"
#include "uint40.h"
#include "half_block_info.h"


namespace psascan_private {

void pSAscan(std::string input_filename, std::string output_filename,
    std::string gap_filename, long ram_use, long max_threads,
    bool verbose, long gap_buf_size = (1L << 21)) {
  long n_gap_buffers = 2 * max_threads;
  if (ram_use < 6L) {
    fprintf(stderr, "Error: not enough memory to run pSAscan.\n");
    std::exit(EXIT_FAILURE);
  }
  
  // Turn paths absolute.
  input_filename = utils::absolute_path(input_filename);
  output_filename = utils::absolute_path(output_filename);
  gap_filename = utils::absolute_path(gap_filename);
  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input filename = %s\n", input_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Gap filename = %s\n", gap_filename.c_str());
  fprintf(stderr, "Input length = %ld (%.1LfMiB)\n", length, 1.L * length / (1L << 20));
  fprintf(stderr, "\n");

  long ram_for_threads = n_gap_buffers * gap_buf_size;  // for buffers
  if (ram_use / 5.2L < (long double)(1L << 31))  // for oracle
    ram_for_threads += max_threads * gap_buf_size;
  else ram_for_threads += ((4.L / 5) * max_threads) * gap_buf_size;
  ram_for_threads += max_threads * gap_buf_size;  // for temp
  ram_for_threads += max_threads * (6L << 20);  // for reader/writer buffers

  long ram_use_excluding_threads = ram_use - ram_for_threads;
  if (ram_use_excluding_threads < 6L) {
    long required_MiB = (ram_for_threads + (1L << 20) - 1) / (1L << 20);
    fprintf(stderr, "Error: not enough memory to start threads. You need "
        "at least %ldMiB\n", required_MiB + 1);
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "RAM budget = %ld (%.1LfMiB)\n", ram_use, 1.L * ram_use / (1L << 20));
  fprintf(stderr, "RAM budget (excluding threads) = %ld (%.1LfMiB)\n",
      ram_use_excluding_threads, 1.L * ram_use_excluding_threads / (1L << 20));
  long max_block_size = std::max(2L, (long)(ram_use_excluding_threads / 5.2L));

  fprintf(stderr, "Max block size = %ld (%.1LfMiB)\n\n", max_block_size, 1.L * max_block_size / (1L << 20));
  fprintf(stderr, "Parallel settings:\n");
  fprintf(stderr, "  #streaming threads = %ld\n", max_threads);
  fprintf(stderr, "  #gap buffers = %ld\n", n_gap_buffers);
  fprintf(stderr, "  gap buffer size = %ld\n\n", gap_buf_size);

  // Check if the maximum number of open files
  // is large enough for the merging to work.
  long n_half_blocks_estimated = 2L * (length / max_block_size + 1);
  long merge_max_open_files_estimated = 2L * n_half_blocks_estimated;
  long stream_max_open_files_estimated = 3L * max_threads + 1;
  long max_open_files_estimated = std::max(merge_max_open_files_estimated, stream_max_open_files_estimated);
  rlimit rlimit_res;
  if (!getrlimit(RLIMIT_NOFILE, &rlimit_res) &&
      (long)rlimit_res.rlim_cur < max_open_files_estimated) {
    fprintf(stderr,
"\nError: the limit on the maximum number of open files is too small\n"
"(current limit = %ld, required limit = %ld). See the README for\n"
"more information.\n",
        (long)rlimit_res.rlim_cur, max_open_files_estimated);
    std::exit(EXIT_FAILURE);
  }

  long double start = utils::wclock();
  if (max_block_size < (1L << 31)) {
    std::vector<half_block_info<int> > hblock_info = partial_sufsort<int>(input_filename,
        output_filename, gap_filename, length, max_block_size, ram_use, max_threads, gap_buf_size, verbose);
    merge<int>(output_filename, ram_use, hblock_info);
  } else {
    std::vector<half_block_info<uint40> > hblock_info = partial_sufsort<uint40>(input_filename,
        output_filename, gap_filename, length, max_block_size, ram_use, max_threads, gap_buf_size, verbose);
    merge<uint40>(output_filename, ram_use, hblock_info);
  }
  long double total_time = utils::wclock() - start;

  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  elapsed time: %.2Lfs (%.4Lfs/MiB)\n", total_time, total_time / ((1.L * length) / (1L << 20)));
  fprintf(stderr, "  speed: %.2LfMiB/s\n", ((1.L * length) / (1L << 20)) / total_time);
}

}  // namespace psascan_private


// The main function.
void pSAscan(std::string input_filename, std::string output_filename,
    std::string gap_filename, long ram_use, long max_threads, bool verbose) {
  psascan_private::pSAscan(input_filename, output_filename,
      gap_filename, ram_use, max_threads, verbose);
}

#endif  // __PSASCAN_SRC_PSASCAN_H_INCLUDED
