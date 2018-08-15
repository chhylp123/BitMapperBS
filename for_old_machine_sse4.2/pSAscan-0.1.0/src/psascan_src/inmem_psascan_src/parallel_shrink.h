/**
 * @file    src/psascan_src/inmem_psascan_src/parallel_shrink.h
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_SHRINK_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_SHRINK_H_INCLUDED

#include <algorithm>
#include <thread>


namespace psascan_private {
namespace inmem_psascan_private {

template<typename T, typename S>
void parallel_shrink_aux(T *src, S *dest, long length) {
  for (long i = 0; i < length; ++i)
    dest[i] = (S)src[i];
}


// Requires sizeof(T) > sizeof(S).
template<typename T, typename S>
S *parallel_shrink(T *tab, long length, long max_threads) {
  S *result = (S *)tab;

  long diff = (long)sizeof(T) - (long)sizeof(S);
  if (!diff) {
    fprintf(stderr, "Error: shrinking requires sizeof(T) > sizeof(S)\n");
    std::exit(EXIT_FAILURE);
  }

  // long threshold = (sizeof(T) + diff - 1) / diff;
  if (length < (1L << 20)/*threshold*/) {
    // Move the elelements sequentially.
    for (long i = 0; i < length; ++i)
      result[i] = (S)tab[i];

    return result;
  }

  // Compute the index of the smallest element (of type T)
  // that lies past the end of the last element of tab
  // (after converting all elemeents to type S).
  long bytes_after_shrinking = length * sizeof(S);
  long split = (bytes_after_shrinking + sizeof(T) - 1) / sizeof(T);

  // Recursively shrink the part up to (but excluding) split.
  parallel_shrink<T, S>(tab, split, max_threads);

  // Move the elements in the range [split, length) in parallel.
  // This is safe (no element overwriting) because of how we
  // computed the split.
  long elems = length - split;
  long max_block_size = (elems + max_threads - 1) / max_threads;
  long n_blocks = (elems + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = split + i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(parallel_shrink_aux<T, S>,
        tab + block_beg, result + block_beg, block_size);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;

  return result;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_SHRINK_H_INCLUDED
