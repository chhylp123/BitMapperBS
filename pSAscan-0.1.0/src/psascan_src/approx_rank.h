/**
 * @file    src/psascan_src/approx_rank.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section DESCRIPTION
 *
 * The approximate rank data structure. Based on the 'sparse-LF'
 * data structure described in:
 *
 *   Dominik Kempa, Simon J. Puglisi:
 *   Lempel-Ziv Factorization: Simple, Fast, Practical.
 *   In Proc. ALENEX 2013, p. 103-112.
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

#ifndef __PSASCAN_SRC_APPROX_RANK_H_INCLUDED
#define __PSASCAN_SRC_APPROX_RANK_H_INCLUDED

#include <thread>
#include <algorithm>


namespace psascan_private {

template<long k_sampling_rate_log>
class approx_rank {
  private:
    long *m_list_size;
    long **m_list;

    static const long k_sampling_rate;
    static const long k_sampling_rate_mask;

  public:
    long *m_count;

  private:
    static void compute_symbol_count_aux(const unsigned char *text, long beg,
        long end, long *symbol_count) {
      for (long j = beg; j < end; ++j)
        ++symbol_count[text[j]];
    }

    static void compute_occ_list_aux(const unsigned char *text, long beg,
        long end, long *symbol_count, long **list) {
      // Compute where to start writing positions for each symbol.
      long *ptr = new long[256];
      for (long c = 0; c < 256; ++c)
        ptr[c] = (symbol_count[c] + k_sampling_rate - 1) / k_sampling_rate;

      // Add occurrences in the block to the lists.
      for (long j = beg; j < end; ++j) {
        unsigned char c = text[j];
        if (!((symbol_count[c]++) & k_sampling_rate_mask))
          list[c][ptr[c]++] = j;
      }

      // Clean up.
      delete[] ptr;
    }

  public:
    approx_rank(const unsigned char *text, long length, long max_threads) {
      // Compute symbol counts in each block.
      long max_block_size = (length + max_threads - 1) / max_threads;
      long n_threads = (length + max_block_size - 1) / max_block_size;
      long **symbol_count = new long*[n_threads];
      for (long j = 0; j < n_threads; ++j) {
        symbol_count[j] = new long[256];
        std::fill(symbol_count[j], symbol_count[j] + 256, 0L);
      }

      std::thread **threads = new std::thread*[n_threads];
      for (long t = 0; t < n_threads; ++t) {
        long block_beg = t * max_block_size;
        long block_end = std::min(block_beg + max_block_size, length);

        threads[t] = new std::thread(compute_symbol_count_aux,
            text, block_beg, block_end, symbol_count[t]);
      }

      for (long t = 0; t < n_threads; ++t) threads[t]->join();
      for (long t = 0; t < n_threads; ++t) delete threads[t];

      // Compute (exclusive) partial sums over symbol counts.
      m_count = new long[256];
      std::fill(m_count, m_count + 256, 0L);
      long *temp_count = new long[256];
      for (long i = 0; i < n_threads; ++i) {
        std::copy(symbol_count[i], symbol_count[i] + 256, temp_count);
        std::copy(m_count, m_count + 256, symbol_count[i]);
        for (long j = 0; j < 256; ++j)
          m_count[j] += temp_count[j];
      }
      delete[] temp_count;

      // Compute sizes and allocate occurrences lists.
      m_list_size = new long[256];
      m_list = new long*[256];
      for (long i = 0; i < 256; ++i) {
        m_list_size[i] = (m_count[i] + k_sampling_rate - 1) / k_sampling_rate;
        if (m_list_size[i]) m_list[i] = new long[m_list_size[i]];
        else m_list[i] = NULL;
      }

      for (long t = 0; t < n_threads; ++t) {
        long block_beg = t * max_block_size;
        long block_end = std::min(block_beg + max_block_size, length);

        threads[t] = new std::thread(compute_occ_list_aux, text,
            block_beg, block_end, symbol_count[t], m_list);
      }

      for (long t = 0; t < n_threads; ++t) threads[t]->join();
      for (long t = 0; t < n_threads; ++t) delete threads[t];
      delete[] threads;


      // Clean up.
      for (long j = 0; j < n_threads; ++j)
        delete[] symbol_count[j];
      delete[] symbol_count;
    }

    inline long rank(long i, unsigned char c) const {
      if (i <= 0 || (!m_list_size[c]) || m_list[c][0] >= i)
        return 0L;

      long left = 0, right = m_list_size[c];
      while (left + 1 != right) {
        // Invariant: the answer is in range [left..right).
        long mid = (left + right) / 2;
        if (m_list[c][mid] <= i) left = mid;
        else right = mid;
      }
      return (left << k_sampling_rate_log);
    }

    ~approx_rank() {
      delete[] m_count;
      delete[] m_list_size;
      for (long j = 0; j < 256; ++j) {
        if (m_list[j])
          delete[] m_list[j];
      }
      delete[] m_list;
    }
};

template<long k_sampling_rate_log>
const long approx_rank<k_sampling_rate_log>::k_sampling_rate = (1L << k_sampling_rate_log);

template<long k_sampling_rate_log>
const long approx_rank<k_sampling_rate_log>::k_sampling_rate_mask = (1L << k_sampling_rate_log) - 1;

}  // namespace psascan_private

#endif // __PSASCAN_SRC_APPROX_RANK_H_INCLUDED
