/**
 * @file    src/psascan_src/inmem_psascan_src/sparse_isa.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section DESCRIPTION
 *
 * Sparse ISA encoding based on the ISAs algorithm computing
 * Lempel-Ziv (LZ77) factorization described in
 *
 *   Dominik Kempa, Simon J. Puglisi:
 *   Lempel-Ziv factorization: Simple, fast, practical.
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_H_INCLUDED

#include <algorithm>
#include <thread>


namespace psascan_private {
namespace inmem_psascan_private {

template<typename pagearray_type, typename rank_type, unsigned isa_sampling_rate_log>
struct sparse_isa {
  static const unsigned isa_sampling_rate = (1U << isa_sampling_rate_log);
  static const unsigned isa_sampling_rate_mask = isa_sampling_rate - 1;
  static const long k_sigma = 256;

  static void compute_sparse_isa_aux(const pagearray_type &bwtsa, long block_beg,
      long block_end, long psa_size, long *sparse_isa, long &last) {
    for (long j = block_beg; j < block_end; ++j) {
      long sa_j = bwtsa[j].sa;
      if (!(sa_j & isa_sampling_rate_mask))
        sparse_isa[sa_j >> isa_sampling_rate_log] = j;
      if (sa_j == psa_size - 1) last = j;
    }
  }

  sparse_isa(const pagearray_type *bwtsa, const unsigned char *text,
      const rank_type *rank, long length, long i0, long max_threads) {
    m_bwtsa = bwtsa;
    m_length = length;
    m_rank = rank;
    m_text = text;
    m_i0 = i0;

    long elems = (m_length + isa_sampling_rate - 1) / isa_sampling_rate + 1;
    m_sparse_isa = (long *)malloc(elems * sizeof(long));

    long max_block_size = (m_length + max_threads - 1) / max_threads;
    long n_blocks = (m_length + max_block_size - 1) / max_block_size;

    std::thread **threads = new std::thread*[n_blocks];
    for (long t = 0; t < n_blocks; ++t) {
      long block_beg = t * max_block_size;
      long block_end = std::min(block_beg + max_block_size, m_length);

      threads[t] = new std::thread(compute_sparse_isa_aux, std::ref(*m_bwtsa),
          block_beg, block_end, m_length, m_sparse_isa, std::ref(m_last_isa));
    }

    for (long t = 0; t < n_blocks; ++t) threads[t]->join();
    for (long t = 0; t < n_blocks; ++t) delete threads[t];
    delete[] threads;

    m_count = (long *)malloc(k_sigma * sizeof(long));
    std::copy(rank->m_count, rank->m_count + k_sigma, m_count);
    ++m_count[text[length - 1]];
    --m_count[0];

    for (long i = 0, s = 0; i < k_sigma; ++i) {
      long t = m_count[i];
      m_count[i] = s;
      s += t;
    }
  }

  inline long query(long j) const {
    long isa_i;
    long i = ((j + isa_sampling_rate - 1) >> isa_sampling_rate_log);
    if ((i << isa_sampling_rate_log) < m_length) {
      isa_i = m_sparse_isa[i];
      i <<= isa_sampling_rate_log;
    } else {
      isa_i = m_last_isa;
      i = m_length - 1;
    }

    while (i != j) {
      // Compute ISA[i - 1] from ISA[i].
      // Invariant:
      //   isa_i = ISA[i]
      //   j <= i
      unsigned char c = m_text[i - 1];
      int delta = (isa_i > m_i0 && c == 0);

      isa_i = m_count[c] + m_rank->rank(isa_i, c) - delta;
      if (isa_i < 0 || ((long)((*m_bwtsa)[isa_i].sa)) != i - 1)
        ++isa_i;

      --i;
    }

    return isa_i;
  }

  ~sparse_isa() {
    free(m_sparse_isa);
    free(m_count);
  }


private:
  long m_length;
  long m_last_isa;
  long m_i0;

  long *m_count;
  long *m_sparse_isa;
  
  const unsigned char *m_text;
  const pagearray_type *m_bwtsa;
  const rank_type *m_rank;
};

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_H_INCLUDED
