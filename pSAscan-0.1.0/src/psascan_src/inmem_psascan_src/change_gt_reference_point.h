/**
 * @file    src/psascan_src/inmem_psascan_src/change_gt_reference_point.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section DESCRIPTION
 *
 * In-place computation of gt_begin bitvector from gt_end bitvector
 * (reversed). The procedure uses the string range matching algorithm
 * described in
 *
 *   Juha Karkkainen, Dominik Kempa, Simon J. Puglisi:
 *   String Range Matching.
 *   In Proc. CPM 2014, p. 232-241.
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_CHANGE_GT_REFERENCE_POINT_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_CHANGE_GT_REFERENCE_POINT_H_INCLUDED

#include <cstring>
#include <algorithm>
#include <thread>

#include "../bitvector.h"
#include "srank_aux.h"


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// Compute range [microblock_beg..microblock_end) of bits in the output
// bitvector gt_out.
//==============================================================================
void gt_end_to_gt_begin_aux(const unsigned char *text, long text_length,
    long block_beg, long block_end, bitvector *gt) {
  long block_size = block_end - block_beg;
  const unsigned char *pat = text + block_beg, *txt = pat;

  long i = 1, el = 0L, s = 0L, p = 0L;
  long i_max = i, el_max = 0L, s_max = 0L, p_max = 0L;

  long rev_end = text_length - block_beg;
  while (i < block_size) {
    // Compute lcp(text[left_block_beg..), text[left_block_beg+i..),
    // but compare not more than left_block_size symbols (we have gt
    // to resolve the long comparisons).
    while (block_beg + i + el < block_end && txt[i + el] == pat[el])
      update_ms(pat, ++el, s, p);

    if (((block_beg + i + el != block_end && txt[i + el] > pat[el]) ||
         (block_beg + i + el == block_end && !gt->get(rev_end - i))))
      gt->set(rev_end - i);
    else gt->reset(rev_end - i);

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      i_max = i;
    }

    if (el < 100) {
      ++i;
      el = 0;
    } else if (p > 0L && (p << 2) <= el && !memcmp(pat, pat + p, s)) {
      long maxk = std::min(block_size - i, p);
      for (long k = 1L; k < maxk; ++k) {
        if (gt->get(rev_end - (j + k))) gt->set(rev_end - (i + k));
        else gt->reset(rev_end - (i + k));
      }

      i += p;
      el -= p;
    } else {
      long h = (el >> 2) + 1L;
      long maxk = std::min(h, block_size - i);
      for (long k = 1L; k < maxk; ++k) {
        if (gt->get(rev_end - (j + k))) gt->set(rev_end - (i + k));
        else gt->reset(rev_end - (i + k));
      }

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }
}


//==============================================================================
// Change gt_end bitvector into gt_begin using string range matching.
//==============================================================================
void gt_end_to_gt_begin(const unsigned char *text, long text_length,
    bitvector *gt, long max_block_size) {
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: Compute the last bit in every block.
  //----------------------------------------------------------------------------
  for (long i = 0; i < n_blocks; ++i) {
    long block_end = text_length - (n_blocks - 1 - i) * max_block_size;
    long rev_beg = text_length - block_end;
    gt->flip(rev_beg);
  }

  //----------------------------------------------------------------------------
  // STEP 2: compute remaining bits in every block.
  //----------------------------------------------------------------------------
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_end = text_length - (n_blocks - 1 - i) * max_block_size;
    long block_beg = std::max(0L, block_end - max_block_size);

    threads[i] = new std::thread(gt_end_to_gt_begin_aux,
        text, text_length, block_beg, block_end, gt);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_CHANGE_GT_REFERENCE_POINT_H_INCLUDED
