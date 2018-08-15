/**
 * @file    src/psascan_src/em_compute_initial_ranks.h
 * @author  Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *          Dominik Kempa <dominik.kempa (at) gmail.com>
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

#ifndef __PSASCAN_SRC_EM_COMPUTE_INITIAL_RANKS_INCLUDED
#define __PSASCAN_SRC_EM_COMPUTE_INITIAL_RANKS_INCLUDED

#include <string>
#include <vector>
#include <algorithm>
#include <thread>

#include "approx_rank.h"
#include "sparse_isa.h"
#include "background_block_reader.h"
#include "background_chunk_reader.h"
#include "multifile_bit_stream_reader.h"
#include "utils.h"


namespace psascan_private {

// #define EM_STARTING_POS_MODULE_DEBUG_MODE

inline int lcp_compare(
    const unsigned char *text,  // only text[block_suf_beg..block_end) will be accessed
    long text_length,
    long block_end,             // wrt to text beg
    long block_suf_beg,         // wrt to text beg
    const unsigned char *pat,   // only pat[lcp..pat_length) will be accessed
    long pat_beg,               // wrt to text beg
    long pat_length,
    multifile_bit_stream_reader &gt_reader,
    long &lcp) {
  while (block_suf_beg + lcp < block_end && lcp < pat_length &&
      text[block_suf_beg + lcp] == pat[lcp]) ++lcp;
  if (block_suf_beg + lcp >= block_end) {
    if (gt_reader.access(text_length - (pat_beg + (block_end - block_suf_beg)))) return 1;
    else return -1;
  } else if (lcp == pat_length) {
    if (pat_beg + pat_length >= text_length) return -1;
    else return 0;
  } else {
    if (pat[lcp] > text[block_suf_beg + lcp]) return 1;
    else return -1;
  } 
}

template<typename saidx_t>
void refine_range(
    const unsigned char *block,
    const saidx_t *block_psa,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long pat_beg,    // same here
    long text_length,
    long left,
    long right,
    long old_lcp,
    long new_lcp,
    const unsigned char *pat,  // only pat[old_lcp..new_lcp) can and will be accessed
    multifile_bit_stream_reader &gt_reader,
    long &newleft,
    long &newright) {
  long low = left - 1;
  long high = right;
  long llcp = old_lcp;
  long rlcp = old_lcp;

#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  long min_discrepancy = utils::random_long(0L, 10L);
  long balancing_factor = utils::random_long(1L, 10L);
#else
  static const long min_discrepancy = (1L << 16);
  static const long balancing_factor = 64L;
#endif

  const unsigned char *text = block - block_beg;
  while (low + 1 != high) {
    // Invariant: newleft is in the range (low, high].
    long lcp = std::min(llcp, rlcp);
    long mid = 0L;
    if (llcp + min_discrepancy < rlcp) {
      long d = rlcp - llcp;
      long logd = utils::log2ceil(d);
      mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      long d = llcp - rlcp;
      long logd = utils::log2ceil(d);
      mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else mid = (low + high) / 2;

    if (lcp_compare(text, text_length, block_end, block_beg + (long)block_psa[mid],
        pat, pat_beg, new_lcp, gt_reader, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }
  newleft = high;

  if (rlcp >= new_lcp) {
    high = right;
    rlcp = old_lcp;

    while (low + 1 != high) {
      // Invariant: newright is in the range (low, high].
      long lcp = std::min(llcp, rlcp);
      long mid = 0L;
      if (llcp + min_discrepancy < rlcp) {
        long d = rlcp - llcp;
        long logd = utils::log2ceil(d);
        mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        long d = llcp - rlcp;
        long logd = utils::log2ceil(d);
        mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare(text, text_length, block_end, block_beg + (long)block_psa[mid],
            pat, pat_beg, new_lcp, gt_reader, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }
  }
  newright = high;
}

template<typename saidx_t>
void em_compute_single_initial_rank(
    const unsigned char *block,
    const saidx_t *block_psa,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long pat_beg,    // same here
    long text_length,
    long max_lcp,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::pair<long, long> &result) {
  if (pat_beg == text_length) {
    result = std::make_pair(0, 0);
    return;
  }

  long block_size = block_end - block_beg;
  long pat_end = pat_beg + max_lcp;

  multifile_bit_stream_reader gt_reader(tail_gt_begin_reversed);

  // Reads text[pat_beg..pat_end) in chunks.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  long chunk_length = utils::random_long(1L, 10L); 
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end, chunk_length);
#else
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end);
#endif

  // The current range is [left, right).
  long left = 0;
  long right = block_size;
  long lcp = 0;

  while (left != right && lcp < max_lcp) {
    long this_chunk_length = std::min(max_lcp - lcp, chunk_reader->get_chunk_size());
    long new_lcp = lcp + this_chunk_length;
    chunk_reader->wait(pat_beg + new_lcp);

    // Invariant:
    //   reader->chunk[0..chunk_length) = pattern[lcp..new_lcp).
    long newleft = 0;
    long newright = 0;
    refine_range(block, block_psa, block_beg, block_end, pat_beg, text_length, left,
        right, lcp, new_lcp, chunk_reader->m_chunk - lcp, gt_reader, newleft, newright);
    left = newleft;
    right = newright;
    lcp = new_lcp;
  }

  delete chunk_reader;

  result = std::make_pair(left, right);
}

template<typename saidx_t>
void em_compute_initial_ranks(
    const unsigned char *block,
    const saidx_t *block_psa,
    const unsigned char *block_pbwt,
    long i0,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long text_length,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::vector<long> &result,
    long max_threads,
    long tail_end,
    long initial_rank_after_tail) {
  // Note, that bits of tail_gt_begin_reversed are indexed in the
  // range [text_length - tail_end.. text_length - block_end). This
  // is because the same multifile is then used in the streaming and
  // for streaming is much more natural to use this indexing.
  long block_length = block_end - block_beg;
  long tail_length = tail_end - block_end;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  std::vector<std::pair<long, long> > ranges(n_threads);
  std::thread **threads = new std::thread*[n_threads];

  for (int t = n_threads - 1; t >= 0; --t) {
    long stream_block_beg = block_end + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, tail_end);
    long stream_block_size = stream_block_end - stream_block_beg;

    threads[t] = new std::thread(em_compute_single_initial_rank<saidx_t>,
        block, block_psa, block_beg, block_end, stream_block_beg, text_length,
        stream_block_size, text_filename, tail_gt_begin_reversed, std::ref(ranges[t]));
  }

  for (int t = 0; t < n_threads; ++t) threads[t]->join();
  for (int t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;

  // Refine ranges until all are single elements.
  result.resize(n_threads);

  bool nontrivial_range = false;
  for (long t = 0; t < n_threads; ++t)
    if (ranges[t].first != ranges[t].second)
      nontrivial_range = true;

  if (nontrivial_range) {
    multifile_bit_stream_reader *gt_reader =
      new multifile_bit_stream_reader(tail_gt_begin_reversed);

#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
    typedef approx_rank<1L> rank_type;
    typedef sparse_isa<rank_type, saidx_t, 1L> isa_type;
#else
    typedef approx_rank<8L> rank_type;
    typedef sparse_isa<rank_type, saidx_t, 8L> isa_type;
#endif
    rank_type *pbwt_rank = new rank_type(block_pbwt, block_length, max_threads);
    isa_type *block_sparse_isa = new isa_type(block_psa, block, block_length, i0, pbwt_rank, max_threads);

    long prev_rank = initial_rank_after_tail;
    for (long t = n_threads - 1; t >= 0; --t) {
      long stream_block_beg = block_end + t * stream_max_block_size;
      long stream_block_end = std::min(stream_block_beg + stream_max_block_size, tail_end);
      long stream_block_size = stream_block_end - stream_block_beg;

      long left = ranges[t].first;
      long right = ranges[t].second;

      while (left != right) {
        // Valid values for mid are in [left..right).
        long mid = (left + right) / 2;

        if ((long)block_psa[mid] + stream_block_size >= block_length) {
          if (gt_reader->access(text_length - (stream_block_beg + (block_length - (long)block_psa[mid])))) left = mid + 1;
          else right = mid; 
        } else {
          long j = (long)block_psa[mid] + stream_block_size;
          if (block_sparse_isa->query(j) < prev_rank) left = mid + 1;
          else right = mid;
        }
      }

      result[t] = left;
      prev_rank = result[t];
    }

    delete pbwt_rank;
    delete block_sparse_isa;
    delete gt_reader;
  } else {
    for (long t = 0; t < n_threads; ++t)
      result[t] = ranges[t].first;
  }
}

int lcp_compare_2(
    const unsigned char *text,  // only text[block_suf_beg..block_end) will be accessed
    long text_length,
    long block_end,             // wrt to text beg
    long block_suf_beg,         // wrt to text beg
    const unsigned char *pat,   // only pat[lcp..pat_length) will be accessed
    long pat_beg,               // wrt to text beg
    long pat_length,
    long tail_begin,            // wrt to text beg
    background_block_reader *mid_block_reader,
    multifile_bit_stream_reader &gt_reader,
    long &lcp) {
  while (block_suf_beg + lcp < block_end && lcp < pat_length &&
      text[block_suf_beg + lcp] == pat[lcp]) ++lcp;
  if (block_suf_beg + lcp < block_end && lcp < pat_length) {
    if (pat[lcp] > text[block_suf_beg + lcp]) return 1;
    else return -1;
  }

  if (block_suf_beg + lcp >= block_end && block_end < tail_begin && lcp < pat_length) {
    // To finish the comparison, we need to access symbols from the mid block.
    // First, wait until enough symbols are available.
    mid_block_reader->wait(std::min(tail_begin, block_suf_beg + pat_length) - block_end);

    // Now continue the comparison.
    const unsigned char *text2 = mid_block_reader->m_data - block_end;
    while (block_suf_beg + lcp < tail_begin && lcp < pat_length &&
        text2[block_suf_beg + lcp] == pat[lcp]) ++lcp;
    if (block_suf_beg + lcp < tail_begin && lcp < pat_length) {
      if (pat[lcp] > text2[block_suf_beg + lcp]) return 1;
      else return -1;
    }
  }

  if (block_suf_beg + lcp >= tail_begin) {
    // Use gt to resolve comparison.
    if (gt_reader.access(text_length -  (pat_beg + (tail_begin - block_suf_beg)))) return 1;
    else return -1;
  } else {  // lcp == pat_length
    if (pat_beg + pat_length >= text_length) return -1;
    else return 0;
  }
}

template<typename saidx_t>
void refine_range_2(
    const unsigned char *block,
    const saidx_t *block_psa,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long pat_beg,    // same here
    long tail_begin,
    background_block_reader *mid_block_reader,
    long text_length,
    long left,
    long right,
    long old_lcp,
    long new_lcp,
    const unsigned char *pat,  // only pat[old_lcp..new_lcp) can and will be accessed
    multifile_bit_stream_reader &gt_reader,
    long &newleft,
    long &newright) {
  long low = left - 1;
  long high = right;
  long llcp = old_lcp;
  long rlcp = old_lcp;

#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  long min_discrepancy = utils::random_long(0L, 10L);
  long balancing_factor = utils::random_long(1L, 10L);
#else
  static const long min_discrepancy = (1L << 16);
  static const long balancing_factor = 64L;
#endif

  const unsigned char *text = block - block_beg;
  while (low + 1 != high) {
    // Invariant: newleft is in the range (low, high].
    long lcp = std::min(llcp, rlcp);
    long mid = 0L;
    if (llcp + min_discrepancy < rlcp) {
      long d = rlcp - llcp;
      long logd = utils::log2ceil(d);
      mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      long d = llcp - rlcp;
      long logd = utils::log2ceil(d);
      mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else mid = (low + high) / 2;

    if (lcp_compare_2(text, text_length, block_end, block_beg + (long)block_psa[mid],
        pat, pat_beg, new_lcp, tail_begin, mid_block_reader, gt_reader, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }
  newleft = high;

  if (rlcp >= new_lcp) {
    high = right;
    rlcp = old_lcp;

    while (low + 1 != high) {
      // Invariant: newright is in the range (low, high].
      long lcp = std::min(llcp, rlcp);
      long mid = 0L;
      if (llcp + min_discrepancy < rlcp) {
        long d = rlcp - llcp;
        long logd = utils::log2ceil(d);
        mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        long d = llcp - rlcp;
        long logd = utils::log2ceil(d);
        mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare_2(text, text_length, block_end, block_beg + (long)block_psa[mid],
        pat, pat_beg, new_lcp, tail_begin, mid_block_reader, gt_reader, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }
  }
  newright = high;
}

template<typename saidx_t>
void em_compute_single_initial_rank_2(
    const unsigned char *block,
    const saidx_t *block_psa,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long pat_beg,    // same here
    long text_length,
    long max_lcp,
    long tail_begin,
    background_block_reader *mid_block_reader,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    long &result) {
  if (pat_beg == text_length) {
    result = 0;
    return;
  }

  long block_size = block_end - block_beg;
  long pat_end = std::min(text_length, pat_beg + max_lcp);

  multifile_bit_stream_reader gt_reader(tail_gt_begin_reversed);

  // Reads text[pat_beg..pat_end) in chunks.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  long chunk_length = utils::random_long(1L, 10L); 
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end, chunk_length);
#else
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end);
#endif

  // The current range is [left, right).
  long left = 0;
  long right = block_size;
  long lcp = 0;

  while (left != right && lcp < max_lcp) {
    long this_chunk_length = std::min(max_lcp - lcp, chunk_reader->get_chunk_size());
    long new_lcp = lcp + this_chunk_length;
    chunk_reader->wait(pat_beg + new_lcp);

    // Invariant:
    //   reader->chunk[0..chunk_length) = pattern[lcp..new_lcp).
    long newleft = 0;
    long newright = 0;
    refine_range_2(block, block_psa, block_beg, block_end, pat_beg, tail_begin,
        mid_block_reader, text_length, left, right, lcp, new_lcp,
        chunk_reader->m_chunk - lcp, gt_reader, newleft, newright);
    left = newleft;
    right = newright;
    lcp = new_lcp;
  }
  result = left;

  delete chunk_reader;
}

template<typename saidx_t>
void em_compute_initial_ranks(
    const unsigned char *block,
    const saidx_t *block_psa,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long text_length,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::vector<long> &result,
    long max_threads,
    long tail_begin) {
  // Compute some initial parameters.
  long block_length = block_end - block_beg;
  long tail_length = text_length - tail_begin;
  long mid_block_beg = block_end;
  long mid_block_end = tail_begin;
  long mid_block_size = mid_block_end - mid_block_beg;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  // Start reading the text between the block and the tail in the backgrond.
  background_block_reader *mid_block_reader =
    new background_block_reader(text_filename, mid_block_beg, mid_block_size);

  // Compute the initial ranks.
  std::vector<long> res(n_threads);
  std::thread **threads = new std::thread*[n_threads];

  for (int t = 0; t < n_threads; ++t) {
    long stream_block_beg = tail_begin + t * stream_max_block_size;
    long max_lcp = std::min(block_length + mid_block_size, text_length - stream_block_beg);

    threads[t] = new std::thread(em_compute_single_initial_rank_2<saidx_t>,
        block, block_psa, block_beg, block_end, stream_block_beg, text_length,
        max_lcp, tail_begin, mid_block_reader, text_filename,
        tail_gt_begin_reversed, std::ref(res[t]));
  }

  for (int t = 0; t < n_threads; ++t) threads[t]->join();
  for (int t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;

  mid_block_reader->stop();
  delete mid_block_reader;

  result = res;
}

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_EM_COMPUTE_INITIAL_RANKS_INCLUDED
