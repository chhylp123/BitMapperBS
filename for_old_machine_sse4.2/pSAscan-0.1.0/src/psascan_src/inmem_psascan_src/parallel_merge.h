/**
 * @file    src/psascan_src/inmem_psascan_src/parallel_merge.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section DESCRIPTION
 *
 * Parallel version of almost in-place stable merging described in
 * the Appending B of
 *
 *   Juha Karkkainen, Peter Sanders, Stefan Burkhardt:
 *   Linear work suffix array construction.
 *   J. ACM 53(6), p. 918-936 (2006).
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_MERGE_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_MERGE_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>

#include "../utils.h"
#include "pagearray.h"
#include "inmem_gap_array.h"


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// Compute the range [res_beg..res_beg+res_size) of the output (i.e., the
// sequence after merging). The range is guaranteed to be aligned with page
// boundaries.
//==============================================================================
template<typename pagearray_type>
void parallel_merge_aux(
    const pagearray_type *l_pagearray,
    const pagearray_type *r_pagearray,
    pagearray_type *output,
    const inmem_gap_array *gap,
    long left_idx, long right_idx,
    long remaining_gap,
    long page_range_beg,
    long page_range_end,
    long what_to_add) {

  typedef typename pagearray_type::value_type value_type;
  static const unsigned pagesize = pagearray_type::pagesize;

  long res_beg = std::max(0L, output->get_page_addr(page_range_beg) - output->m_origin);
  long res_end = std::min(output->m_length, output->get_page_addr(page_range_end) - output->m_origin);
  long res_size = res_end - res_beg;

  long lpage_read = 0L;
  long lpage_id = l_pagearray->get_page_id(left_idx);
  long lpage_offset = l_pagearray->get_page_offset(left_idx);
  value_type *lpage = l_pagearray->m_pageindex[lpage_id++];

  long rpage_read = 0L;
  long rpage_id = r_pagearray->get_page_id(right_idx);
  long rpage_offset = r_pagearray->get_page_offset(right_idx);
  value_type *rpage = r_pagearray->m_pageindex[rpage_id++];

  long pageid = output->get_page_id(res_beg);
  long filled = output->get_page_offset(res_beg);
  value_type *dest = new value_type[pagesize];
  output->m_pageindex[pageid++] = dest;

  std::stack<value_type*> freepages;
  size_t excess_ptr = std::lower_bound(gap->m_excess.begin(),
      gap->m_excess.end(), left_idx + 1) - gap->m_excess.begin();

  for (long i = 0; i < res_size; ++i) {
    if (filled == pagesize) {
      if (freepages.empty()) dest = new value_type[pagesize];
      else { dest = freepages.top(); freepages.pop(); }
      output->m_pageindex[pageid++] = dest;
      filled = 0L;
    }
    if (remaining_gap > 0) {
      --remaining_gap;
      // The next element comes from the right subarray.
      dest[filled] = rpage[rpage_offset++];
      dest[filled++].sa += what_to_add;
      rpage_read++;
      if (rpage_offset == pagesize) {
        // We reached the end of page in the right subarray.
        // We put it into free pages if we read exactly
        // pagesize elements from it. This means the no other
        // thread will attemp to read from it in the future.
        if (rpage_read == pagesize) freepages.push(r_pagearray->m_pageindex[rpage_id - 1]);

        // Note: we don't have to check, if the page below exists, because we have
        // a sentinel page in the page index of every pagearray.
        rpage = r_pagearray->m_pageindex[rpage_id++];
        rpage_offset = 0L;
        rpage_read = 0L;
      }
    } else {
      // Next elem comes from the left subarray.
      dest[filled++] = lpage[lpage_offset++];
      left_idx++;
      lpage_read++;

      // Compute gap[left_idx].
      long gap_left_idx = gap->m_count[left_idx];
      while (excess_ptr < gap->m_excess.size() &&
          gap->m_excess[excess_ptr] == left_idx) {
        gap_left_idx += 256L;
        ++excess_ptr;
      }

      remaining_gap = gap_left_idx;
      if (lpage_offset == pagesize) {
        // We reached the end of page in the left
        // subarray, proceed analogously.
        if (lpage_read == pagesize) freepages.push(l_pagearray->m_pageindex[lpage_id - 1]);

        // Note: we don't have to check, if the page below exists, because we have
        // a sentinel page in the page index of every pagearray.
        lpage = l_pagearray->m_pageindex[lpage_id++];
        lpage_offset = 0L;
        lpage_read = 0L;
      }
    }
  }

  // Release the unused auxiliary pages.
  while (!freepages.empty()) {
    value_type* p = freepages.top();
    freepages.pop();
    if (!output->owns_page(p))
      delete[] p;
  }
}

template<typename pagearray_type>
pagearray_type *parallel_merge(pagearray_type *l_pagearray,
    pagearray_type *r_pagearray, const inmem_gap_array *gap, long max_threads,
    long i0, long &aux_result, long what_to_add) {
  static const unsigned pagesize_log = pagearray_type::pagesize_log;
  static const unsigned pagesize = pagearray_type::pagesize;
  typedef typename pagearray_type::value_type value_type;
  typedef pagearray<value_type, pagesize_log> output_type;

  //----------------------------------------------------------------------------
  // STEP 1: compute the initial parameters for each thread.
  //----------------------------------------------------------------------------
  fprintf(stderr, "queries: ");
  long double start = utils::wclock();
  long length = l_pagearray->m_length + r_pagearray->m_length;
  long n_pages = (length + pagesize - 1) / pagesize;
  long pages_per_thread = (n_pages + max_threads - 1) / max_threads;
  long n_threads = (n_pages + pages_per_thread - 1) / pages_per_thread;
  output_type *result = new output_type(l_pagearray->m_origin, length);

  long *left_idx = new long[n_threads];
  long *right_idx = new long[n_threads];
  long *remaining_gap = new long[n_threads];

  // Prepare gap queries.
  long *gap_query = new long[n_threads];
  long *gap_answer_a = new long[n_threads];
  long *gap_answer_b = new long[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    long page_range_beg = i * pages_per_thread;
    long res_beg = std::max(0L, result->get_page_addr(page_range_beg) - result->m_origin);
    gap_query[i] = res_beg;
  }

  // Answer these queries in parallel and convert the answers
  // to left_idx, right_idx and remaining_gap values.
  aux_result = gap->answer_queries(n_threads, gap_query, gap_answer_a, gap_answer_b, max_threads, i0);
  for (long i = 0; i < n_threads; ++i) {
    long page_range_beg = i * pages_per_thread;
    long res_beg = std::max(0L, result->get_page_addr(page_range_beg) - result->m_origin);
    long j = gap_answer_a[i], s = gap_answer_b[i];
    left_idx[i] = j;
    right_idx[i] = res_beg - j;
    remaining_gap[i] = j + s - res_beg;
  }
  delete[] gap_query;
  delete[] gap_answer_a;
  delete[] gap_answer_b;
  fprintf(stderr, "%.2Lf ", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 2: merge the arrays.
  //----------------------------------------------------------------------------
  fprintf(stderr, "merge: ");
  start = utils::wclock();

  std::thread **threads = new std::thread*[n_threads];
  for (long t = 0; t < n_threads; ++t) {
    long page_range_beg = t * pages_per_thread;
    long page_range_end = std::min(page_range_beg + pages_per_thread, n_pages);

    threads[t] = new std::thread(parallel_merge_aux<pagearray_type>,
        l_pagearray, r_pagearray, result, gap,  left_idx[t], right_idx[t],
        remaining_gap[t], page_range_beg, page_range_end, what_to_add);
  }
  for (long t = 0; t < n_threads; ++t) threads[t]->join();
  for (long t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;
  delete[] left_idx;
  delete[] right_idx;
  delete[] remaining_gap;

  bool *usedpage = new bool[n_pages];
  std::fill(usedpage, usedpage + n_pages, false);

  // Handle the page that was not full
  // manually (if there was one).
  if (length % pagesize) {
    long size = length % pagesize;
    value_type *src = result->m_pageindex[0];
    value_type *dest = result->get_page_addr(0);
    std::copy(src + pagesize - size, src + pagesize, dest + pagesize - size);
    result->m_pageindex[0] = dest;
    usedpage[0] = true;

    // Release the lastpage if it was temporary.
    if (!result->owns_page(src))
      delete[] src;
  }

  // Find unused input pages.
  std::vector<std::pair<long, value_type*> > auxpages;
  for (long i = 0; i < n_pages; ++i) {
    value_type *p = result->m_pageindex[i];
    if (result->owns_page(p)) usedpage[result->get_page_id(p)] = true;
    else auxpages.push_back(std::make_pair(i, p));
  }

  // Assign aux pages to unused pages in any
  // order and release them (aux pages).
  for (long i = 0, ptr = 0; i < n_pages; ++i) {
    if (!usedpage[i]) {
      long id = auxpages[ptr].first;
      value_type *src = auxpages[ptr++].second;
      value_type *dest = result->get_page_addr(i);
      std::copy(src, src + pagesize, dest);
      result->m_pageindex[id] = dest;
      delete[] src;
    }
  }
  delete[] usedpage;
  fprintf(stderr, "%.2Lf ", utils::wclock() - start);

  return result;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_MERGE_H_INCLUDED
