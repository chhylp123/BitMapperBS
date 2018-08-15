/**
 * @file    src/psascan_src/inmem_psascan_src/merge_schedule.h
 * @author  Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_MERGE_SCHEDULE_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_MERGE_SCHEDULE_H_INCLUDED

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>


namespace psascan_private {
namespace inmem_psascan_private {

class MergeSchedule {
private:
  float rl_ratio;
  std::vector<int> split;
  std::vector<int> left_cost;
  std::vector<int> right_cost;

public:
  MergeSchedule(int no_of_blocks, float right_left_ratio,
                int max_left_size = 0)
  { reset(no_of_blocks, right_left_ratio, max_left_size); }

  int left_size(int n) const {
    assert(n < (long)split.size());
    return split[n];
  }
  int right_size(int n) const {
    assert(n < (long)split.size());
    return n - split[n];
  }
  float cost(int n) const {
    assert(n < (long)split.size());
    return (left_cost[n] + rl_ratio * right_cost[n]) / n;
  }
  float n_left_merges(int n) const {
    assert(n < (long)split.size());
    return left_cost[n] / (1.0*n);
  }
  float n_right_merges(int n) const {
    assert(n < (long)split.size());
    return right_cost[n] / (1.0*n);
  }

  void reset(int no_of_blocks, float right_left_ratio,
        int max_left_size = 0)
  {
    int n = no_of_blocks;
    rl_ratio = right_left_ratio;
    if (max_left_size == 0) {
      max_left_size = n-1;
    }

    split.resize(n+1);
    left_cost.resize(n+1);
    right_cost.resize(n+1);
    
    split[1] = 0;
    left_cost[1] = 0;
    right_cost[1] = 0;

    for (int i=2; i<=n; ++i) {
      //int min_l = std::min((i+1)/2, max_left_size);
      int max_l = std::min(i-1, max_left_size);
      float min_cost = 1E40;
      for (int l=1; l<=max_l; ++l) {
        int r = i-l;
        int l_cost = l + left_cost[l] + left_cost[r];
        int r_cost = r + right_cost[l] + right_cost[r];
        float total_cost = l_cost + rl_ratio * r_cost;
        if (total_cost < min_cost) {
          min_cost = total_cost;
          split[i] = l;
          left_cost[i] = l_cost;
          right_cost[i] = r_cost;
        }
      }
    }
  }
};

void print_schedule(const MergeSchedule & sched, int n, std::string indent) {
  if (n == 1) {
    std::cerr << "1\n";
    return;
  }
  std::cerr << n << "\t";
  int l = sched.left_size(n);
  print_schedule(sched, l, indent + ":\t");
  std::cerr << indent;
  print_schedule(sched, n-l, indent + "\t");
}

void print_schedule(const MergeSchedule & sched, int n) {
  std::string intend = "\t";
  print_schedule(sched, n, intend);
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_MERGE_SCHEDULE_H_INCLUDED
