/**
 * @file    src/psascan_src/bitvector.h
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

#ifndef __PSASCAN_SRC_BITVECTOR_H_INCLUDED
#define __PSASCAN_SRC_BITVECTOR_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <stdint.h>

#include "utils.h"


namespace psascan_private {

struct bitvector {
  private:
    long m_alloc_bytes;
    unsigned char *m_data;

  public:
    bitvector(std::string filename) {
      utils::read_objects_from_file<unsigned char>(m_data, m_alloc_bytes, filename);
    }

    bitvector(long length) {
      m_alloc_bytes = (length + 7) / 8;
      m_data = (unsigned char *)calloc(m_alloc_bytes, sizeof(unsigned char));
    }

    inline bool get(long i) const {
      return m_data[i >> 3] & (1 << (i & 7));
    }

    inline void set(long i) {
      m_data[i >> 3] |= (1 << (i & 7));
    }

    inline void reset(long i) {
      m_data[i >> 3] &= (~(1 << (i & 7)));
    }

    inline void flip(long i) {
      if (get(i)) reset(i);
      else set(i);
    }

    inline void save(std::string filename) const {
      utils::write_objects_to_file<unsigned char>(m_data, m_alloc_bytes, filename);
    }

    // Number of 1 bits in the range [beg..end).
    long range_sum(long beg, long end) const {
      long result = 0L;
    
      long j = beg;
      while (j < end && (j & 63))
        result += get(j++);

      uint64_t *ptr64 = (uint64_t *)(m_data + (j >> 3));
      while (j + 64 <= end) {
        result += __builtin_popcountll(*ptr64++);
        j += 64;
      }

      while (j < end)
        result += get(j++);

      return result;
    }

    ~bitvector() {
      free(m_data);
    }
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_BITVECTOR_H_INCLUDED
