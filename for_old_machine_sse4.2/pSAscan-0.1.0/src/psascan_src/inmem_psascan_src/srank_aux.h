/**
 * @file    src/psascan_src/inmem_psascan_src/srank_aux.h
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_H_INCLUDED


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// Compute ms-decomposition of text[0..length) from ms-decomposition of
// text[0..length - 1). The result is returned via updated values s, p, r.
//==============================================================================
template<typename T>
inline void update_ms(const unsigned char *text, T length, T &s, T &p) {
  if (length == 1) { s = 0; p = 1; return; }

  T i = length - 1;
  while (i < length) {
    unsigned char a = text[i - p];
    unsigned char b = text[i];

    if (a > b) p = i - s + 1;
    else if (a < b) {
      long r = (i - s);
      while (r >= p) r -= p;
      i -= r;
      s = i;
      p = 1;
    }

    ++i;
  }
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_H_INCLUDED
