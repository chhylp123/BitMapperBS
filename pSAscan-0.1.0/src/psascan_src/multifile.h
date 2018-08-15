/**
 * @file    src/psascan_src/multifile.h
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

#ifndef __PSASCAN_SRC_MULTIFILE_H_INCLUDED
#define __PSASCAN_SRC_MULTIFILE_H_INCLUDED

#include <vector>
#include <string>

#include "utils.h"


namespace psascan_private {

struct single_file_info {
  long m_beg;
  long m_end;
  std::string m_filename;

  single_file_info(long beg, long end, std::string filename) {
    m_beg = beg;
    m_end = end;
    m_filename = filename;
  }
};

struct multifile {
  std::vector<single_file_info> files_info;

  void add_file(long beg, long end, std::string filename) {
    files_info.push_back(single_file_info(beg, end, filename));
  }

  ~multifile() {
    for (size_t i = 0; i < files_info.size(); ++i)
      utils::file_delete(files_info[i].m_filename);
  }
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_MULTIFILE_H_INCLUDED
