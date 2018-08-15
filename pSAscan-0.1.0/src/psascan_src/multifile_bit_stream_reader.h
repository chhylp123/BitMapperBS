/**
 * @file    src/psascan_src/multifile_bit_stream_reader.h
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

#ifndef __PSASCAN_SRC_MULTIFILE_BIT_STREAM_READER_H_INCLUDED
#define __PSASCAN_SRC_MULTIFILE_BIT_STREAM_READER_H_INCLUDED

#include <cstdio>
#include <vector>

#include "utils.h"
#include "multifile.h"


namespace psascan_private {

struct multifile_bit_stream_reader {
private:
  static const long k_bufsize;

  // info for currently accessed file.
  std::FILE *m_file;
  long m_file_beg;
  long m_file_end;

  unsigned char *m_buffer;
  long m_offset; // number of the first bit in the buffer
  long m_filled; // how many bits we have in a buffer

  long cur_bit_absolute;
  long cur_bit_buffer;
  long cur_bit;
  long cur_byte;

  std::vector<single_file_info> files_info;

public:
  multifile_bit_stream_reader(const multifile *m) {
    m_file = NULL;
    m_file_beg = 0;
    m_file_end = 0;
    m_buffer = new unsigned char[k_bufsize];

    if (m != NULL)
      files_info = m->files_info;
  }

  // Subsequent access operations are quaranteed
  // to be with increasing argument.
  bool access(long i) {
    if (i < m_file_beg || m_file_end <= i) {
      open_file_for_index(i);
      i -= m_file_beg;
    } else {
      i -= m_file_beg;

      if (i < m_offset || m_offset + m_filled <= i) {
        refill(i);
      }
    }

    i -= m_offset;
    return (m_buffer[i >> 3] & (1 << (i & 7)));
  }

  void initialize_sequential_reading(long i) {
    open_file_for_index(i);

    cur_bit_absolute = i;
    cur_bit_buffer = cur_bit_absolute - (m_file_beg + m_offset);
    cur_byte = (cur_bit_buffer >> 3);
    cur_bit = (cur_bit_buffer & 7);
  }

  inline bool read() {
    if (cur_bit_absolute == m_file_end) open_file_for_index(m_file_end);
    if (cur_bit_buffer == m_filled) refill(m_offset + m_filled);

    bool ans = (m_buffer[cur_byte] & (1 << cur_bit));
    ++cur_bit;
    if (cur_bit == 8) {
      cur_bit = 0;
      ++cur_byte;
    }
    
    ++cur_bit_buffer;
    ++cur_bit_absolute;
    return ans;
  }

  ~multifile_bit_stream_reader() {
    if (m_file)
      std::fclose(m_file);
    delete[] m_buffer;
  }

private:
  void refill(long offset) {
    offset -= (offset & 7);
    if (m_offset + m_filled != offset)
      std::fseek(m_file, (offset >> 3), SEEK_SET);
    long bytes_read = std::fread(m_buffer, 1, k_bufsize, m_file);
    m_filled = std::min(m_file_end - offset, 8L * bytes_read); // in bits
    m_offset = offset; // in bits

    cur_byte = 0; // in the buffer
    cur_bit = 0; // in the current byte
    cur_bit_buffer = 0;
  }

  void open_file_for_index(long i) {
    // Close current file (if any is open).
    if (m_file) std::fclose(m_file);

    // First find the right file.
    long id = 0;
    while (i < files_info[id].m_beg || files_info[id].m_end <= i)
      ++id;

    m_file = utils::open_file(files_info[id].m_filename, "r");
    m_file_beg = files_info[id].m_beg;
    m_file_end = files_info[id].m_end;

    cur_bit_absolute = m_file_beg;
    cur_bit_buffer = 0;
    cur_bit = 0;
    cur_byte = 0;

    m_offset = 0;
    m_filled = 0;

    refill(i - m_file_beg);
  }
};

const long multifile_bit_stream_reader::k_bufsize = (1L << 20);

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_MULTIFILE_BIT_STREAM_READER_H_INCLUDED
