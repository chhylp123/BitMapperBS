/**
 * @file    tools/delete-bytes-255/main.cpp
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

#include <cstdio>
#include <cstdlib>
#include <sys/time.h>


long double wallclock() {
  timeval tim;
  gettimeofday(&tim, NULL);
  return tim.tv_sec + (tim.tv_usec / 1000000.L);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr, "Usage: %s FILE\nErase all bytes with value 255 "
        "from FILE. Write result on standard output.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  // Open the input file.
  std::FILE *f = std::fopen(argv[1], "r");
  if (f == NULL) {
    std::perror(argv[1]);
    std::exit(EXIT_FAILURE);
  }

  // Get the file size.
  std::fseek(f, 0L, SEEK_END);
  long size = std::ftell(f);
  std::rewind(f);

  // Allocate the buffer.
  static const long bufsize = (2L << 20);
  unsigned char *buffer = new unsigned char[bufsize];

  // Do the filtering.
  long double start = wallclock();
  std::size_t elems, count = 0, total = 0;
  while ((elems = std::fread(buffer, 1, bufsize, f)) > 0) {
    total += elems;
    count += elems;

    // Filter the buffer.
    std::size_t ptr = 0;
    for (std::size_t j = 0; j < elems; ++j)
      if (buffer[j] != 255)
        buffer[ptr++] = buffer[j];

    // Write filtered buffer to stdout.
    if (ptr > 0)
      std::fwrite(buffer, 1, ptr, stdout);

    // Print progress message.
    if (count > (64L << 20)) {
      count = 0;
      long double elapsed = wallclock() - start;
      long double mib = (long double)total / (1L << 20);
      std::fprintf(stderr, "Processed %.0LfMiB (%.1Lf%%). Speed: %.2LfMiB/s.\r",
          mib, (100.L * total) / size, mib / elapsed);
    }
  }

  // Clean up.
  delete[] buffer;
  std::fclose(f);

  // Print summary.
  long double elapsed = wallclock() - start;
  long double mib = (long double)size / (1L << 20);
  std::fprintf(stderr, "Processed %.0LfMiB (100.0%%). Speed: %.2LfMiB/s.\n",
      mib, mib / elapsed);
}
