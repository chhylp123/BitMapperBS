/**
 * @file    src/psascan_src/gap_array.h
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

#ifndef __PSASCAN_SRC_GAP_ARRAY_H_INCLUDED
#define __PSASCAN_SRC_GAP_ARRAY_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <mutex>
#include <string>
#include <thread>
#include <algorithm>
#include <parallel/algorithm>

#include "utils.h"
#include "bitvector.h"
#include "parallel_utils.h"
#include "async_stream_writer.h"


namespace psascan_private {

struct buffered_gap_array {
  buffered_gap_array(long length, std::string storage_fname = std::string("")) {
    if (length <= 0L) {
      fprintf(stderr, "\nError: attempting to construct empty gap array.\n");
      std::exit(EXIT_FAILURE);
    }

    m_length = length;
    m_count = (unsigned char *)malloc(m_length);
    std::fill(m_count, m_count + m_length, 0);

    m_excess = new long[k_excess_limit];

    // File used to store excess values.
    m_storage_filename = storage_fname;
    if (!m_storage_filename.length())
      m_storage_filename = ".excess." + utils::random_string_hash();

    m_excess_filled = 0L;
    m_excess_disk = 0L;
    m_sorted_excess = NULL;
    m_sequential_read_initialized = false;
  }

  void add_excess(long x) {
    m_excess[m_excess_filled++] = x;
    if (m_excess_filled == k_excess_limit) {
      m_gap_writing_mutex.lock();
      m_excess_disk += m_excess_filled;
      utils::add_objects_to_file(m_excess, m_excess_filled, m_storage_filename);
      m_excess_filled = 0L;
      m_gap_writing_mutex.unlock();
    }
  }

  void flush_excess_to_disk() {
    if (m_excess_filled > 0) {
      utils::add_objects_to_file(m_excess, m_excess_filled, m_storage_filename);
      m_excess_disk += m_excess_filled;
      m_excess_filled = 0L;
    }
  }

  void start_sequential_access() {
    if (!m_sequential_read_initialized) {
      m_sequential_read_initialized = true;
      m_total_excess = m_excess_filled + m_excess_disk;
      m_sorted_excess = (long *)malloc(m_total_excess * sizeof(long));
      std::copy(m_excess, m_excess + m_excess_filled, m_sorted_excess);
      if (m_excess_disk > 0L) {
        long *dest = m_sorted_excess + m_excess_filled;
        long toread = m_excess_disk;
        utils::read_n_objects_from_file(dest, toread, m_storage_filename.c_str());
      }
      std::sort(m_sorted_excess, m_sorted_excess + m_total_excess);
    }

    m_excess_ptr = 0;
    m_current_pos = 0;
  }

  inline long get_next() {
    long c = 0;
    while (m_excess_ptr < m_total_excess && m_sorted_excess[m_excess_ptr] == m_current_pos)
      ++m_excess_ptr, ++c;
    long result = c * 256L + m_count[m_current_pos];

    ++m_current_pos;
    return result;
  }

  void stop_sequential_access() {
    if (m_sequential_read_initialized) {
      free(m_sorted_excess);
      m_sequential_read_initialized = false;
    } else {
      fprintf(stderr, "\nError: attempting to stop sequential "
          "access to the gap array before it was initialized.\n");
      std::exit(EXIT_FAILURE);
    }
  }

  std::mutex m_excess_mutex;
  std::mutex m_gap_writing_mutex;

  ~buffered_gap_array() {
    if (m_sequential_read_initialized) {
      fprintf(stderr, "\nError: sequential access to gap was not terminated.");
      std::exit(EXIT_FAILURE);
    }

    free(m_count);
    delete[] m_excess;
  }

  void erase_disk_excess() {
    if (utils::file_exists(m_storage_filename))
      utils::file_delete(m_storage_filename);
  }
  
  // Write to a given file using v-byte encoding.
  void save_to_file(std::string fname) {
    fprintf(stderr, "    Write gap to file: ");
    long double gap_write_start = utils::wclock();
    long bytes_written = 0L;

    start_sequential_access();
    typedef async_stream_writer<unsigned char> stream_writer_type;
    stream_writer_type *writer = new stream_writer_type(fname);

    for (long j = 0; j < m_length; ++j) {
      long val = get_next();
      while (val > 127) {
        writer->write((val & 0x7f) | 0x80);
        val >>= 7;
        ++bytes_written;
      }
      writer->write(val);
    }

    bytes_written += m_length;
    stop_sequential_access();
    delete writer;

    long double gap_write_time = utils::wclock() - gap_write_start;
    long double io_speed = (bytes_written / (1024.L * 1024)) / gap_write_time;
    fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", gap_write_time, io_speed);
  }
  
  
  //==============================================================================
  // Note about the input:
  // - j is the maximal integer such that gapsum[j] + j <= beg.
  // - S contains value gapsum[j] + j.
  //==============================================================================
  static void convert_gap_to_bitvector_aux(long beg, long end, long j, long S, buffered_gap_array *gap, bitvector *bv) {
    // Initialize pointer to sorted excess values.
    long excess_pointer = std::lower_bound(gap->m_sorted_excess,
        gap->m_sorted_excess + gap->m_total_excess, j) - gap->m_sorted_excess;

    // Compute gap[j].
    long gap_j = gap->m_count[j];
    while (excess_pointer < gap->m_total_excess && gap->m_sorted_excess[excess_pointer] == j) {
      gap_j += 256L;
      ++excess_pointer;
    }

    long p = beg;
    long ones = std::min(end - p, gap_j - (beg - S));
    for (long k = 0; k < ones; ++k) bv->set(p++);
    ++j;

    while (p < end) {
      ++p;

      // Compute gap[j].
      gap_j = gap->m_count[j];
      while (excess_pointer < gap->m_total_excess && gap->m_sorted_excess[excess_pointer] == j) {
        gap_j += 256L;
        ++excess_pointer;
      }

      ones = std::min(end - p, gap_j);

      for (long k = 0; k < ones; ++k) bv->set(p++);
      ++j;
    }
  }

  static void compute_j_aux(long range_beg, long n_chunks, long max_chunk_size,
      const long *sparse_gapsum, long &initial_gap_ptr, long &initial_gapsum_value, const buffered_gap_array *gap) {
    // Fast forward through as many chunks as possible.
    long j = 0L;
    long gapsum_j = 0L;  // At any time gapsum_j = gap[0] + .. + gap[j - 1].
    while (j + 1 < n_chunks && sparse_gapsum[j + 1] + (max_chunk_size * (j + 1)) <= range_beg) ++j;
    gapsum_j = sparse_gapsum[j];
    j = (j * max_chunk_size);

    // Slowly find the right place in a single chunk.
    long excess_ptr = std::lower_bound(gap->m_sorted_excess, gap->m_sorted_excess + gap->m_total_excess, j) - gap->m_sorted_excess;
    while (j < gap->m_length) {
      long gap_j = gap->m_count[j];
      while (excess_ptr < gap->m_total_excess && gap->m_sorted_excess[excess_ptr] == j) {
        gap_j += 256L;
        ++excess_ptr;
      }

      if (gapsum_j + gap_j + j + 1 <= range_beg) {
        gapsum_j += gap_j;
        ++j;
      } else break;
    }

    // Store the answer.
    initial_gap_ptr = j;
    initial_gapsum_value = gapsum_j + j;
  }


  static void compute_gapsum_for_chunk_group(long group_beg, long group_end, long max_chunk_size,
      long *sparse_gapsum, const buffered_gap_array *gap) {
    for (long chunk_id = group_beg; chunk_id < group_end; ++chunk_id) {
      long chunk_beg = chunk_id * max_chunk_size;
      long chunk_end = std::min(chunk_beg + max_chunk_size, gap->m_length);

      // Compute sum of gap values inside chunk. We assume that
      // the excess values are in RAM and were sorted.
      long occ = std::upper_bound(gap->m_sorted_excess, gap->m_sorted_excess + gap->m_total_excess, chunk_end - 1)
        - std::lower_bound(gap->m_sorted_excess, gap->m_sorted_excess + gap->m_total_excess, chunk_beg);
      long gap_sum_inside_chunk = 256L * std::max(0L, occ);
      for (long j = chunk_beg; j < chunk_end; ++j)
        gap_sum_inside_chunk += gap->m_count[j];

      // Store the result.
      sparse_gapsum[chunk_id] = gap_sum_inside_chunk;
    }
  }

  bitvector* convert_to_bitvector(long max_threads) {
    // 1
    //
    // The term chunks is used to compute sparse gapsum array.
    // Chunk is a length such that
    // gapsum[k] = gap[0] + gap[1] + .. + gap[k * max_chunk_size - 1]
    long max_chunk_size = std::min(4L << 20, (m_length + max_threads - 1) / max_threads);
    long n_chunks = (m_length + max_chunk_size - 1) / max_chunk_size;
    long *sparse_gapsum = (long *)malloc(n_chunks * sizeof(long));


    // 2
    //
    // Compute the sum of gap value inside each chunk. Since there can be
    // more chunks than threads, we split chunks into groups and let each
    // thread compute the sum of gap values inside the group of chunks.
    long chunk_group_size = (n_chunks + max_threads - 1) / max_threads;
    long n_chunk_groups = (n_chunks + chunk_group_size - 1) / chunk_group_size;

    start_sequential_access();
    std::thread **threads = new std::thread*[n_chunk_groups];
    for (long t = 0; t < n_chunk_groups; ++t) {
      long chunk_group_beg = t * chunk_group_size;
      long chunk_group_end = std::min(chunk_group_beg + chunk_group_size, n_chunks);

      threads[t] = new std::thread(compute_gapsum_for_chunk_group, chunk_group_beg,
          chunk_group_end, max_chunk_size, sparse_gapsum, this);
    }

    for (long t = 0; t < n_chunk_groups; ++t) threads[t]->join();
    for (long t = 0; t < n_chunk_groups; ++t) delete threads[t];
    delete[] threads;


    // 3
    //
    // Compute comulative sum over sparse_gapsum array.
    long double gap_total_sum = 0L;
    for (long i = 0L; i < n_chunks; ++i) {
      long temp = sparse_gapsum[i];
      sparse_gapsum[i] = gap_total_sum;
      gap_total_sum += temp;
    }


    // 4
    //
    // Compute all initial gap pointers. For a thread handling range [beg..end), the
    // initial_gap_ptr values is the largest j, such that gapsum[j] + j <= beg.
    // After we find j, we store the value of gapsum[j] + j in initial_gapsum_value.
    long result_length = (m_length + gap_total_sum) - 1;
    bitvector *result = new bitvector(result_length + 1);  // +1 is to make room for sentinel

    long max_range_size = (result_length + max_threads - 1) / max_threads;
    while (max_range_size & 7) ++max_range_size;
    long n_ranges = (result_length + max_range_size - 1) / max_range_size;

    long *initial_gap_ptr = new long[n_ranges];
    long *initial_gapsum_value = new long[n_ranges];

    threads = new std::thread*[n_ranges];
    for (long t = 0; t < n_ranges; ++t) {
      long range_beg = t * max_range_size;
      threads[t] = new std::thread(compute_j_aux, range_beg, n_chunks, max_chunk_size,
          sparse_gapsum, std::ref(initial_gap_ptr[t]), std::ref(initial_gapsum_value[t]), this);
    }
    for (long t = 0; t < n_ranges; ++t) threads[t]->join();
    for (long t = 0; t < n_ranges; ++t) delete threads[t];


    // 5
    //
    // Compute the bitvector. Each thread fills in the range of bits.
    for (long t = 0; t < n_ranges; ++t) {
      long range_beg = t * max_range_size;
      long range_end = std::min(range_beg + max_range_size, result_length);

      threads[t] = new std::thread(convert_gap_to_bitvector_aux, range_beg,
          range_end, initial_gap_ptr[t], initial_gapsum_value[t], this, result);
    }

    for (long t = 0; t < n_ranges; ++t) threads[t]->join();
    for (long t = 0; t < n_ranges; ++t) delete threads[t];
    delete[] threads;

    delete[] initial_gap_ptr;
    delete[] initial_gapsum_value;
    stop_sequential_access();
    free(sparse_gapsum);

    return result;
  }
  
  static const long k_excess_limit = (1L << 22);

  unsigned char *m_count;
  long m_length;
  long m_excess_filled;
  long m_excess_disk;
  long *m_excess;

  std::string m_storage_filename;

  bool m_sequential_read_initialized;
  long m_excess_ptr;
  long m_current_pos;

public:
  long *m_sorted_excess;
  long m_total_excess; 
};


struct gap_array_2n {
  gap_array_2n(const buffered_gap_array *gap, long max_threads) {
    m_length = gap->m_length;
    m_count = (uint16_t *)malloc(m_length * sizeof(uint16_t));
    parallel_utils::parallel_copy<unsigned char, uint16_t>(gap->m_count, m_count, m_length, max_threads);
    m_storage_filename = gap->m_storage_filename;
    m_excess_disk = gap->m_excess_disk;
  }

  gap_array_2n(long length) {
    m_length = length;
    m_count = (uint16_t *)malloc(m_length * sizeof(uint16_t));
  }

  ~gap_array_2n() {
    if (m_count)
      free(m_count);
  }

  static void apply_excess_aux(gap_array_2n *gap, const long *tab,
      long block_beg, long block_end, uint64_t &initial_run_length) {
    long block_size = block_end - block_beg;

    // Each thread gathers excess values in a buffer and at the end
    // copies then to the gap array's mutex-protected m_excess vector.
    std::vector<long> excess_buffer;

    // Compute the length of initial run.
    initial_run_length = 1UL;
    while (initial_run_length < (uint64_t)block_size && tab[block_beg] ==
        tab[block_beg + initial_run_length]) ++initial_run_length;

    // Update count values.
    for (long i = block_beg + initial_run_length; i < block_end; ++i) {
      long x = tab[i];
      uint64_t value = (uint64_t)gap->m_count[x] + 256UL;
      if (value >= (1UL << 16)) {
        value -= (1UL << 16);
        excess_buffer.push_back(x);
      }
      gap->m_count[x] = value;
    }

    // Copy the excess values to the gap array's mutex-protected vector.
    std::unique_lock<std::mutex> lk(gap->m_excess_mutex);
    for (long i = 0; i < (long)excess_buffer.size(); ++i)
      gap->m_excess.push_back(excess_buffer[i]);
    lk.unlock();
  }

  void apply_excess_from_disk(long ram_budget, long max_threads) {
    if (!m_excess_disk) return;

    // We only use half of the RAM for buffer, because we will use parallel
    // merge sort for sorting the buffer (which requires double the space
    // for the input).
    long elems = std::max(1L, ram_budget / (2L * (long)sizeof(long)));
    long *buffer = (long *)malloc(elems * sizeof(long));

    std::FILE *f = utils::open_file(m_storage_filename.c_str(), "r");
    std::thread **threads = new std::thread*[max_threads];

    // After sorting the buffer, when we split it equally between threads
    // we obey the rule, the every thread only counts the number of 
    // elements equal to the first element in the handled range, but does
    // not do any updates for these elements. This prevents two threads
    // trying to update the same elements in the m_count array. The
    // length of the first run is computed and returned by each thread.
    // It is then updated sequentially.
    uint64_t *first_run_length = new uint64_t[max_threads];
 
    while (m_excess_disk > 0) {
      // Read a portion of excess values from disk.
      long toread = std::min(m_excess_disk, elems);
      utils::read_n_objects_from_file(buffer, toread, f);

      // Sort excess values in parallel.
      __gnu_parallel::sort(buffer, buffer + toread);

      // Update m_count and m_excess with elements from the buffer.
      // The buffer is dividied into blocks, each blocks handles one
      // block. Each thread updates the values except the first run
      // of the block, which is handled separatelly (sequentially).
      long max_block_size = (toread + max_threads - 1) / max_threads;
      long n_blocks = (toread + max_block_size - 1) / max_block_size;

      for (long t = 0; t < n_blocks; ++t) {
        long block_beg = t * max_block_size;
        long block_end = std::min(block_beg + max_block_size, toread);

        threads[t] = new std::thread(apply_excess_aux, this, buffer,
            block_beg, block_end, std::ref(first_run_length[t]));
      }

      for (long t = 0; t < n_blocks; ++t) threads[t]->join();
      for (long t = 0; t < n_blocks; ++t) delete threads[t];

      // Sequentially handle the elements in the first run of each block.
      for (long t = 0; t < n_blocks; ++t) {
        long block_beg = t * max_block_size;
        long first = buffer[block_beg];  // first elements in the block

        uint64_t freq = (uint64_t)m_count[first] + (first_run_length[t] * 256L);
        while (freq >= (1UL << 16)) {
          freq -= (1UL << 16);
          m_excess.push_back(first);
        }
        m_count[first] = freq;
      }

      m_excess_disk -= toread;
    }

    __gnu_parallel::sort(m_excess.begin(), m_excess.end());

    delete[] threads;
    delete[] first_run_length;

    std::fclose(f);
    free(buffer);
  }

  void set_count(long pos, long value) {
    while (value >= (1L << 16)) {
      m_excess.push_back(pos);
      value -= (1L << 16);
    }
    m_count[pos] = (uint64_t)value;
  }

  void erase_disk_excess() {
    if (utils::file_exists(m_storage_filename))
      utils::file_delete(m_storage_filename);
  }

  uint16_t *m_count;

  long m_length;
  long m_excess_disk;

  std::mutex m_excess_mutex;
  std::string m_storage_filename;
  std::vector<long> m_excess;  // all excess values are in RAM
};

}  // psascan_private

#endif  // __PSASCAN_SRC_GAP_ARRAY_H_INCLUDED
