pSAscan - parallel external memory suffix array construction algorithm.
=======================================================================


Description
-----------

This package contains implementation of the parallel external-memory
suffix array construction algorithm called pSAscan described in the paper

  Juha Karkkainen, Dominik Kempa, and Simon J. Puglisi,
  Parallel External Memory Suffix Sorting.
  In Proc. 26th Annual Symposium on Combinatorial Pattern Matching (CPM) 2015.

The algorithm is based on the sequential external-memory suffix array
construction algorithm called SAscan described in

  Juha Karkkainen and Dominik Kempa,
  Engineering a Lightweight External Memory Suffix Array Construction Algorithm.
  In Proc. 2nd International Conference on Algorithms for Big Data (ICABD) 2014.

The latest version of SAscan/pSAscan is available at:
  http://www.cs.helsinki.fi/group/pads/



Compilation and usage
---------------------

1. Download http://libdivsufsort.googlecode.com/files/libdivsufsort-2.0.1.tar.gz
   and install. Make sure to compile libdivsufsort to static 64-bit libraries,
   i.e. set options in the main CMakeLists.txt to

   option(BUILD_SHARED_LIBS "Set to OFF to build static libraries" OFF)
   option(BUILD_DIVSUFSORT64 "Build libdivsufsort64" ON)

2. Compile pSAscan using the provided Makefile

   $ cd src
   $ make

This will produce the executable 'psascan' that allows computing the suffix
array of a given file. For usage, run the 'psascan' program without any
arguments.

Example
~~~~~~~

To compute the suffix array of a file in.txt located in /data01/ using 8GiB
of RAM run the 'psascan' command (assuming you are in the src/ directory) as:

   $ ./psascan /data01/in.txt -m 8192

By default, the resulting suffix array is written to a file matching the
filename of the input text with the .sa5 extension (/data01/in.txt.sa5
in this case). To write the suffix array to a different file, use the
-o flag, e.g.,

   $ ./psascan /data01/in.txt -m 8192 -o /data02/in.txt.suf

The current implementation encodes the output suffix array using unsigned
40-bit integers. For further processing of the suffix array, one should use
the same or compatible encoding. The class implementing the unsigned 40-bit
integers is located in the src/psascan_src/uint40.h file.



Disk space requirements
-----------------------

To compute the suffix array of an n-byte input text, pSAscan needs about
7.5n bytes of disk space. This includes the input (n bytes) and output
(5n bytes).

In the default mode, the 'psascan' program assumes, that there is 6.5n bytes
of free disk space available in the location used as the destination for the
suffix array. This space is used for auxiliary files created during the
computation and to accommodate the output.

The above disk space requirement may in some cases prohibit the use of
algorithm, e.g., if there is enough space (5n) on one physical disk to hold
the suffix array, but not enough (6.5n) to run the algorithm. To still
allow the computation in such cases, the 'psascan' program implements the
-g flag. With this flag, one can force pSAscan to use disk space from two
physically different locations (e.g., on two disks).

More precisely, out of 6.5n bytes of disk space used by pSAscan, about n
bytes is used to store the so-called "gap array". By default, the gap array
is stored along with the suffix array. The -g flag allows explicitly
specifying the location of the gap array. This way, it suffices that there
is only 5.5n bytes of disk space in the location specified as the destination
of the suffix array. The remaining n bytes can be allocated in other location
specified with the -g flag.

Example
~~~~~~~

Assume the location of input/output files and RAM usage as in the example
from the previous section. To additionally specify the location of the gap
array as /data03/in.txt.gap run the 'psascan' command as:

   $ ./psascan /data01/in.txt -m 8192 -o /data02/in.txt.suf -g /data03/in.txt.gap



RAM requirements
----------------

The algorithm does not have a fixed memory requirements. In principle, it
can run with any amount of RAM (though there is some minimal per-thread
amount necessary in the streaming phase). However, since the time complexity
(without logarithmic factors) of the algorithm is O(n^2 / M), where M is the
amount of RAM used in the computation, using more RAM decreases the runtime.
Thus, the best performance is achieved when nearly all unused RAM available
in the system (as shown by the Linux 'free' command) is used for the
computation. Leaving about 5% (but not more than 2GiB) of RAM free is
advised to prevent thrashing.

Example
~~~~~~~

On a machine with 12 physical cores and Hyper-Threading (and thus capable
of simultaneously running 24 threads) it takes about a week to compute a
suffix array of a 200GiB file using 3.5GiB of RAM. Using 120GiB of RAM
reduces the time to less than 12 hours.



Troubleshooting
---------------

1. I am getting "Error: the limit on the maximum number of open files
   is too small (...)".

Solution: The error is caused by the operating system imposing a limit
on the maximum number of files opened by a program. The limit (in Linux
referred to as the soft limit) can be increased with the "ulimit -n newlimit"
command. However, in Linux the soft limit cannot be increased beyond the
so-called "hard limit", which is usually only few times larger than the
soft limit. Furthermore, this is a temporary solution that needs to repeated
every time a new session is started. To increase the limits permanently,
edit (as a root) the file /etc/security/limits.conf and add the following
lines at the end (including the asterisks):

* soft nofile 128000
* hard nofile 128000

This increases the limit to 128000 (use larger values if necessary). The
new limits apply (check with ulimit -n) after starting a new session.

2. Program stops without any error message.

Solution: Most likely the problem occurred during internal-memory sorting.
Re-running the program with -v flag should show the error message.



Limitations / known issues
--------------------------

1. The maximum size of input text is 1TiB (2^40 bytes).
2. The current implementation supports only inputs over byte alphabet.
3. Only texts not containing bytes with value 255 are handled correctly.
   The 255-bytes can be removed from the input text using the tool located
   in the directory tools/delete-bytes-255/ of this package.
4. The current internal-memory suffix sorting algorithm used internally
   in pSAscan works only if the input text is split into segments of
   size at most 2GiB each. Therefore, pSAscan will fail, if the memory
   budget X for the computation (specified with the -m flag) satisfies
   X / p > 10 * 2^31, where p is the number of threads used during
   the computation. On most systems, this is not a severe limitation,
   e.g., for a regular 4-core machine supporting Hyper-Threading (and
   thus capable of simultaneously running 8 threads), pSAscan can utilize
   up to 160GiB of RAM.

The above limitations (except possibly 2) are not inherent to the algorithm
but rather the current implementation. Future releases will most likely
overcome these limitations.



Third-party code
----------------

The pSAscan implementation makes use of some third-party code, in particular:
  - the uint40 class was copied (and slightly modified) from the eSAIS-0.5.2
    algorithm (https://panthema.net/2012/1119-eSAIS-Inducing-Suffix-and-
    LCP-Arrays-in-External-Memory/)
  - pSAscan uses the libdivsufsort-2.0.1 algorithm as the internal
    suffix-sorting routine (https://code.google.com/p/libdivsufsort/)



Terms of use
------------

pSAscan is released under the MIT/X11 license. See the file LICENCE for
more details.

If you use this code, please cite the paper mentioned above and publish
the URL from which you downloaded the code.



Helsinki, June 2015.
Written by Dominik Kempa <dominik.kempa (at) gmail.com>
