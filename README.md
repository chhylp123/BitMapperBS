BitMapperBS: a fast and accurate read aligner for whole-genome bisulfite sequencing
============

### Introduction ###
-------  

Here are the implementations of "BitMapperBS: a fast and accurate read aligner for whole-genome bisulfite sequencing". 
BitMapperBS is an ultra-fast and memory-efficient aligner that is designed for WGBS reads from directional protocol. 

- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `(update on August 10, 2019) Please do not use version 1.0.2.0, which has some problems. The current version is 1.0.2.1.`

### Build Requirements ###

(1) G++.

(2) CMake.

(3) CMake-supported build tool.

(4) zlib, libbz2 and liblzma libraries. In Ubuntu, please try: sudo apt-get install liblzma-dev zlib1g-dev libbz2-dev.

### Hardware&software requirements ###

(1) CPU must support AVX2 or SSE4.2 instructions.

(2) When building the index for the human genome, BitMapperBS requires about 10GB RAM and 60GB disk space.

(3) When aligning the bs-seq to the human genome, BitMapperBS requires about 7GB RAM. 

### Supported platforms ###

BitMapperBS has been successfully tested using six CPU threads on a computer with a six-core Intel Core i7-8770k processor and 64GB RAM, running Ubuntu 16.04. The indexes, reference genomes and reads were stored in a Solid State Drive (SSD) to minimize the loading time.
It is also actively used by Computational Biology of Aging Group and BGI Genomics to analyze WGBS data.

### Docker container ###

BitMapper can be run as binary or as a docker container. For example, if we assume that reference genomes and samples are in /data/indexes and /data/samples folders, then:
* to build index of the reference genome:
```
docker run -v /data:/data quay.io/comp-bio-aging/bit_mapper_bs:latest /opt/BitMapperBS/bitmapperBS --index /data/indexes/HUMAN/29/GRCh38.primary_assembly.genome.fa  --index_folder human_bs_index
```
* to align sequence to bitmapper index:
```
docker run -v /data:/data quay.io/comp-bio-aging/bit_mapper_bs:latest /opt/BitMapperBS/bitmapperBS --search /data/indexes/human_bs_index  --seq1 /data/samples/SRR948855/SRR948855_1.fastq.gz --seq2 /data/samples/SRR948855/SRR948855_2.fastq.gz --sensitive --pe -t 8 --mapstats --bam -o SRR948855.bam
```

### Installation ###
(1) Download the source code from Github

    git clone https://github.com/chhylp123/BitMapperBS.git

(2) Build (CPU must support AVX2 instructions or SSE4.2 instructions)
    
    cd BitMapperBS
    make
    
- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `(update on October 10, 2018) If system reports: "cmake: not found" or "psascan_src/inmem_psascan_src/divsufsort_template.h:42:24: fatal error: divsufsort.h: not found", please install CMake in your system.`

- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `(update on November 28, 2018) If system reports: "fatal error: zlib.h: no such file or directory" or "fatal error: bzlib.h: No such file or directory" or "fatal error: lzma.h: No such file or directory", please install zlib, libbz2 and liblzma libraries. In Ubuntu, please try: sudo apt-get install liblzma-dev zlib1g-dev libbz2-dev. `



(3) (update on October 10, 2018) In most cases, BitMapperBS can be compiled from source code automatically, and is able to be implemented successfully. However, in some rare cases (e.g, old version of Linux operating system), BitMapperBS may report error message when building index. For example, report: "sh: 1: ./psascan: not found". This is because BitMapperBS utlizes psascan to build FM-index, and psascan (binary file) cannot be compiled from source code automatically. In this case, please compile psascan manually (https://www.cs.helsinki.fi/group/pads/pSAscan.html), and copy it (binary file of psascan) to the folder of BitMapperBS. 

### - ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `Please Note!!!` ###

>1. (update on October 10, 2018) In most cases, BitMapperBS can be compiled from source code automatically, and is able to be implemented successfully. However, in some rare cases (e.g, old version of Linux operating system), BitMapperBS may report error message when building index. For example, report: "sh: 1: ./psascan: not found". This is because BitMapperBS utlizes psascan to build FM-index, and psascan (binary file) cannot be compiled from source code automatically. In this case, please compile psascan manually, and copy it (binary file of psascan) to the folder of BitMapperBS. The detailed steps are listed as follows: 

 >> (1) Download psascan from https://www.cs.helsinki.fi/group/pads/pSAscan.html, and complie it from source code.

 >> (2) Copy psascan (binary file) to the folder of BitMapperBS.


>2. (update on November 28, 2018)  When compiling BitMapperBS, if you get the error message "fatal error: zlib.h: no such file or directory" or "fatal error: bzlib.h: No such file or directory" or "fatal error: lzma.h: No such file or directory", please install zlib, libbz2 and liblzma libraries. In Ubuntu, please try:
   >> sudo apt-get install liblzma-dev zlib1g-dev libbz2-dev


>3. (update on November 28, 2018) Although BitMapperBS itself is significantly faster than other methods, the slow disk I/O cannot be accelerated. In practice, the most serious bottleneck of BitMapperBS is the poor performance of disk I/O, especially when using multiple CPU threads. Thus, if you want to run BitMapperBS using many CPU threads, we suggest you to adopt at least one of the following strategies: 
   >> (1) To reduce the amount of disk I/O, you can use the compressed fastq files (.fastq.gz or .fq.gz format) rather than the uncompressed raw files (.fastq or .fq format).

   >> (2) To reduce the amount of disk I/O, you can output the mapping results in BAM format (using the option --bam) rather than SAM format.

   
   >> (3) The input files and output files of BitMapperBS (e.g., the read files and the output SAM or BAM files) can be saved in fast solid state drives (SSD) storage devices, rather than slow hard disk drive (HDD) storage devices.


If you have problem with the "make" part described above, please contact us (chhy@mail.ustc.edu.cn).

Usage
=====

### Indexing Genome ###
    
    ./bitmapperBS --index <genome file name>

The suffix of the index file should be '.index.*'.

- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `(update on October 10, 2018) If BitMapperBS reports: "sh: 1: ./psascan: not found" when building index, this means that psascan did not compiled and installed successfully. BitMapperBS utlizes psascan to build FM-index, and psascan (binary file) cannot be compiled from source code automatically. In this case, please compile psascan manually (https://www.cs.helsinki.fi/group/pads/pSAscan.html), and copy it (binary file of psascan) to the folder of BitMapperBS.`


### Quality and Adapter Trimming ###

We recommed users to use Trim_Galore to perform the quality and adapter trimming. Please see https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. 

### Bisulfite Mapping ###

single-end reads

    ./bitmapperBS --search <genome file name> --seq <read file name> [options]

paired-end reads

    ./bitmapperBS --search <genome file name> --seq1 <read1 file name> --seq2 <read2 file name> [options]


output mapping results in BAM format

    ./bitmapperBS --search <genome file name> --seq <read file name> --bam -o output.bam [options]


### Methylation Extracting ###


We recommend users to first remove the duplicates by Picard or samtools, and then use MethylDackel to extract methylation information. Please see https://github.com/dpryan79/methyldackel. Please note that MethylDackel maybe quite slow when using too many CPU threads. So we personally recommend users to use 2, 3 or 4 CPU threads for MethylDackel. Too many threads cannot accelerate methylation extracting step.


### General Options ###


| Option | Short Option | Type | Default | Brief Description |
| :------: | :---------------: | :-----:|:-----:| :-----|
| --help | -h | String | NULL | Show the help information. |
| --version | -v | String | NULL | Show current version of BitMapperBS. |

### Index Options ###

| Option | Short Option | Type | Default | Brief Description |
| :------: | :---------------: | :-----:|:-----:| :-----|
| --index | -i | String | NULL | Generate an index from the specified fasta file. |
| --index_folder | NULL | String | NULL | Set the folder that stores the genome indexes. If this option is not set, the indexes would be stores in the same folder of genome (input fasta file). |


### Mapping Options ###



| Option | Short Option | Type | Default | Brief Description |
| :------: | :---------------: | :-----:|:-----:| :-----|
| --search | NULL| String | NULL | Search in the specified genome. If the indexes of this genome are built without "--index_folder", please provide the path to the fasta file when aligning. Otherwise please provide the path to the index folder (set by "--index_folder" during indexing).|
| --fast | NULL| NULL | NULL | Set bitmapperBS in fast mode (default). Only available for paired-end mode.|
| --sensitive | NULL| NULL | NULL | Set bitmapperBS in sensitive mode. Only available for paired-end mode.|
| --seq | NULL| String | NULL | Provide the name of single-end read file (.fastq/.fq/.fastq.gz/.fq.gz format). |
| --seq1 | NULL| String | NULL | Provide the name of paired-end read_1 file (.fastq/.fq/.fastq.gz/.fq.gz format). |
| --seq2 | NULL| String | NULL | Provide the name of paired-end read_2 file (.fastq/.fq/.fastq.gz/.fq.gz format). |
| -o | -o | String | stdout (Standard output) | Provide the name of output file (SAM or BAM format). |
| --sam | NULL| NULL | NULL | Output mapping results in SAM format (default). |
| --bam | NULL| NULL | NULL | Output mapping results in BAM format. |
| -e | -e | Double | 0.08 | Set the edit distance rate of read length, which is between 0 and 1. |
| --min | NULL | Int | 0 | Minimum observed template length between a pair of end sequences. |
| --max | NULL | Int | 500 | Maximum observed template length between a pair of end sequences. |
| --threads | -t | Int | 1 | Set the number of CPU threads. |
| --pbat | NULL | NULL | NULL | Map the bs-seq from pbat protocol. |
| --unmapped_out | NULL | NULL | NULL | Report unmapped reads. |
| --ambiguous_out | NULL | NULL | NULL | Random report one of hit of each ambiguous mapped read. |
| --mapstats | NULL | String | NULL | Output the statistical information of read alignment into file.|
| --phred33 | NULL | NULL | NULL | Input read qualities are encoded by Phred33 (default). |
| --phred64 | NULL | NULL | NULL | Input read qualities are encoded by Phred64. |
| --mp_max | NULL | Int | 6 | Maximum mismatch penalty. |
| --mp_min | NULL | Int | 2 | Minimum mismatch penalty. |
| --np | NULL | Int | 1 | Ambiguous character (e.g., N) penalty. |
| --gap_open | NULL | Int | 5 | Gap open penalty. |
| --gap_extension | NULL | Int | 3 | Gap extension penalty. |



### Examples ###

(1) **Indexing Genome**

For example, to make an index for human genome (GRCH38):

	./bitmapperBS --index human_genome.fa
   
The suffix of the index file should be '.index.*'.

- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `(update on October 10, 2018) If BitMapperBS reports: "sh: 1: ./psascan: not found" when building index, this means that psascan did not compiled and installed successfully. BitMapperBS utlizes psascan to build FM-index, and psascan (binary file) cannot be compiled from source code automatically. In this case, please compile psascan manually (https://www.cs.helsinki.fi/group/pads/pSAscan.html), and copy it (binary file of psascan) to the folder of BitMapperBS.`

(2) **Quality and Adapter Trimming**

For example, to trim paired-end reads:

	trim_galore --paired read_1.fq read_2.fq
    
(3) **Bisulfite Mapping**

For example, to map reads to human genome (GRCH38) in single-end mode using 6 CPU threads:

	./bitmapperBS --search ../../ssd/human_genome.fa --seq ../../ssd/read.fq -t 6
    
If map the reads from the *_2 reads file or the pbat protocol, the --pbat option should be set:

    	./bitmapperBS --search ../../ssd/human_genome.fa --seq ../../ssd/read_2.fq --pbat -t 6

to map reads to human genome (GRCH38) in paired-end mode using 6 CPU threads:

	./bitmapperBS --search ../../ssd/human_genome.fa --seq1 ../../ssd/read_1.fq --seq2 ../../ssd/read_1.fq -t 6



to output mapping results to the file named "output.bam" in BAM format

    ./bitmapperBS --search ../../ssd/human_genome.fa --seq ../../ssd/read.fq -t 6 --bam -o output.bam


The default maximum allowed edit distance is 8% of read length (0.08). This option can be changed using -e option. In this example, the maximum allowed edit distance is set to 4% of read length:

    ./bitmapperBS --search ../../ssd/human_genome.fa --seq ../../ssd/read.fq -t 6 -e 0.04


(4) **Methylation Extracting**

The output file of BitMapperBS must be first sorted into a coordinate-sorted BAM file by utilizing samtools. After that, duplicate alignments should be removed by Picard or samtools. At last, the methylation information is extracted using MethylDackel (see https://github.com/dpryan79/methyldackel). Please note that MethylDackel maybe quite slow when using too many CPU threads. So we personally recommend users to use 2, 3 or 4 CPU threads for MethylDackel. Too many threads cannot accelerate methylation extracting step.


### Changelog ###

(16) August 24, 2019: version 1.0.2.1 released. 

    >> Fix the bug in version 1.0.2.0.
    >> Add support of BitMapperBS to output MAPQ like Bowtie2.
    >> Recalculate alignment with large or complex gaps.
    >> For the '--mapstats' option, the output file name can be defined by users.
    >> Add '--phred33', '--phred64', '--mp_max', '--mp_min', '--np', '--gap_open', and '--gap_extension' options.

(15) June 29, 2019: version 1.0.1.6 released. 

    >> Remove the `--pe` option. The paired-end mode or the single-end mode can be determined automatically.
    >> All logging information is printed to standard error.
    >> If `-o` is not specified, the alignment results would be printed to standard output.
    >> Fix a minor bug of the paired-end read name.

(14) June 28, 2019: version 1.0.1.5 released. 

    >> Fix a minor bug when mapping.

(13) June 18, 2019: version 1.0.1.4 released. 

    >> Fix the bug for methylation extraction when loading index.

(12) April 13, 2019: version 1.0.1.3 released. 

    >> Update the option for methylation extraction.

(12) March 10, 2019: version 1.0.1.2 released. 

    >> Fix the bug when indexing in extreme cases.

(11) March 9, 2019: version 1.0.1.1 released. 

    >> Fix the bug when indexing in extreme cases.

(10) March 9, 2019: version 1.0.1.0 released. 

    >> Fix the bug about the optopn "--index_folder" when indexing.

(9) March 8, 2019: version 1.0.0.9 released. 

    >> Fix the bug that may be happened with old version of G++.

(8) January 17, 2019: version 1.0.0.8 released. 

    >>Add support of BitMapperBS to report the statistical information of read alignment into file.
    >>When the input read files are in compressed fastq files (.fq.gz or .fastq.gz), BitMapperBS is sligntly faster at the expense of sligntly higher CPU usage. 

(7) January 8, 2019: version 1.0.0.7 released. 

    >>Add support of BitMapperBS to report unmapped reads and ambiguous mapped read.


(6) December 27, 2018: version 1.0.0.6 released. 

    >>BitMapperBS fixed the bug of the option --min when aligning paired-end reads. 
    >>Add support of BitMapperBS to report mismatch and indel rate in alignment.

(5) December 22, 2018: version 1.0.0.5 released. 

    >>BitMapperBS recalculated CIGAR in SAM file in some extreme cases. 

(4) December 11, 2018: version 1.0.0.4 released. 

    >>BitMapperBS fixed the bug of the QNAME in SAM file when aligning paired-end reads. 
    >>Remove the reads which are mapped off the end of the reference.

(3) November 28, 2018: version 1.0.0.3 released. 

    >>Add support of BitMapperBS to output mapping results in BAM format.


(2) November 16, 2018: version 1.0.0.2 released. 

    >>BitMapperBS fixed the bug of the TLEN field of SAM format in output SAM file when aligning the paired-end reads.


(1) November 9, 2018: version 1.0.0.1 released. 
     
     >>Add support of BitMapperBS to accept read files compressed by gzip (.fq.gz or .fastq.gz).
     >>BitMapperBS fixed the bug of the SAM flag when aligning the single-end reads from pbat protocol.
     >>BitMapperBS can check if the AVX2 instructions are supported by CPU automatically. If AVX2 is supported, the AVX2 version of BitMapperBS will be compiled, otherwise the SSE4.2 version of BitMapperBS will be compiled. Please note that the AVX2 version of BitMapperBS may be slightly slower than SSE4.2 version of BitMapperBS.




### Contacts ###

Haoyu Cheng: hcheng@jimmy.harvard.edu


### Note ###
-------
* We adopt the pSAscan algorithm [1] to build the suffix array, and build BWT from suffix array.



### References ###
-------


[1] Kärkkäinen J, Kempa D, Puglisi S J. Parallel external memory suffix sorting[C]//Annual Symposium on Combinatorial Pattern Matching. Springer, Cham, 2015: 329-342.
