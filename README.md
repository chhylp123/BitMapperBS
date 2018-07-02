BitMapperBS: a fast and accurate read aligner for whole-genome bisulfite sequencing
============






Introduction
-------  

Here are the implementations of "BitMapperBS: a fast and accurate read aligner for whole-genome bisulfite sequencing". 
BitMapperBS is an ultra-fast and memory-efficient aligner that is designed for WGBS reads
from directional protocol. 




### Installation ###
(1) Download the source code from Github

    git clone https://github.com/chhylp123/BitMapperBS.git

(2) Build and Install
    
    cd BitMapperBS
    make


### Indexing Genome ###
    
    ./bitmapperBS --index <index file>

### Bisulfite Mapping ###

single-end reads

    walt -i <index file> -r <read files> -o <output file> [options]

paired-end reads

    walt -i <index file> -1 <mate_1 read files> -2 <mate_2 read files> -o <output file> [options]


### Mapping Options ###


| Option | Long Tag | Type | Default | Brief Description |
| :-------------: |:-------------:|:-----:|:-----:| :-----|
| -i      | -index | String | NULL | index file created by ***makedb*** command ( .dbindex) |
| -r      | -reads | String | NULL | list of single-end read files (.fastq or .fq) |
| -1      | -reads1 | String | NULL | list of paired-end read _1 files (.fastq or .fq) |
| -2      | -reads2 | String | NULL | list of paired-end read _2 files (.fastq or .fq) |
| -o      | -output | String | NULL | output file name (.sam or .mr) |
| -m      | -mismatch | Integer | 6 | maximum allowed mismatches |
| -N      | -number | Integer | 1000000 | number of reads to map in one loop |
| -a      | -ambiguous | Boolean | false | randomly output one mapped position for ambiguous reads |
| -u      | -unmapped | Boolean | false | output unmapped reads |
| -C      | -clip | String | empty | clip the specified adaptor |
| -A      | -ag-wild | Boolean | false | map using A/G bisulfite wildcards |
| -P      | -pbat | Boolean | false | map post-bisulfite adaptor tagging reads |
| -b      | -bucket | Integer | 5000 | maximum candidates for a seed |
| -k      | -topk | Integer | 50 | maximum allowed mappings for a read in paired-end mapping |
| -L      | -fraglen | Integer | 1000 | max fragment length in paired-end mapping |
| -t      | -thread | Integer | 1 | number of threads for mapping |

To see the list of options, use "-?" or "-help".

### Examples ###

(1) **Indexing Genome**

For example, to make an index for UCSC hg19

	makedb -c hg19/ -o hg19.dbindex
   
or to make an index for chromosome 2

	makedb -c chr2.fa -o chr2.dbindex

The suffix of the index file should be '.dbindex'.
    
(2) **Bisulfite Mapping**

For example, to mapping reads to human genome hg19

	walt -i hg19.dbindex -r read_1.fq -o reads_1_mapping.sam
    
If mapping the reads from the *_2 reads file, the -A option should be set. This means that all Gs in the reads and genome are converted to As. If -A option is not set, all Cs in the reads and genome are converted to Ts.

    walt -i hg19.dbindex -r read_2.fq -A -o reads_2_mapping.sam

    
If mapping post-bisulfite adaptor tagging reads (PBAT), the -P option should be set. 

    walt -i hg19.dbindex -r read_2.fq -P -o reads_2_mapping.sam
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -P -o paired_mapping_PBAT.sam
    
Additionally, WALT supports comma-separated list of read files. WALT produces one mapping output file for each read file. For single-end mapping, the output file names will be appended "_s1", "_s2", and so on. Notice: except the first file path, all other file paths cannot use '~'. For example, "-r ~/read_file1.fq,~/read_file2.fq" is not allowed. It should be "-r ~/read_file1.fq,/home/read_file2.fq", since linux system doesn't know it is a path except the first one.
	 
	 walt -i hg19.dbindex -r read_file1.fq,read_file2.fq,read_file3.fq -o reads_1_mapping.sam
	 
For paired-end reads, -1 and -2 options are used for the mate read files.
    
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -o paired_mapping.sam
    
Similarly, WALT supports comma-separated list of read files for paired-end mapping. WALT produces one mapping output file for each pair of read files. The output file names will be appended "_p1", "_p2", and so on. One other thing to note is the mate 1 and mate 2 paired files should be in the same order.

	 walt -i hg19.dbindex -1 read_file1_1.fq,read_file2_1.fq,read_file3_1.fq \ 
	                      -2 read_file1_2.fq,read_file2_2.fq,read_file3_2.fq \
	                      -o paired_mapping.sam

The default number of maximum allowed mismatches is 6. The maximum allowed mismatches can be set using -m option.

    walt -i hg19.dbindex -r read_1.fq -m 4 -o reads_1_mapping.sam
    
The option -N sets the number of reads to mapping in each loop. If N is larger, the program takes large memory, especially for paired-end mapping. For Human Genome, if N is 1000000, both single-end and paired-end mapping take about 16 Gb memory. If N is 5000000, single-end mapping takes about 18 Gb memory, while paired-end mapping takes about 28 Gb memory. If N is set to be larger than 10000000, the program will set N to be 10000000 since when N is too large the program will take large memory but it will not be significantly faster. The estimate memory for single-end mapping is 15 + N * (2 * rl + rnl + 16) / (1024^3) Gb, and for paired-end mapping is 15 + N * (4 * rl + 2 * rnl + 24 * k + 16) / (1024^3) Gb. The index size of human genome is about 15 Gb. N is the number of reads to map in one loop. rl is the length of reads (WALT supports mix of different lengths reads, so here rl is estimation of the average length). rnl is the length of read names. k is maximum allowed mappings for a read in paired-end mapping.
    
    walt -i hg19.dbindex -r read_1.fq -N 1000000 -o reads_1_mapping.sam
    
To output the ambiguous mapped reads or unmapped reads, -u and -a options should be set. If the output format is MR, the ambigous mapped reads and unmapped reads will output to a separated file, respectively. For paired-end mapping, there will be a unmapped file and an ambiguous file for each of the mate 1 and mate 2 reads file. If the output format is SAM, uniquely mapped, ambiguous mapped, unmapped reads are output to the same file with different FLAGs.
    
    walt -i hg19.dbindex -r read_1.fq -u -a -o reads_1_mapping.sam
    
To trim 3' end adaptor sequence, -C option should be set. For paired-end read mapping, using T_adaptor[:A_adaptor] format to set the adaptor. If only one adaptor seqeunce is given, the same adaptor sequence will be used for both T-rech and A-rich reads.

    walt -i hg19.dbindex -r read_1.fq -C AGATCGGAAGAGC -o reads_1_mapping.sam
    walt -i hg19.dbindex -r read_2.fq -A -C AGATCGGAAGAGC -o reads_2_mapping.sam
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -C AGATCGGAAGAGC -o paired_mapping.sam
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -C AGATCGGAAGAGC:AGATCGG -o paired_mapping.sam
    
By default, WALT produces SAM format files. If you plan to get MR file, please make sure the suffix of the output file is ".mr".

	 walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -o paired_mapping.mr
	 
WALT produces a separate mapping statistics file (.mapstats) for each run. It includes percentages of uniquely mapped, ambiguously mapped and unmapped reads.
	 
	 
### Output Format ###

WALT supports Tab-delimited SAM and MR output formats. By default, WALT produces SAM format output files. To get MR format, the suffix of the output file should be ".mr".

**SAM Format**

* QNAME (read name)
* FLAG (bitwise FLAG)
* RNAME (chromosome name)
* POS (start position, 1-based, 0 means unmapped)
* MAPQ (255 in WALT)
* CIGAR (CIGAR string)
* RNEXT (chromosome name of the mate read)
* PNEXT (start position of the mate read)
* TLEN (observed segment length)
* SEQ (read sequence)
* QUAL (quality sequence read from fastq file)
* NM-tag (number of mismatches)

By default, WALT only outputs uniquely mapped reads. If -u or -a option is set, the ambiguous mapped or unmapped reads will be in the same SAM file with different FLAGs. The FLAG 0x100 will be set for ambiguous mapped reads, and the FLAG 0x4 will be set for unmapped reads. 

**MR Format**

* RNAME (chromosome name)
* SPOS (start position, 0-based)
* EPOS (end position, 0-based)
* QNAME (read name)
* MISMATCH (number of mismatches)
* STRAND (forward or reverse strand)
* SEQ
* QUAL

If paired-end reads are mapped in proper pair, the QNAME is added "FRAG:" in the beginning of the read name, the STRAND is the strand of the first mate mapped and SEQ and QUAL is merged according to their mapping positions. The overlap segment of SEQ and QUAL is from the mate 1 or mate 2 and it is the one with less number of 'N' in the read sequence. MISMATCH is the sum of mismatches in the mate 1 and mismatches in the mate 2. If paired-end reads are not mapped in proper pair, they are treated as single-end reads. If -u and/or -a option is set, the ambiguous mapped or unmapped reads are output in separate files.

    
### Contacts ###

***Haifeng Chen***   *haifengc@usc.edu*

***Andrew Smith***   *andrewds@usc.edu*

***Ting Chen***   *tingchen@usc.edu*


### Citation ###
WALT: fast and accurate read mapping for bisulfite sequencing. Bioinformatics, 32(22), pp.3507-3509. [Supplementary](https://github.com/smithlabcode/walt/blob/master/doc/Supplementary%20Data.pdf)

### Copyright ###

    Copyright (C) 2015 University of Southern California
                   Andrew D. Smith and Ting Chen
                   
    Authors: Haifeng Chen, Andrew D. Smith and Ting Chen

    WALT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    WALT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with WALT.  If not, see <http://www.gnu.org/licenses/>.

Note
-------
* We adopt the SACA-K algorithm [1] to build the suffix array, and build BWT from suffix array. As such when building the index, the memory requirement of FMtree, Original_s and Original_v is about 5 times larger than that of the input text.

* Please note that FMtree, Original_s and Original_v do not share a same index. For each method, users should build its own index.

References
-------


[1] Nong G. Practical linear-time O (1)-workspace suffix sorting for constant alphabets[J]. ACM Transactions on Information Systems (TOIS), 2013, 31(3): 15.
