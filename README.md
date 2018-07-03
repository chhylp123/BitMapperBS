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
    
    ./bitmapperBS --index <genome file name>

### Bisulfite Mapping ###

single-end reads

    ./bitmapperBS --search <genome file name> --seq <read file name> [options]

paired-end reads

    ./bitmapperBS --search <genome file name> --seq1 <read1 file name> --seq2 <read2 file name> --pe [options]


### General Options ###


| Option | Short Option | Type | Default | Brief Description |
| :------: | :---------------: | :-----:|:-----:| :-----|
| --help | -h | String | NULL | Show the help information. |
| --version | -v | String | NULL | Show current version of BitMapperBS. |

### Index Options ###

| Option | Short Option | Type | Default | Brief Description |
| :------: | :---------------: | :-----:|:-----:| :-----|
| --index | -i | String | NULL | Generate an index from the specified fasta file. |


### Mapping Options ###



| Option | Short Option | Type | Default | Brief Description |
| :------: | :---------------: | :-----:|:-----:| :-----|
| --search | NULL| String | NULL | Provide the path to the fasta file when aligning.|
| --pe | NULL| NULL | NULL | Searching will be done in paired-end mode. |
| --seq | NULL| String | NULL | Provide the name of single-end read file. |
| --seq1 | NULL| String | NULL | Provide the name of paired-end read_1 file. |
| --seq2 | NULL| String | NULL | Provide the name of paired-end read_2 file. |
| -o | -o | String | output | Provide the name of output file. |
| -e | -e | Double | 0.08 | Set the edit distance rate of read length, which is between 0 and 1. |
| --min | NULL | Int | 0 | Minimum distance between a pair of end sequences. |
| --max | NULL | Int | 0 | Maximum distance between a pair of end sequences. |
| --threads | -t | Int | 1 | Set the number of CPU threads. |
| --pbat | NULL | NULL | NULL | Map the bs-seq from pbat protocol. |






### Examples ###

(1) **Indexing Genome**

For example, to make an index for UCSC hg19

	./bitmapperBS --index Homo_sapiens.GRCh38.dna.primary_assembly.fa
   
The suffix of the index file should be '.dbindex.*'.
    
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






Note
-------
* We adopt the pSAscan algorithm [1] to build the suffix array, and build BWT from suffix array.


References
-------


[1] Kärkkäinen J, Kempa D, Puglisi S J. Parallel external memory suffix sorting[C]//Annual Symposium on Combinatorial Pattern Matching. Springer, Cham, 2015: 329-342.
