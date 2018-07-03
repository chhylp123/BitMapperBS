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

For example, to make an index for human genome (GRCH38):

	./bitmapperBS --index human_genome.fa
   
The suffix of the index file should be '.index.*'.
    
(2) **Bisulfite Mapping**

For example, to map reads to human genome (GRCH38) in single-end mode using 6 CPU threads:

	./bitmapperBS --search ../../ssd/human_genome.fa --seq ../../ssd/read.fq -t 6
    
If map the reads from the *_2 reads file or the pbat protocol, the --pbat option should be set:

    	./bitmapperBS --search ../../ssd/human_genome.fa --seq ../../ssd/read_2.fq --pbat -t 6

to map reads to human genome (GRCH38) in paired-end mode using 6 CPU threads:

	./bitmapperBS --search ../../ssd/human_genome.fa --seq1 ../../ssd/read_1.fq --seq2 ../../ssd/read_1.fq --pe -t 6


The default maximum allowed edit distance is 8% of read length (0.08). This option can be changed using -e option. In this example, the maximum allowed edit distance is set to 4% of read length:

    ./bitmapperBS --search ../../ssd/human_genome.fa --seq ../../ssd/read.fq -t 6 -e 0.04






Note
-------
* We adopt the pSAscan algorithm [1] to build the suffix array, and build BWT from suffix array.


References
-------


[1] Kärkkäinen J, Kempa D, Puglisi S J. Parallel external memory suffix sorting[C]//Annual Symposium on Combinatorial Pattern Matching. Springer, Cham, 2015: 329-342.
