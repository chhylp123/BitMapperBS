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

### Mapping Options ###



| Option | Type | Default | Brief Description |
| :-------------: |:-----:|:-----:| :-----|
| -h | String | NULL | Show the help information. |
| -v/--version | String | NULL | Show current version of BitMapperBS. |
| -i/--index | String | NULL | Generate an index from the specified fasta file. |





#### General Options ####

 -v|--version		Current Version.

 -h			Show the help file.



#### Indexing Options ####

 --index [file]		Generate an index from the specified fasta file. 


#### Searching Options ####

 --search [file]	Search in the specified genome. Provide the path to the fasta file. Index file should be in the same directory.


 --pe 			Search will be done in paired-end mode.


 --seq [file]		Input sequences in fastq format [file]. This option is used for single-end reads.


 --seq1 [file]		Input sequences in fastq format [file] (First file). Use this option to indicate the first file of paired-end reads. 


 --seq2 [file]		Input sequences in fastq format [file] (Second file). Use this option to indicate the second file of paired-end reads.  

 -o [file]		Output of the mapped sequences. The default is "output".


 -e [float]		Maximum allowed edit distance (default 8% of the read length).


 --min [int]		Min distance allowed between a pair of end sequences (default: 0).


 --max [int]		Max distance allowed between a pair of end sequences (default: 500).



 --threads, -t [int]	Set the number of CPU threads (default: 1).


 --pbat 		Mapping the BS-seq from pbat protocol.



To see the list of options, use "-?" or "-help".

Note
-------
* We adopt the pSAscan algorithm [1] to build the suffix array, and build BWT from suffix array.


References
-------


[1] Kärkkäinen J, Kempa D, Puglisi S J. Parallel external memory suffix sorting[C]//Annual Symposium on Combinatorial Pattern Matching. Springer, Cham, 2015: 329-342.
