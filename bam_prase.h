/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/


#ifndef __BAM_PRASE__
#define __BAM_PRASE__

#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <htslib/klist.h>
#include <htslib/thread_pool.h>
#include <htslib/bgzf.h>


typedef struct
{
	samFile* bam_file;
	bam_hdr_t* bam_header;
} bam_operations;



void init_bam_header(char* tmp_sam_file_name, char* bam_file_name);
void close_bam_file();
void write_alignment(const char* alignment);
void write_alignment_muti_thread(const char* alignment);

#endif
