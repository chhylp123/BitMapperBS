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
#include <htslib/sam.h>
#include "Schema.h"


typedef struct
{
	samFile* bam_file;
	bam_hdr_t* bam_header;
	char* output_name; 
} bam_operations;



typedef struct
{
	kstring_t line;
	bam1_t *b;
} bam_output_cell;


typedef struct
{
	bam_hdr_t* bam_header;
	kstring_t line;
	bam1_t **b;
	long long size;
	long long number;
} bam_output_group;


typedef struct
{
	kstring_t line;
	bam1_t *b;
	BGZF* bgzf;
	bgzf_buffer result;
} bam_phrase;








void init_bam_buffer(bam_phrase* buf);

void init_bam_output_group(bam_output_group* group, long long total_size);

void phrase_alignment_only(char* alignment, long long alignment_length, bam_output_group* bam_groups);
void write_alignment_only(bam_output_group* bam_groups);

void get_bam_header(bam_hdr_t** bam_header);

void init_bam_output_cell(bam_output_cell* cell);
void write_alignment_directly(char* alignment, long long alignment_length, bam_output_cell* cell);


void close_bam_file_rename(char* old_name, char* new_name);

void init_bam_header(char* bam_file_name, _rg_name_l  *_ih_refGenName, int refChromeCont, int argc, char *argv[]);
void init_bam_file_from_sam(char* file_name, char* outputFileName,
_rg_name_l  *_ih_refGenName, int refChromeCont, int argc, char *argv[]);


void close_bam_file();
void write_alignment(const char* alignment);
void write_alignment_muti_thread(const char* alignment);
void write_alignment_group(char* alignment, const long long length);
void write_alignment_group_pre_allocate(char* alignment, const long long length, const long long number);

void convert_string_to_bam(char* alignment, long long alignment_length, bam_phrase* bam_groups);


void flush_bam_buffer(bam_phrase* bam_groups);

#endif