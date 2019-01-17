/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#ifndef __READ__
#define __READ__

#include <math.h>
#include<stdint.h>
#include <zlib.h>

/**
typedef struct
{
  char name[SEQ_MAX_LENGTH/2];
  char seq[SEQ_MAX_LENGTH/2];
  char rseq[SEQ_MAX_LENGTH/2];
  char qual[SEQ_MAX_LENGTH/2];
  int hits[1];
  int length;
} Read;
**/
#define READ_NAME_INIT_LENGTH 20
#define READ_SEQ_INIT_LENGTH 155
#define E_FORMAT 0
#define FASTQ 1
#define FASTQGZ 2



typedef struct
{
	char* name;
	char* seq;
	char* rseq;
	char* qual;
	uint16_t length;
	uint16_t name_size;
	uint16_t seq_size;
	///char hits[1];  ///这个到时候要改了
} Read;




typedef struct
{
	Read* read1;
	Read* read2;
	long long sub_block_read_number;

} Read_buffer_pe_sub_block;


typedef struct Read_buffer_pe
{
	Read_buffer_pe_sub_block* sub_block;
	///这个值的单位是单个read
	long long sub_block_inner_size;
	///这两个单位是subblock
	long long sub_block_size;
	long long sub_block_number;
	///int all_read_end = 0;
	int all_read_end;

	Read_buffer_pe() :all_read_end(0)
	{};

} Read_buffer_pe;







typedef struct
{
	Read* read;
	long long sub_block_read_number;

} Read_buffer_single_sub_block;


typedef struct Read_buffer_single
{
	Read_buffer_single_sub_block* sub_block;
	///这个值的单位是单个read
	long long sub_block_inner_size;
	///这两个单位是subblock
	long long sub_block_size;
	long long sub_block_number;
	///int all_read_end = 0;
	int all_read_end;


	Read_buffer_single() :all_read_end(0)
	{};

} Read_buffer_single;














void init_Single_Seq_input_buffer(int thread_number);
void* input_single_reads_muti_threads(void*);
void* input_single_reads_muti_threads_pbat(void*);

int get_single_reads_mul_thread_pbat(Read_buffer_single_sub_block* curr_sub_block);


void init_single_sub_block_single(Read_buffer_single_sub_block* tmp_sub_block);

void init_Pair_Seq_input_buffer(int thread_number);
void* input_pe_reads_muti_threads(void*);
void init_single_sub_block_pe(Read_buffer_pe_sub_block* tmp_sub_block);
int get_pe_reads_mul_thread(Read_buffer_pe_sub_block* curr_sub_block);
int get_single_reads_mul_thread(Read_buffer_single_sub_block* curr_sub_block);

void init_single_read(Read* read);

inline int inputReads_single_pure(Read *seqList1);

int initiReadAllReads(char *fileName1, char *fileName2, unsigned char pairedEnd, int* format);
void Unmapped_process(char *fileName);
int inputReads_single( Read *seqList1);
int inputReads_paired( Read *seqList1, Read *seqList2);
int inputReads_single_directly(Read *seqList1);
int inputReads_single_directly_pbat(Read *seqList1);
int inputReads_paired_directly(Read *seqList1, Read *seqList2);
int exchange_two_reads();

#endif
