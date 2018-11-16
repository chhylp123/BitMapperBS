/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#ifndef __OUTPUT__
#define __OUTPUT__

#define	FORWARD 0
#define REVERSE 1

#include "Schema.h"


typedef struct
{
  char		*tag;
  char		type;
  char		cVal;
  int		iVal;
  float		fVal;
  char		*sVal;
} OPT_FIELDS;

typedef struct
{
  char			*QNAME;
  short			FLAG;
  char			*RNAME;
  unsigned int			POS;
  unsigned char  	MAPQ;
  char			*CIGAR;
  char          *CIGAR_CHAR;
  int           *CIGAR_SIZE;
  char			*MRNAME;
  int			MPOS;
  int			ISIZE;
  char			*SEQ;
  char			*QUAL;
  int			optSize;
  OPT_FIELDS	        *optFields;
} SAM;

int Output_gene(char *fileName);
///void (*finalizeOutput)();
void finalizeOutput();
void multi_outputQ(SAM map, int thread_id);
FILE* get_Ouput_Dec();
void OutPutSAM_Nounheader(_rg_name_l  *_ih_refGenName, int refChromeCont, int argc, char *argv[]);

///是双端的read
#define IS_PAIRED 1
///双端read符合要求的比对上(满足一定距离)
#define MATE_MAPPED 2
///read本身没有比对上
#define SINGLE_UNMAPPED 4
///双端read没有符合要求的比对上(满足一定距离)
#define MATE_UNMAPPED 8
///read本身比对到反向互补链
#define RC_MAPPED 16        ///OX10
///双端read的另一条比对到反向互补链，也就是该read本身比对到正向链了
#define RC_MATE_MAPPED 32  ///OX20
///双端read中的read1
#define READ1 64
///双端read中的read2
#define READ2 128
///有mutiple mapping
#define MUTIREAD 128


typedef struct
{
	char* buffer;
	long long size;
	long long length;

} Output_buffer_sub_block;


typedef struct
{
	Output_buffer_sub_block* sub_buffer;

	///这两个单位是subblock
	long long sub_block_size;
	long long sub_block_number;
	int all_buffer_end = 0;

} Output_buffer;


typedef struct
{
	Methylation* sub_buffer;

	///这两个单位是subblock
	long long sub_block_size;
	long long sub_block_number;
	int all_buffer_end = 0;

} Output_methy_buffer;

void init_buffer_sub_block(Output_buffer_sub_block* sub_block);
void init_output_buffer(int thread_number);
void* pop_buffer(void*);
void push_results_to_buffer(Output_buffer_sub_block* sub_block);
void finish_output_buffer();


void init_output_methy_buffer(int thread_number);
void* pop_methy_buffer(void*);
void push_methy_to_buffer(Methylation* methy);


#define OFFSET 32
#define INT_LENGTH 20 //2^64 -1 = 18446744073709551615, 占了20位

///扩大空间留的都有余量，所以向后加一点\t之类的不需要再检查了
///如果输入变量是char*,且长度不知
inline void output_to_buffer_char_no_length(Output_buffer_sub_block* sub_block, char* input)
{
	


	int input_length = strlen(input);



	///更新了length
	sub_block->length = sub_block->length + input_length;

	///如果空间足够的,则size不用更新
	if (sub_block->length + OFFSET< sub_block->size)
	{
		memcpy(sub_block->buffer + sub_block->length - input_length, input, input_length);
		///sub_block->buffer[sub_block->length] = '\0';
	}
	else
	{
		///空间不够，size需要扩大
		sub_block->size = (sub_block->length + OFFSET)* 2 + 1;
		char* tmp = (char*)malloc(sub_block->size);
		memcpy(tmp, sub_block->buffer, sub_block->length - input_length);
		free(sub_block->buffer);
		sub_block->buffer = tmp;
		memcpy(sub_block->buffer + sub_block->length - input_length, input, input_length);
		///sub_block->buffer[sub_block->length] = '\0';

	}


}



///扩大空间留的都有余量，所以向后加一点\t之类的不需要再检查了
///如果输入变量是char*,且长度已知
inline void output_to_buffer_char_length(Output_buffer_sub_block* sub_block, char* input, int input_length)
{


	///更新了length
	sub_block->length = sub_block->length + input_length;



	///如果空间足够的,则size不用更新
	if (sub_block->length + OFFSET< sub_block->size)
	{
		memcpy(sub_block->buffer + sub_block->length - input_length, input, input_length);
		///sub_block->buffer[sub_block->length] = '\0';

	}
	else
	{
		///空间不够，size需要扩大
		sub_block->size = (sub_block->length + OFFSET) * 2 + 1;
		char* tmp = (char*)malloc(sub_block->size);
		memcpy(tmp, sub_block->buffer, sub_block->length - input_length);
		free(sub_block->buffer);
		sub_block->buffer = tmp;
		memcpy(sub_block->buffer + sub_block->length - input_length, input, input_length);
		///sub_block->buffer[sub_block->length] = '\0';

	}


}






///扩大空间留的都有余量，所以向后加一点\t之类的不需要再检查了
///注意这个数是
inline void output_to_buffer_int(Output_buffer_sub_block* sub_block, bitmapper_bs_iter input)
{


	int input_int_length;
	//2^64 -1 = 18446744073709551615, 占了20位
	///如果空间足够的,则size不用更新
	if (sub_block->length + INT_LENGTH + OFFSET< sub_block->size)
	{
		sprintf(sub_block->buffer + sub_block->length, "%llu", input);
		input_int_length = strlen(sub_block->buffer + sub_block->length);
		sub_block->length = sub_block->length + input_int_length;
	}
	else
	{
		///空间不够，size需要扩大
		sub_block->size = (sub_block->length + INT_LENGTH + OFFSET) * 2 + 1;
		char* tmp = (char*)malloc(sub_block->size);
		memcpy(tmp, sub_block->buffer, sub_block->length);
		free(sub_block->buffer);
		sub_block->buffer = tmp;

		sprintf(sub_block->buffer + sub_block->length, "%llu", input);
		input_int_length = strlen(sub_block->buffer + sub_block->length);
		sub_block->length = sub_block->length + input_int_length;

	}


	
}

int init_output_methy(char *fileName);


void output_single_methy_CpG(bitmapper_bs_iter tmp_pos, char* chrome_name, int nmethyl, int total);

void output_single_methy_CHG(bitmapper_bs_iter tmp_pos, char* chrome_name, int nmethyl, int total);

void output_single_methy_CHH(bitmapper_bs_iter tmp_pos, char* chrome_name, int nmethyl, int total);



#endif
