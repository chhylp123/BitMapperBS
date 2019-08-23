/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
///#include <zlib.h>
#include "Auxiliary.h"
#include "Process_Reads.h"
#include "Index.h"
#include <pthread.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>



pthread_mutex_t i_readinputMutex;
pthread_mutex_t i_queueMutex;
pthread_mutex_t i_terminateMutex;
pthread_cond_t i_flushCond;
pthread_cond_t i_readinputflushCond;
pthread_cond_t i_stallCond;
pthread_cond_t i_readinputstallCond;
pthread_mutex_t i_doneMutex;

FILE *_r_fp1;
FILE *_r_fp2;

gzFile g_r_fp1;
gzFile g_r_fp2;

Read *_r_seq;
int _r_seqCnt;

int fastq_init;

char rc_table[128];

Read_buffer_pe buffer_pe;
Read_buffer_single buffer_single;

Read_buffer_pe_sub_block tmp_sub_block;
Read_buffer_single_sub_block tmp_sub_block_single;

char tmp_read_buffer1[SEQ_MAX_LENGTH];
char tmp_read_buffer2[SEQ_MAX_LENGTH];

int file_format;
char flag_array[2];





inline char *Input_line_from_file1( char *seq )
{
	if (file_format == FASTQ)
	{
		return fgets(seq, SEQ_MAX_LENGTH, _r_fp1);
	}
	else
	{
		seq = gzgets(g_r_fp1, seq, SEQ_MAX_LENGTH);

		return (!gzeof(g_r_fp1)) ? flag_array : NULL;

	}
  
}

inline char *Input_line_from_file2(char *seq)
{
	if (file_format == FASTQ)
	{
		return fgets(seq, SEQ_MAX_LENGTH, _r_fp2);
	}
	else
	{
		seq = gzgets(g_r_fp2, seq, SEQ_MAX_LENGTH);

		return (!gzeof(g_r_fp2)) ? flag_array : NULL;
	}
}









inline void output_to_read_name(Read* current_read, char* input)
{
	int name_length = strlen(input);


	///如果空间足够的,则size不用更新
	if (name_length < current_read->name_size)
	{
		memcpy(current_read->name, input, name_length);
		current_read->name[name_length] = '\0';
	}
	else
	{
		///空间不够，size需要扩大
		current_read->name_size = name_length + 1;
		free(current_read->name);
		current_read->name = (char *)malloc(current_read->name_size);


		memcpy(current_read->name, input, name_length);
		current_read->name[name_length] = '\0';

	}

}



inline void resize_read_seq_req_qual(Read* current_read, int seq_length)
{


	///如果空间足够的,则size不用更新
	if (seq_length >= current_read->seq_size)
	{
		///空间不够，size需要扩大
		current_read->seq_size = seq_length + 1;
		free(current_read->qual);
		free(current_read->seq);
		free(current_read->rseq);


		current_read->qual = (char *)malloc(current_read->seq_size);
		current_read->seq = (char *)malloc(current_read->seq_size);
		current_read->rseq = (char *)malloc(current_read->seq_size);
	}

}






/**********************************************/
int inputReads_paired_directly(
	Read *seqList1,
	Read *seqList2
	)
{

	int read1_length;
	int read2_length;
	int j;

	////读入read1和read2的名称
	if (Input_line_from_file1(tmp_read_buffer1) && Input_line_from_file2(tmp_read_buffer2))
	{


		//****************************读第一个read*************************************
		///写到read1的名称中去
		output_to_read_name(seqList1, tmp_read_buffer1);



		Input_line_from_file1(tmp_read_buffer1);
		read1_length = strlen(tmp_read_buffer1);
		resize_read_seq_req_qual(seqList1, read1_length);

		///因为read本身被读到tmp_read_buffer1中去了
		memcpy(seqList1->seq, tmp_read_buffer1, read1_length);
		seqList1->seq[read1_length] = '\0';
		///第三行
		//Input_line_from_file1(seqList1->qual);  
		Input_line_from_file1(tmp_read_buffer1);
		///第四行
		Input_line_from_file1(seqList1->qual);
		//****************************读第一个read*************************************

		/*************************处理第一个read********************************/
		///处理read seq
		read1_length = strlen(seqList1->seq);
		///去回车
		seqList1->seq[read1_length - 1] = '\0';
		read1_length = read1_length - 1;

		for (j = 0; j < read1_length; j++)
		{
			seqList1->seq[j] = toupper(seqList1->seq[j]);
			///这是反向互补
			seqList1->rseq[read1_length - j - 1] = rc_table[seqList1->seq[j]];
		}
		seqList1->rseq[read1_length] = '\0';


		///处理read qual
		///去回车
		seqList1->qual[read1_length] = '\0';



		///curr_sub_block->read1[i].hits[0] = 0;
		seqList1->length = read1_length;

		/*************************处理第一个read********************************/
















		//****************************读第二个read*************************************
		///写到read2的名称中去
		output_to_read_name(seqList2, tmp_read_buffer2);


		Input_line_from_file2(tmp_read_buffer2);
		read2_length = strlen(tmp_read_buffer2);
		resize_read_seq_req_qual(seqList2, read2_length);

		///因为read本身被读到tmp_read_buffer2中去了
		memcpy(seqList2->rseq, tmp_read_buffer2, read2_length);
		seqList2->rseq[read2_length] = '\0';


		///第三行    ///
		///Input_line_from_file2(seqList2->qual);
		Input_line_from_file2(tmp_read_buffer2);
		///第四行
		Input_line_from_file2(seqList2->qual);
		//****************************读第二个read*************************************



		/*************************处理第二个read********************************/
		///处理read seq
		read2_length = strlen(seqList2->rseq);
		///去回车
		seqList2->rseq[read2_length - 1] = '\0';
		read2_length = read2_length - 1;

		for (j = 0; j < read2_length; j++)
		{
			seqList2->rseq[j] = toupper(seqList2->rseq[j]);
			seqList2->seq[read2_length - j - 1] = rc_table[seqList2->rseq[j]];
		}
		seqList2->seq[read2_length] = '\0';


		///处理read qual
		///去回车
		seqList2->qual[read2_length] = '\0';



		///curr_sub_block->read2[i].hits[0] = 0;
		seqList2->length = read2_length;


		/*************************处理第二个read********************************/



		///处理read name
		read1_length = strlen(seqList1->name);
		///去回车
		seqList1->name[read1_length - 1] = '\0';

		///处理read name
		read2_length = strlen(seqList2->name);
		///去回车
		seqList2->name[read2_length - 1] = '\0';



		for (j = 0; j < read1_length && j < read2_length; j++)
		{
			if (seqList1->name[j] != seqList2->name[j]||
				seqList1->name[j] == ' ' || 
				seqList1->name[j] == '/'
				)
			{
				seqList1->name[j] = '\0';
				seqList2->name[j] = '\0';
				break;
			}
		}

		return 1;

	}
	else
	{
		return 0;
	}

}



inline int post_process_paired_reads(Read_buffer_pe_sub_block* curr_sub_block)
{
	long long i, j;
	int seq_length1, seq_length2;
	for (i = 0; i < curr_sub_block->sub_block_read_number; i++)
	{
		/*************************处理第一个read********************************/
		///处理read seq
		seq_length1 = strlen(curr_sub_block->read1[i].seq);
		///去回车
		curr_sub_block->read1[i].seq[seq_length1 - 1] = '\0';
		seq_length1 = seq_length1 - 1;

		for (j = 0; j < seq_length1; j++)
		{
			curr_sub_block->read1[i].seq[j] = toupper(curr_sub_block->read1[i].seq[j]);
			///这是反向互补
			curr_sub_block->read1[i].rseq[seq_length1 - j - 1] = rc_table[curr_sub_block->read1[i].seq[j]];
		}
		curr_sub_block->read1[i].rseq[seq_length1] = '\0';


		///处理read qual
		///去回车
		curr_sub_block->read1[i].qual[seq_length1] = '\0';



		///curr_sub_block->read1[i].hits[0] = 0;
		curr_sub_block->read1[i].length = seq_length1;

		/*************************处理第一个read********************************/




		/*************************处理第二个read********************************/
		///处理read seq
		seq_length2 = strlen(curr_sub_block->read2[i].rseq);
		///去回车
		curr_sub_block->read2[i].rseq[seq_length2 - 1] = '\0';
		seq_length2 = seq_length2 - 1;

		for (j = 0; j < seq_length2; j++)
		{
			curr_sub_block->read2[i].rseq[j] = toupper(curr_sub_block->read2[i].rseq[j]);
			curr_sub_block->read2[i].seq[seq_length2 - j - 1] = rc_table[curr_sub_block->read2[i].rseq[j]];
		}
		curr_sub_block->read2[i].seq[seq_length2] = '\0';


		///处理read qual
		///去回车
		curr_sub_block->read2[i].qual[seq_length2] = '\0';



		///curr_sub_block->read2[i].hits[0] = 0;
		curr_sub_block->read2[i].length = seq_length2;


		/*************************处理第二个read********************************/

		///处理read name
		seq_length1 = strlen(curr_sub_block->read1[i].name);
		///去回车
		curr_sub_block->read1[i].name[seq_length1 - 1] = '\0';
		///处理read name
		seq_length2 = strlen(curr_sub_block->read2[i].name);
		///去回车
		curr_sub_block->read2[i].name[seq_length2 - 1] = '\0';

		for (j = 0; j < seq_length1; j++)
		{
			if (
				curr_sub_block->read1[i].name[j] != curr_sub_block->read2[i].name[j]||
				curr_sub_block->read1[i].name[j] == ' ' || 
				curr_sub_block->read1[i].name[j] == '/'				
				)
			{
				curr_sub_block->read1[i].name[j] = '\0';
				curr_sub_block->read2[i].name[j] = '\0';
				break;
			}
		}

	}
}












inline int post_process_single_reads(Read_buffer_single_sub_block* curr_sub_block)
{
	long long i, j;
	int seq_length1;
	for (i = 0; i < curr_sub_block->sub_block_read_number; i++)
	{
		/*************************处理第一个read********************************/
		///处理read name
		seq_length1 = strlen(curr_sub_block->read[i].name);
		///去回车
		curr_sub_block->read[i].name[seq_length1 - 1] = '\0';
		for (j = 0; j < seq_length1; j++)
		{
			if (curr_sub_block->read[i].name[j] == ' ' || curr_sub_block->read[i].name[j] == '/')
			{
				curr_sub_block->read[i].name[j] = '\0';
				break;
			}
		}



		///处理read seq
		seq_length1 = strlen(curr_sub_block->read[i].seq);
		///去回车
		curr_sub_block->read[i].seq[seq_length1 - 1] = '\0';
		seq_length1 = seq_length1 - 1;

		for (j = 0; j < seq_length1; j++)
		{
			curr_sub_block->read[i].seq[j] = toupper(curr_sub_block->read[i].seq[j]);
			///这是反向互补
			curr_sub_block->read[i].rseq[seq_length1 - j - 1] = rc_table[curr_sub_block->read[i].seq[j]];
		}
		curr_sub_block->read[i].rseq[seq_length1] = '\0';


		///处理read qual
		///去回车
		curr_sub_block->read[i].qual[seq_length1] = '\0';



		///curr_sub_block->read1[i].hits[0] = 0;
		curr_sub_block->read[i].length = seq_length1;

		/*************************处理第一个read********************************/



	}
}



inline int post_process_single_reads_pbat(Read_buffer_single_sub_block* curr_sub_block)
{
	long long i, j;
	int seq_length2;
	for (i = 0; i < curr_sub_block->sub_block_read_number; i++)
	{


		/*************************处理第二个read********************************/
		///处理read name
		seq_length2 = strlen(curr_sub_block->read[i].name);
		///去回车
		curr_sub_block->read[i].name[seq_length2 - 1] = '\0';
		for (j = 0; j < seq_length2; j++)
		{
			if (curr_sub_block->read[i].name[j] == ' ' || curr_sub_block->read[i].name[j] == '/')
			{
				curr_sub_block->read[i].name[j] = '\0';
				break;
			}
		}



		///处理read seq
		seq_length2 = strlen(curr_sub_block->read[i].rseq);
		///去回车
		curr_sub_block->read[i].rseq[seq_length2 - 1] = '\0';
		seq_length2 = seq_length2 - 1;

		for (j = 0; j < seq_length2; j++)
		{
			curr_sub_block->read[i].rseq[j] = toupper(curr_sub_block->read[i].rseq[j]);
			curr_sub_block->read[i].seq[seq_length2 - j - 1] = rc_table[curr_sub_block->read[i].rseq[j]];
		}
		curr_sub_block->read[i].seq[seq_length2] = '\0';


		///处理read qual
		///去回车
		curr_sub_block->read[i].qual[seq_length2] = '\0';



		///curr_sub_block->read2[i].hits[0] = 0;
		curr_sub_block->read[i].length = seq_length2;


		/*************************处理第二个read********************************/





	}
}




/**********************************************/
inline int inputReads_paired_pure(
	Read *seqList1,
	Read *seqList2
	)
{
	int read1_length;
	int read2_length;

	////读入read1和read2的名称
	if (Input_line_from_file1(tmp_read_buffer1) && Input_line_from_file2(tmp_read_buffer2))
	{


		//****************************读第一个read*************************************
		///写到read1的名称中去
		output_to_read_name(seqList1, tmp_read_buffer1);

		

		Input_line_from_file1(tmp_read_buffer1);
		read1_length = strlen(tmp_read_buffer1);
		resize_read_seq_req_qual(seqList1, read1_length);

		///因为read本身被读到tmp_read_buffer1中去了
		memcpy(seqList1->seq, tmp_read_buffer1, read1_length);
		seqList1->seq[read1_length] = '\0';
		///第三行
		//Input_line_from_file1(seqList1->qual);  
		Input_line_from_file1(tmp_read_buffer1);
		///第四行
		Input_line_from_file1(seqList1->qual);
		//****************************读第一个read*************************************


		//****************************读第二个read*************************************
		///写到read2的名称中去
		output_to_read_name(seqList2, tmp_read_buffer2);


		Input_line_from_file2(tmp_read_buffer2);
		read2_length = strlen(tmp_read_buffer2);
		resize_read_seq_req_qual(seqList2, read2_length);

		///因为read本身被读到tmp_read_buffer2中去了
		memcpy(seqList2->rseq, tmp_read_buffer2, read2_length);
		seqList2->rseq[read2_length] = '\0';


		///第三行    ///
		///Input_line_from_file2(seqList2->qual);
		Input_line_from_file2(tmp_read_buffer2);
		///第四行
		Input_line_from_file2(seqList2->qual);
		//****************************读第二个read*************************************
		return 1;

	}
	else
	{
		return 0;
	}

}





/**********************************************/
inline int inputReads_single_pure(
	Read *seqList1
	)
{
	int read1_length;

	////读入read1和read2的名称
	if (Input_line_from_file1(tmp_read_buffer1))
	{


		//****************************读第一个read*************************************
		///写到read1的名称中去
		output_to_read_name(seqList1, tmp_read_buffer1);



		Input_line_from_file1(tmp_read_buffer1);
		read1_length = strlen(tmp_read_buffer1);
		resize_read_seq_req_qual(seqList1, read1_length);

		///因为read本身被读到tmp_read_buffer1中去了
		memcpy(seqList1->seq, tmp_read_buffer1, read1_length);
		seqList1->seq[read1_length] = '\0';
		///第三行
		//Input_line_from_file1(seqList1->qual);  
		Input_line_from_file1(tmp_read_buffer1);
		///第四行
		Input_line_from_file1(seqList1->qual);
		//****************************读第一个read*************************************


		return 1;

	}
	else
	{
		return 0;
	}

}








/**********************************************/
inline int inputReads_single_pure_pbat(
	Read *seqList2
	)
{
	int read2_length;

	////读入read1和read2的名称
	if (Input_line_from_file2(tmp_read_buffer2))
	{



		//****************************读第二个read*************************************
		///写到read2的名称中去
		output_to_read_name(seqList2, tmp_read_buffer2);


		Input_line_from_file2(tmp_read_buffer2);
		read2_length = strlen(tmp_read_buffer2);
		resize_read_seq_req_qual(seqList2, read2_length);

		///因为read本身被读到tmp_read_buffer2中去了
		memcpy(seqList2->rseq, tmp_read_buffer2, read2_length);
		seqList2->rseq[read2_length] = '\0';


		///第三行    ///
		///Input_line_from_file2(seqList2->qual);
		Input_line_from_file2(tmp_read_buffer2);
		///第四行
		Input_line_from_file2(seqList2->qual);
		//****************************读第二个read*************************************
		return 1;

	}
	else
	{
		return 0;
	}

}





















/**********************************************/
int inputReads_single_directly_back(
	Read *seqList1
	)
{

	int nCnt;
	int i;
	int seq_length;

	if (Input_line_from_file1(seqList1->name))
	{

		///这一块是处理read的名称
		seq_length = strlen(seqList1->name);
		///这个就是去回车
		seqList1->name[seq_length - 1] = '\0';

		for (i = 0; i<seq_length; i++)
		{
			if (seqList1->name[i] == ' ')
			{
				seqList1->name[i] = '\0';
				break;
			}
		}




		//dummy纯粹就是个打杂的...，暂存用的
		//读质量值什么的
		Input_line_from_file1(seqList1->seq);

		if (seqList1->name[0] != '>')
		{
			Input_line_from_file1(seqList1->qual);
			Input_line_from_file1(seqList1->qual);
			///这个也是去回车
			seqList1->qual[strlen(seqList1->qual) - 1] = '\0';
		}



		//nCnt统计N的个数，
		//这仅统计了N的个数而已，并没用对其做任何处理
		//同时确定read的边界
		nCnt = 0;
		seq_length = strlen(seqList1->seq);
		///去回车
		seqList1->seq[seq_length - 1] = '\0';
		seq_length = seq_length - 1;

		for (i = 0; i < seq_length; i++)
		{
			seqList1->seq[i] = toupper(seqList1->seq[i]);

			if (seqList1->seq[i] == 'N')
			{
				nCnt++;
			}
			
			seqList1->rseq[seq_length - i - 1] = rc_table[seqList1->seq[i]];
			///seqList1->rseq[seq_length - i - 1] = seqList1->seq[i];
		}
		seqList1->rseq[seq_length] = '\0';

		if (nCnt > seq_length*0.7)
		{
			return 3;
		}

		///seqList1->hits[0] = 0;
		seqList1->length = seq_length;


		return 1;
	}
	else
	{
		return 0;
	}



}





/**********************************************/
int inputReads_single_directly(
	Read *seqList1
	)
{

	int read1_length;
	int j;

	if (Input_line_from_file1(tmp_read_buffer1))
	{
		//****************************读第一个read*************************************
		///写到read1的名称中去
		output_to_read_name(seqList1, tmp_read_buffer1);



		Input_line_from_file1(tmp_read_buffer1);
		read1_length = strlen(tmp_read_buffer1);
		resize_read_seq_req_qual(seqList1, read1_length);

		///因为read本身被读到tmp_read_buffer1中去了
		memcpy(seqList1->seq, tmp_read_buffer1, read1_length);
		seqList1->seq[read1_length] = '\0';
		///第三行
		//Input_line_from_file1(seqList1->qual);  
		Input_line_from_file1(tmp_read_buffer1);
		///第四行
		Input_line_from_file1(seqList1->qual);
		//****************************读第一个read*************************************

		/*************************处理第一个read********************************/
		///处理read name
		read1_length = strlen(seqList1->name);
		///去回车
		seqList1->name[read1_length - 1] = '\0';
		for (j = 0; j < read1_length; j++)
		{
			if (seqList1->name[j] == ' ' || seqList1->name[j] == '/')
			{
				seqList1->name[j] = '\0';
				break;
			}
		}



		///处理read seq
		read1_length = strlen(seqList1->seq);
		///去回车
		seqList1->seq[read1_length - 1] = '\0';
		read1_length = read1_length - 1;

		for (j = 0; j < read1_length; j++)
		{
			seqList1->seq[j] = toupper(seqList1->seq[j]);
			///这是反向互补
			seqList1->rseq[read1_length - j - 1] = rc_table[seqList1->seq[j]];
		}
		seqList1->rseq[read1_length] = '\0';


		///处理read qual
		///去回车
		seqList1->qual[read1_length] = '\0';



		///curr_sub_block->read1[i].hits[0] = 0;
		seqList1->length = read1_length;

		/*************************处理第一个read********************************/

		return 1;
	}
	else
	{
		return 0;
	}


}




int inputReads_single_directly_pbat_back(Read *seqList1)
{

	int nCnt;
	int i;
	int seq_length;

	if (Input_line_from_file1(seqList1->name))
	{

		///这一块是处理read的名称
		seq_length = strlen(seqList1->name);
		///这个就是去回车
		seqList1->name[seq_length - 1] = '\0';

		for (i = 0; i<seq_length; i++)
		{
			if (seqList1->name[i] == ' ')
			{
				seqList1->name[i] = '\0';
				break;
			}
		}


		//dummy纯粹就是个打杂的...，暂存用的
		//读质量值什么的
		////Input_line_from_file1(seqList1->seq);
		Input_line_from_file1(seqList1->rseq);

		if (seqList1->name[0] != '>')
		{
			Input_line_from_file1(seqList1->qual);
			Input_line_from_file1(seqList1->qual);
			///这个也是去回车
			seqList1->qual[strlen(seqList1->qual) - 1] = '\0';
		}



		//nCnt统计N的个数，
		//这仅统计了N的个数而已，并没用对其做任何处理
		//同时确定read的边界
		nCnt = 0;
		///seq_length = strlen(seqList1->seq);
		seq_length = strlen(seqList1->rseq);
		///去回车
		///seqList1->seq[seq_length - 1] = '\0';
		seqList1->rseq[seq_length - 1] = '\0';
		seq_length = seq_length - 1;

		for (i = 0; i < seq_length; i++)
		{
			///seqList1->seq[i] = toupper(seqList1->seq[i]);
			seqList1->rseq[i] = toupper(seqList1->rseq[i]);

			///if (seqList1->seq[i] == 'N')
			if (seqList1->rseq[i] == 'N')
			{
				nCnt++;
			}

			///seqList1->rseq[seq_length - i - 1] = rc_table[seqList1->seq[i]];
			seqList1->seq[seq_length - i - 1] = rc_table[seqList1->rseq[i]];
		}
		///seqList1->rseq[seq_length] = '\0';
		seqList1->seq[seq_length] = '\0';

		if (nCnt > seq_length*0.7)
		{
			return 3;
		}

		///seqList1->hits[0] = 0;
		seqList1->length = seq_length;

		return 1;
	}
	else
	{
		return 0;
	}
}








int inputReads_single_directly_pbat(Read *seqList2)
{

	int read2_length;
	int j;



	////读入read1和read2的名称
	if (Input_line_from_file2(tmp_read_buffer2))
	{




		//****************************读第二个read*************************************
		///写到read2的名称中去
		output_to_read_name(seqList2, tmp_read_buffer2);


		Input_line_from_file2(tmp_read_buffer2);
		read2_length = strlen(tmp_read_buffer2);
		resize_read_seq_req_qual(seqList2, read2_length);

		///因为read本身被读到tmp_read_buffer2中去了
		memcpy(seqList2->rseq, tmp_read_buffer2, read2_length);
		seqList2->rseq[read2_length] = '\0';


		///第三行    ///
		///Input_line_from_file2(seqList2->qual);
		Input_line_from_file2(tmp_read_buffer2);
		///第四行
		Input_line_from_file2(seqList2->qual);
		//****************************读第二个read*************************************



		/*************************处理第二个read********************************/
		///处理read name
		read2_length = strlen(seqList2->name);
		///去回车
		seqList2->name[read2_length - 1] = '\0';
		for (j = 0; j < read2_length; j++)
		{
			
			if (seqList2->name[j] == ' ' || seqList2->name[j] == '/')
			{
				seqList2->name[j] = '\0';
				break;
			}
		}



		///处理read seq
		read2_length = strlen(seqList2->rseq);
		///去回车
		seqList2->rseq[read2_length - 1] = '\0';
		read2_length = read2_length - 1;

		for (j = 0; j < read2_length; j++)
		{
			seqList2->rseq[j] = toupper(seqList2->rseq[j]);
			seqList2->seq[read2_length - j - 1] = rc_table[seqList2->rseq[j]];
		}
		seqList2->seq[read2_length] = '\0';


		///处理read qual
		///去回车
		seqList2->qual[read2_length] = '\0';



		///curr_sub_block->read2[i].hits[0] = 0;
		seqList2->length = read2_length;


		/*************************处理第二个read********************************/




		return 1;

	}
	else
	{
		return 0;
	}

}












/**********************************************/
int inputReads_single(
		 Read *seqList1
		 )
{

	 char seq1[SEQ_MAX_LENGTH];
	 char rseq1[SEQ_MAX_LENGTH];
	 char name1[SEQ_MAX_LENGTH];
	 char qual1[SEQ_MAX_LENGTH];
	 char dummy[SEQ_MAX_LENGTH];
	 char ch;
	 int err1, err2;
	 int nCnt;
	 int maxCnt = 0;
	 int i;
	 int clipped = 0;
	 int count_reads_avli=0;
	 int seq_length;

	if(Input_line_from_file1(name1))
	{

		///这一块是处理read的名称
		err1 = 0;
		seq_length = strlen(name1);
		///这个就是去回车
		name1[seq_length - 1] = '\0';

		for (i = 0; i<seq_length; i++)
		{
			if (name1[i] == ' ')
			{
				name1[i] = '\0';
				break;
			}
		}




		//dummy纯粹就是个打杂的...，暂存用的
		//读质量值什么的
		Input_line_from_file1(seq1);

                if(name1[0]!='>')
                {
		      Input_line_from_file1(dummy);
		      Input_line_from_file1(qual1);
		      ///这个也是去回车
		      qual1[strlen(qual1) - 1] = '\0';
                }



		//nCnt统计N的个数，
		//这仅统计了N的个数而已，并没用对其做任何处理
		//同时确定read的边界
		nCnt = 0;
		seq_length = strlen(seq1);
		///去回车
		seq1[seq_length - 1] = '\0';
		seq_length = seq_length - 1;

		for (i = 0; i < seq_length; i++)
		{
		  seq1[i] = toupper (seq1[i]);

		  if (seq1[i] == 'N')
		  {
			  nCnt++;
		  }
		}
                


		///if (nCnt <=thread_e)
                if (nCnt <=seq_length*0.7)
		{


		  seq_length = strlen(seq1);
		  //这是read编号吧，确定这是当前第几个read
		  //rseq1这个字符数组中存储了反向互补链。既要反向，又要互补
		  ///reverseComplement(seq1, rseq1, seq_length);
		  //rseq1这个字符数组中存储了反向链, 只反向，不互补
		  reverse_pattern(seq1, rseq1, seq_length);

		  int i;
		  //hits就分配了一个字符啊，貌似仅有这一个hits[0]吧
		  //这个应该是存这个read匹配上了多少个位置，这里当然置为0了
		 /// seqList1->hits[0] = 0;

		  //这是把终止符也复制进去了的节奏

		  strcpy(seqList1->seq, seq1);
		  strcpy(seqList1->rseq, rseq1);
		  strcpy(seqList1->qual, qual1);
		  strcpy(seqList1->name, name1);
		  seqList1->length = seq_length;
		}
		else
		{
    		  //这个应该是N太多了，超过阈值了，那这个read就直接丢弃了
			return 3;
		}

			return 1;
		}
		else
		{
			return 0;
		}



}



/**********************************************/
int inputReads_paired(
		 Read *seqList1,
		 Read *seqList2
		 )
{
  char seq1[SEQ_MAX_LENGTH];
  char rseq1[SEQ_MAX_LENGTH];
  char name1[SEQ_MAX_LENGTH];
  char qual1[SEQ_MAX_LENGTH];
  char seq2[SEQ_MAX_LENGTH];
  char rseq2[SEQ_MAX_LENGTH];
  char name2[SEQ_MAX_LENGTH];
  char qual2[SEQ_MAX_LENGTH];
  char dummy[SEQ_MAX_LENGTH];
  char ch;
  int err1, err2;
  int nCnt;
  int maxCnt = 0;
  int i;
  int clipped = 0;
  int count_reads_avli=0;

    if( Input_line_from_file1(name1) )
    {
      err1 = 0;

      err2=0;

      Input_line_from_file1(seq1);

      name1[strlen(name1)-1] = '\0';

      for (i=0; i<strlen(name1);i++)
     {
          if (name1[i] == ' ')
          {
              name1[i] = '\0';
              break;
          }
     }




      //这个应该是读的质量值什么的
      if ( fastq_init )
	{
        //dummy纯粹就是个打杂的...，暂存用的
          Input_line_from_file1(dummy);
          Input_line_from_file1(qual1);
          qual1[strlen(qual1)-1] = '\0';
	}
      else
	{
        //如果是fasta格式，直接向qual1里存一个*
          sprintf(qual1, "*");
	}

      if (cropSize > 0)
	{
	  seq1[cropSize] = '\0';
	  //如果有质量值把质量值也裁剪
	  if ( fastq_init )
	    qual1[cropSize] = '\0';
	}






	 Input_line_from_file2(name2);
	  Input_line_from_file2(seq2);
	  name2[strlen(name2)-1] = '\0';
	  for (i=0; i<strlen(name2);i++)
	    {
	      if (name2[i] == ' ')
		{
		  name2[i] = '\0';
		  break;
		}

	    }

	  if ( fastq_init )
	    {
	      Input_line_from_file2(dummy);
	      Input_line_from_file2(qual2);

	      qual2[strlen(qual2)-1] = '\0';
	    }
	  else
	    {
	      sprintf(qual2, "*");
	    }


	  if (cropSize > 0)
	    {
	      seq2[cropSize] = '\0';
	      if ( fastq_init )
		qual2[cropSize] = '\0';
	    }



      nCnt = 0;
      for (i=0; i<strlen(seq1); i++)
	{
	  seq1[i] = toupper (seq1[i]);
	  if (seq1[i] == 'N')
	    {
	      nCnt++;
	    }
	  else if (isspace(seq1[i]))
	    {
	      seq1[i] = '\0';
	      break;
	    }
	}

	 if (nCnt > thread_e)
	    {
	      err1 = 1;
	    }




	  nCnt = 0;
	  for (i=0; i<strlen(seq2); i++)
	    {
	      seq2[i] = toupper (seq2[i]);
	      if (seq2[i] == 'N')
		{
		  nCnt++;

		}
	      else if (isspace(seq2[i]))
		{
		  seq2[i] = '\0';
		}
	    }


	  if (nCnt > thread_e)
	    {
	      err2 = 1;
	    }


	  if (strlen(seq1) < strlen(seq2)) {
	    seq2[strlen(seq1)] = '\0';
	    if ( fastq_init )
	      qual2[strlen(seq1)] = '\0';
	    if (!clipped) clipped = 2;
	  }
	  else if (strlen(seq1) > strlen(seq2)){
	    seq1[strlen(seq2)] = '\0';
	    if ( fastq_init )
	      qual1[strlen(seq2)] = '\0';
	    if (!clipped) clipped = 1;
	  }

	  if (clipped == 1 || clipped == 2){
	    fprintf(stderr, "[PE mode Warning] Sequence lengths are different,  read #%d is clipped to match.\n", clipped);
	    clipped = 3;
	    return 3;
	  }








      ///err2是对read2而言的
      if (err1==0 &&err2==0)
	{


	  int tmplen = strlen(name1);
	  if (strcmp(name1, name2) != 0)
	    {
	      tmplen = strlen(name1)-2;
	    }

	  int _mtmp = strlen(seq1);


	  reverseComplement(seq1, rseq1, _mtmp);
	  int i;

	  ///seqList1->hits[0] = 0;

	  for (i=0; i<=_mtmp; i++)
	    {
	      seqList1->seq[i] = seq1[i];
	      seqList1->rseq[i] = rseq1[i] ;
	      seqList1->qual[i] = qual1[i];
	    }


	  name1[tmplen]='\0';
	  seqList1->rseq[_mtmp]=seqList1->qual[_mtmp]='\0';
	  sprintf(seqList1->name,"%s%c", ((char*)name1)+1,'\0');


	  reverseComplement(seq2, rseq2, _mtmp);

	  ///seqList2->hits[0] = 0;

	  for (i=0; i<=_mtmp; i++)
	    {
	      seqList2->seq[i] = seq2[i];
	      seqList2->rseq[i] = rseq2[i] ;
	      seqList2->qual[i] = qual2[i];
	    }


	  name2[tmplen]='\0';
	  seqList2->rseq[_mtmp]=seqList2->qual[_mtmp]='\0';
	  sprintf(seqList2->name,"%s%c", ((char*)name2)+1,'\0');

	}
    else
	{
    	  //这个应该是N太多了，超过阈值了，那这个read就直接丢弃了
        return 3;
	}

        return 1;
    }
    else
    {
        return 0;
    }



}

int check_file_format(char *fileName)
{
	int i = 0;

	char* file_extention;
	i = strlen(fileName) - 1;
	while (i >= 0 && fileName[i] != '.')
	{
		i--;
	}

	if (i < 0)
	{
		fprintf(stderr, 
			"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
		exit(0);
	}

	if ((strcmp(".fq", fileName+i) == 0)
		||
		strcmp(".fastq", fileName + i) == 0)
	{
		return FASTQ;
	}
	else if (strcmp(".gz", fileName + i) == 0)
	{
		i--;
		while (i >= 0 && fileName[i] != '.')
		{
			i--;
		}


		if (i < 0)
		{
			fprintf(stderr,
				"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
			exit(0);
		}

		if ((strcmp(".fq.gz", fileName + i) == 0)
			||
			strcmp(".fastq.gz", fileName + i) == 0)
		{
			return FASTQGZ;
		}
		else
		{
			fprintf(stderr,
				"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
			exit(0);
		}

	}


	fprintf(stderr,
		"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
	exit(0);
}


/**********************************************/
int initiReadAllReads(char *fileName1,
		 char *fileName2,
		 unsigned char pairedEnd,
		 int* format
		 )
{

	int format1, format2;
	format1 = check_file_format(fileName1);
	
	if (format1 == FASTQ)
	{
		///打开read1
		_r_fp1 = fopen(fileName1, "r");

		if (_r_fp1 == NULL)
		{
			return 0;
		}
	}
	else if (format1 == FASTQGZ)
	{
		int fd = open(fileName1, O_CREAT | O_RDONLY, 0666);
		g_r_fp1 = gzdopen(fd, "rb");
	}
	

    


	///处理read2
    if ( pairedEnd && fileName2 != NULL )
	{
	  format2 = check_file_format(fileName2);


	  if (format2 == FASTQ)
	  {

		  _r_fp2 = fopen(fileName2, "r");
		  if (_r_fp2 == NULL)
		  {
			  return 0;
		  }
	  }
	  else if (format2 == FASTQGZ)
	  {
		  int fd = open(fileName2, O_CREAT | O_RDONLY, 0666);
		  g_r_fp2 = gzdopen(fd, "rb");
	  }


	}
    else
	{
	  format2 = format1;
	  _r_fp2 = _r_fp1;
	  g_r_fp2 = g_r_fp1;
	}



	if (format2 == format1 && format2 > 0)
	{
		if (format2 == FASTQ)
		{
			fprintf(stderr,
				" Read files are in FASTQ format...\n");
		}

		if (format2 == FASTQGZ)
		{
			fprintf(stderr,
				" Read files are in compressed FASTQ format...\n");
		}
		
	}
	else
	{
		fprintf(stderr,
			"ERROR: the input read files must be FASTQ format (file extension: .fq or .fastq) or compressed FASTQ format (file extension: .fq.gz or .fastq.gz)!\n");
		exit(0);
	}

	
	for (size_t i = 0; i < 128; i++)
	{
		rc_table[i] = i;
	}
	


	rc_table['A'] = 'T';
	rc_table['T'] = 'A';
	rc_table['C'] = 'G';
	rc_table['G'] = 'C';

	*format = format1;
	file_format = format1;

	return 1;
}







/**********************************************/
int exchange_two_reads()
{
	FILE *tmp;

	tmp = _r_fp1;
	_r_fp1 = _r_fp2;
	_r_fp2 = tmp;


	gzFile tmp_gz;
	tmp_gz = g_r_fp1;
	g_r_fp1 = g_r_fp2;
	g_r_fp2 = tmp_gz;

	return 1;
}







inline void load_sub_block_paired_read(Read* read_batch1, Read* read_batch2, int batch_read_size,
	int* return_file_flag, long long* obtain_reads_num)
{
	int inner_i = 0;
	int file_flag = 1;




	while (inner_i<batch_read_size)
	{
		///file_flag = inputReads_single_directly(&read_batch[inner_i]);



		///file_flag = inputReads_paired_directly(&read_batch1[inner_i], &read_batch2[inner_i]);
		///这里改了
		file_flag = inputReads_paired_pure(&read_batch1[inner_i], &read_batch2[inner_i]);

		



		if (file_flag == 1)
		{
			inner_i++;
		}
		else if (file_flag == 0)
		{
			break;
		}
		///file_flag还有可能为3，为3就丢弃了
	}

	if (!(inner_i == 0 && file_flag == 0))
	{
		file_flag = 1;
	}

	*return_file_flag = file_flag;
	*obtain_reads_num = inner_i;

}





inline void load_sub_block_single_read(Read* read_batch1, int batch_read_size,
	int* return_file_flag, long long* obtain_reads_num)
{
	int inner_i = 0;
	int file_flag = 1;




	while (inner_i<batch_read_size)
	{
		///file_flag = inputReads_single_directly(&read_batch[inner_i]);



		///这里改了
		///file_flag = inputReads_paired_pure(&read_batch1[inner_i], &read_batch2[inner_i]);
		file_flag = inputReads_single_pure(&read_batch1[inner_i]);
		




		if (file_flag == 1)
		{
			inner_i++;
		}
		else if (file_flag == 0)
		{
			break;
		}
		///file_flag还有可能为3，为3就丢弃了
	}

	if (!(inner_i == 0 && file_flag == 0))
	{
		file_flag = 1;
	}

	*return_file_flag = file_flag;
	*obtain_reads_num = inner_i;

}






inline void load_sub_block_single_read_pbat(Read* read_batch1, int batch_read_size,
	int* return_file_flag, long long* obtain_reads_num)
{
	int inner_i = 0;
	int file_flag = 1;




	while (inner_i<batch_read_size)
	{
		///file_flag = inputReads_single_directly(&read_batch[inner_i]);



		///这里改了
		///file_flag = inputReads_paired_pure(&read_batch1[inner_i], &read_batch2[inner_i]);
		file_flag = inputReads_single_pure_pbat(&read_batch1[inner_i]);





		if (file_flag == 1)
		{
			inner_i++;
		}
		else if (file_flag == 0)
		{
			break;
		}
		///file_flag还有可能为3，为3就丢弃了
	}

	if (!(inner_i == 0 && file_flag == 0))
	{
		file_flag = 1;
	}

	*return_file_flag = file_flag;
	*obtain_reads_num = inner_i;

}






///初始化缓冲区
void init_Pair_Seq_input_buffer(int thread_number)
{
	///if (file_format == FASTQ)
	{
		///每一个sub_block的大小
		sub_block_inner_size = sub_block_inner_size / 2;


		buffer_pe.sub_block_number = 0;
		///总共有多少个sub_block
		buffer_pe.sub_block_size = sub_block_number * thread_number;
		///每个sub_block内部有多少对reads
		buffer_pe.sub_block_inner_size = sub_block_inner_size;



		buffer_pe.sub_block = (Read_buffer_pe_sub_block*)
			malloc(sizeof(Read_buffer_pe_sub_block)*buffer_pe.sub_block_size);

		for (int i = 0; i < buffer_pe.sub_block_size; i++)
		{
			init_single_sub_block_pe(&buffer_pe.sub_block[i]);
		}

		buffer_pe.all_read_end = 0;
	}/**
	else
	{
		///每一个sub_block的大小
		sub_block_inner_size = sub_block_inner_size / 2;
		buffer_pe.sub_block_inner_size = sub_block_inner_size;
	}
	**/

	


}






///初始化缓冲区
void init_Single_Seq_input_buffer(int thread_number)
{
	///if (file_format == FASTQ)
	{
		///每一个sub_block的大小
		///单端是不是不用除以2
		///sub_block_inner_size = sub_block_inner_size / 2;
		sub_block_inner_size = sub_block_inner_size;
		///这个好像是指针，代表有多少有效的block
		buffer_single.sub_block_number = 0;
		///总共有多少个sub_block，这个是最大的block数量
		buffer_single.sub_block_size = sub_block_number * thread_number;
		///每个sub_block内部有多少对reads
		buffer_single.sub_block_inner_size = sub_block_inner_size;



		buffer_single.sub_block = (Read_buffer_single_sub_block*)
			malloc(sizeof(Read_buffer_single_sub_block)*buffer_single.sub_block_size);

		for (int i = 0; i < buffer_single.sub_block_size; i++)
		{
			init_single_sub_block_single(&buffer_single.sub_block[i]);
		}

		buffer_single.all_read_end = 0;
	}/**
	else
	{
		buffer_single.sub_block_inner_size = sub_block_inner_size;
	}
	**/


}



inline int if_full_pe()
{

	if (buffer_pe.sub_block_number >= buffer_pe.sub_block_size)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_full_single()
{

	if (buffer_single.sub_block_number >= buffer_single.sub_block_size)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_empty_pe()
{

	if (buffer_pe.sub_block_number == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_empty_single()
{

	if (buffer_single.sub_block_number == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline void push_back_sub_block_pe(Read_buffer_pe_sub_block* tmp_sub_block)
{
	/**
	fprintf(stderr, "(i)Producer: tmp_sub_block.sub_block_read_number: %lld\n",
	tmp_sub_block->sub_block_read_number);
	fflush(stderr);
	**/

	///只交换指针
	///记得要交换，直接赋值会内存泄露，造成野指针
	Read *k1, *k2;
	k1 = buffer_pe.sub_block[buffer_pe.sub_block_number].read1;
	k2 = buffer_pe.sub_block[buffer_pe.sub_block_number].read2;

	buffer_pe.sub_block[buffer_pe.sub_block_number].read1 = tmp_sub_block->read1;
	buffer_pe.sub_block[buffer_pe.sub_block_number].read2 = tmp_sub_block->read2;

	tmp_sub_block->read1 = k1;
	tmp_sub_block->read2 = k2;

	buffer_pe.sub_block[buffer_pe.sub_block_number].sub_block_read_number = tmp_sub_block->sub_block_read_number;
	tmp_sub_block->sub_block_read_number = 0;


	/**
	fprintf(stderr, "(i)Producer: buffer_pe.sub_block[%d].sub_block_read_number: %lld\n",
	buffer_pe.sub_block_number, buffer_pe.sub_block[buffer_pe.sub_block_number].sub_block_read_number);
	fflush(stderr);
	**/


	buffer_pe.sub_block_number++;
}




inline void push_back_sub_block_single(Read_buffer_single_sub_block* tmp_sub_block)
{
	

	///只交换指针
	///记得要交换，直接赋值会内存泄露，造成野指针
	///Read *k1, *k2;
	Read *k1;
	k1 = buffer_single.sub_block[buffer_single.sub_block_number].read;

	buffer_single.sub_block[buffer_single.sub_block_number].read = tmp_sub_block->read;

	tmp_sub_block->read = k1;

	buffer_single.sub_block[buffer_single.sub_block_number].sub_block_read_number = tmp_sub_block->sub_block_read_number;
	tmp_sub_block->sub_block_read_number = 0;


	/**
	fprintf(stderr, "(i)Producer: buffer_single.sub_block[%d].sub_block_read_number: %lld\n",
	buffer_single.sub_block_number, buffer_single.sub_block[buffer_single.sub_block_number].sub_block_read_number);
	fflush(stderr);
	**/


	buffer_single.sub_block_number++;
}




inline void pop_back_sub_block_pe(Read_buffer_pe_sub_block* curr_sub_block)
{
	buffer_pe.sub_block_number--;
	///只交换指针
	///记得要交换，直接赋值会内存泄露，造成野指针
	Read *k1, *k2;
	k1 = buffer_pe.sub_block[buffer_pe.sub_block_number].read1;
	k2 = buffer_pe.sub_block[buffer_pe.sub_block_number].read2;

	buffer_pe.sub_block[buffer_pe.sub_block_number].read1 = curr_sub_block->read1;
	buffer_pe.sub_block[buffer_pe.sub_block_number].read2 = curr_sub_block->read2;

	curr_sub_block->read1 = k1;
	curr_sub_block->read2 = k2;


	curr_sub_block->sub_block_read_number = buffer_pe.sub_block[buffer_pe.sub_block_number].sub_block_read_number;
	buffer_pe.sub_block[buffer_pe.sub_block_number].sub_block_read_number = 0;



}



inline void pop_back_sub_block_single(Read_buffer_single_sub_block* curr_sub_block)
{
	buffer_single.sub_block_number--;
	///只交换指针
	///记得要交换，直接赋值会内存泄露，造成野指针
	///Read *k1, *k2;
	Read *k1;
	k1 = buffer_single.sub_block[buffer_single.sub_block_number].read;
	///k2 = buffer_single.sub_block[buffer_single.sub_block_number].read2;

	buffer_single.sub_block[buffer_single.sub_block_number].read = curr_sub_block->read;
	///buffer_single.sub_block[buffer_single.sub_block_number].read2 = curr_sub_block->read2;

	curr_sub_block->read = k1;
	///curr_sub_block->read2 = k2;


	curr_sub_block->sub_block_read_number = buffer_single.sub_block[buffer_single.sub_block_number].sub_block_read_number;
	buffer_single.sub_block[buffer_single.sub_block_number].sub_block_read_number = 0;



}




void* input_pe_reads_muti_threads(void*)
{


	int i = 0;
	int file_flag = 1;

	///这是小的输入缓冲区
	/**
	tmp_sub_block.read1 = (Read*)malloc(sizeof(Read)*buffer_pe.sub_block_inner_size);
	tmp_sub_block.read2 = (Read*)malloc(sizeof(Read)*buffer_pe.sub_block_inner_size);
	tmp_sub_block.sub_block_read_number = 0;
	**/
	init_single_sub_block_pe(&tmp_sub_block);


	while (1)
	{





		///首先将一个block的reads放到本地的小缓冲区里
		load_sub_block_paired_read(tmp_sub_block.read1, tmp_sub_block.read2,
			buffer_pe.sub_block_inner_size, &file_flag, &tmp_sub_block.sub_block_read_number);


		///返回0就是读到文件末尾了
		if (file_flag == 0)
		{
			break;
		}


		pthread_mutex_lock(&i_readinputMutex);
		///如果满了，则该线程本身要wait，并通知消费者线程消费数据
		while (if_full_pe())
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&i_readinputstallCond);


			pthread_cond_wait(&i_readinputflushCond, &i_readinputMutex);
		}

		push_back_sub_block_pe(&tmp_sub_block);

		pthread_cond_signal(&i_readinputstallCond);
		pthread_mutex_unlock(&i_readinputMutex);
	}



	pthread_mutex_lock(&i_readinputMutex);
	buffer_pe.all_read_end = 1;
	pthread_cond_signal(&i_readinputstallCond);  //这行非常重要，必须要加！否则会陷入死锁！！！
	pthread_mutex_unlock(&i_readinputMutex);


}







void* input_single_reads_muti_threads(void*)
{


	int i = 0;
	int file_flag = 1;

	///这是小的输入缓冲区
	/**
	tmp_sub_block.read1 = (Read*)malloc(sizeof(Read)*buffer_pe.sub_block_inner_size);
	tmp_sub_block.read2 = (Read*)malloc(sizeof(Read)*buffer_pe.sub_block_inner_size);
	tmp_sub_block.sub_block_read_number = 0;
	**/
	init_single_sub_block_single(&tmp_sub_block_single);


	while (1)
	{





		///首先将一个block的reads放到本地的小缓冲区里
		/**
		load_sub_block_paired_read(tmp_sub_block.read1, tmp_sub_block.read2,
			buffer_pe.sub_block_inner_size, &file_flag, &tmp_sub_block.sub_block_read_number);**/
		load_sub_block_single_read(tmp_sub_block_single.read,
			buffer_single.sub_block_inner_size, &file_flag, &tmp_sub_block_single.sub_block_read_number);


		///返回0就是读到文件末尾了
		if (file_flag == 0)
		{
			break;
		}


		pthread_mutex_lock(&i_readinputMutex);
		///如果满了，则该线程本身要wait，并通知消费者线程消费数据
		while (if_full_single())
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&i_readinputstallCond);


			pthread_cond_wait(&i_readinputflushCond, &i_readinputMutex);
		}

		///push_back_sub_block_pe(&tmp_sub_block);
		push_back_sub_block_single(&tmp_sub_block_single);

		pthread_cond_signal(&i_readinputstallCond);
		pthread_mutex_unlock(&i_readinputMutex);
	}



	pthread_mutex_lock(&i_readinputMutex);
	buffer_single.all_read_end = 1;
	pthread_cond_signal(&i_readinputstallCond);  //这行非常重要，必须要加！否则会陷入死锁！！！
	pthread_mutex_unlock(&i_readinputMutex);


}







void* input_single_reads_muti_threads_pbat(void*)
{


	int i = 0;
	int file_flag = 1;


	init_single_sub_block_single(&tmp_sub_block_single);


	while (1)
	{





		///首先将一个block的reads放到本地的小缓冲区里
		/**
		load_sub_block_paired_read(tmp_sub_block.read1, tmp_sub_block.read2,
		buffer_pe.sub_block_inner_size, &file_flag, &tmp_sub_block.sub_block_read_number);**/
		load_sub_block_single_read_pbat(tmp_sub_block_single.read,
			buffer_single.sub_block_inner_size, &file_flag, &tmp_sub_block_single.sub_block_read_number);


		///返回0就是读到文件末尾了
		if (file_flag == 0)
		{
			break;
		}


		pthread_mutex_lock(&i_readinputMutex);
		///如果满了，则该线程本身要wait，并通知消费者线程消费数据
		while (if_full_single())
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&i_readinputstallCond);


			pthread_cond_wait(&i_readinputflushCond, &i_readinputMutex);
		}

		///push_back_sub_block_pe(&tmp_sub_block);
		push_back_sub_block_single(&tmp_sub_block_single);

		pthread_cond_signal(&i_readinputstallCond);
		pthread_mutex_unlock(&i_readinputMutex);
	}



	pthread_mutex_lock(&i_readinputMutex);
	buffer_single.all_read_end = 1;
	pthread_cond_signal(&i_readinputstallCond);  //这行非常重要，必须要加！否则会陷入死锁！！！
	pthread_mutex_unlock(&i_readinputMutex);


}






void init_single_sub_block_pe(Read_buffer_pe_sub_block* tmp_sub_block)
{
	tmp_sub_block->read1 = (Read*)malloc(sizeof(Read)*buffer_pe.sub_block_inner_size);
	tmp_sub_block->read2 = (Read*)malloc(sizeof(Read)*buffer_pe.sub_block_inner_size);




	for (int j = 0; j < buffer_pe.sub_block_inner_size; j++)
	{
		tmp_sub_block->read1[j].name = (char*)malloc(READ_NAME_INIT_LENGTH);
		tmp_sub_block->read1[j].seq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read1[j].rseq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read1[j].qual = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read1[j].name_size = READ_NAME_INIT_LENGTH;
		tmp_sub_block->read1[j].seq_size = READ_SEQ_INIT_LENGTH;


		tmp_sub_block->read2[j].name = (char*)malloc(READ_NAME_INIT_LENGTH);
		tmp_sub_block->read2[j].seq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read2[j].rseq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read2[j].qual = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read2[j].name_size = READ_NAME_INIT_LENGTH;
		tmp_sub_block->read2[j].seq_size = READ_SEQ_INIT_LENGTH;
	}




	tmp_sub_block->sub_block_read_number = 0;
}



void init_single_sub_block_single(Read_buffer_single_sub_block* tmp_sub_block)
{

	///buffer_single.sub_block_inner_size = sub_block_inner_size;

	tmp_sub_block->read = (Read*)malloc(sizeof(Read)*buffer_single.sub_block_inner_size);




	for (int j = 0; j < buffer_single.sub_block_inner_size; j++)
	{
		tmp_sub_block->read[j].name = (char*)malloc(READ_NAME_INIT_LENGTH);
		tmp_sub_block->read[j].seq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read[j].rseq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read[j].qual = (char*)malloc(READ_SEQ_INIT_LENGTH);
		tmp_sub_block->read[j].name_size = READ_NAME_INIT_LENGTH;
		tmp_sub_block->read[j].seq_size = READ_SEQ_INIT_LENGTH;
	}




	tmp_sub_block->sub_block_read_number = 0;
}






void init_single_read(Read* read)
{
		read->name = (char*)malloc(READ_NAME_INIT_LENGTH);
		read->seq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		read->rseq = (char*)malloc(READ_SEQ_INIT_LENGTH);
		read->qual = (char*)malloc(READ_SEQ_INIT_LENGTH);
		read->name_size = READ_NAME_INIT_LENGTH;
		read->seq_size = READ_SEQ_INIT_LENGTH;
}


int get_pe_reads_mul_thread(Read_buffer_pe_sub_block* curr_sub_block)
{

	///if (file_format == FASTQ)
	{
		pthread_mutex_lock(&i_readinputMutex);


		///如果缓冲区为0且read还没从文件里读完
		///则消费者要等待，并且通知生成者进行生产
		while (if_empty_pe() && buffer_pe.all_read_end == 0)
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&i_readinputflushCond);



			pthread_cond_wait(&i_readinputstallCond, &i_readinputMutex);
		}


		///如果走到这一步的时候，缓冲区里还有东西
		///则明显要取东西了
		if (!if_empty_pe())
		{
			pop_back_sub_block_pe(curr_sub_block);

			pthread_cond_signal(&i_readinputflushCond);
			pthread_mutex_unlock(&i_readinputMutex);



			///这个时候已经拿到了read，并且出了临界区
			///在这里要对read做后处理
			post_process_paired_reads(curr_sub_block);


			return 1;
		}
		else
		{
			curr_sub_block->sub_block_read_number = 0;

			pthread_cond_signal(&i_readinputstallCond);   //这行非常重要，必须要加！否则会陷入死锁！！！

			pthread_mutex_unlock(&i_readinputMutex);
			return 0;
		}

	}/**
	else
	{
		int file_flag = 0;

		pthread_mutex_lock(&i_readinputMutex);

		load_sub_block_paired_read(curr_sub_block->read1, curr_sub_block->read2,
			buffer_pe.sub_block_inner_size, &file_flag, &curr_sub_block->sub_block_read_number);

		pthread_mutex_unlock(&i_readinputMutex);

		if (file_flag == 0)
		{
			return 0;
		}

		///这个时候已经拿到了read，并且出了临界区
		///在这里要对read做后处理
		post_process_paired_reads(curr_sub_block);

		return 1;
	}
	**/



}





int get_single_reads_mul_thread(Read_buffer_single_sub_block* curr_sub_block)
{

	///if (file_format == FASTQ)
	{
		pthread_mutex_lock(&i_readinputMutex);


		///如果缓冲区为0且read还没从文件里读完
		///则消费者要等待，并且通知生成者进行生产
		while (if_empty_single() && buffer_single.all_read_end == 0)
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&i_readinputflushCond);



			pthread_cond_wait(&i_readinputstallCond, &i_readinputMutex);
		}


		///如果走到这一步的时候，缓冲区里还有东西
		///则明显要取东西了
		if (!if_empty_single())
		{
			pop_back_sub_block_single(curr_sub_block);

			pthread_cond_signal(&i_readinputflushCond);
			pthread_mutex_unlock(&i_readinputMutex);



			///这个时候已经拿到了read，并且出了临界区
			///在这里要对read做后处理
			post_process_single_reads(curr_sub_block);

			return 1;
		}
		else
		{
			curr_sub_block->sub_block_read_number = 0;

			pthread_cond_signal(&i_readinputstallCond);   //这行非常重要，必须要加！否则会陷入死锁！！！

			pthread_mutex_unlock(&i_readinputMutex);

			return 0;
		}
	}/**
	else
	{
		int file_flag = 0;

		pthread_mutex_lock(&i_readinputMutex);

		load_sub_block_single_read(curr_sub_block->read,
			sub_block_inner_size, &file_flag, &curr_sub_block->sub_block_read_number);
		                  
		pthread_mutex_unlock(&i_readinputMutex);

		if (file_flag == 0)
		{
			return 0;
		}

		///这个时候已经拿到了read，并且出了临界区
		///在这里要对read做后处理
		post_process_single_reads(curr_sub_block);

		return 1;
	}
	**/

}








int get_single_reads_mul_thread_pbat(Read_buffer_single_sub_block* curr_sub_block)
{

	///if (file_format == FASTQ)
	{


		pthread_mutex_lock(&i_readinputMutex);


		///如果缓冲区为0且read还没从文件里读完
		///则消费者要等待，并且通知生成者进行生产
		while (if_empty_single() && buffer_single.all_read_end == 0)
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&i_readinputflushCond);



			pthread_cond_wait(&i_readinputstallCond, &i_readinputMutex);
		}


		///如果走到这一步的时候，缓冲区里还有东西
		///则明显要取东西了
		if (!if_empty_single())
		{
			pop_back_sub_block_single(curr_sub_block);

			pthread_cond_signal(&i_readinputflushCond);
			pthread_mutex_unlock(&i_readinputMutex);



			///这个时候已经拿到了read，并且出了临界区
			///在这里要对read做后处理
			post_process_single_reads_pbat(curr_sub_block);









			return 1;
		}
		else
		{
			curr_sub_block->sub_block_read_number = 0;

			pthread_cond_signal(&i_readinputstallCond);   //这行非常重要，必须要加！否则会陷入死锁！！！

			pthread_mutex_unlock(&i_readinputMutex);
			return 0;
		}


	}/**
	else
	{


		int file_flag = 0;

		pthread_mutex_lock(&i_readinputMutex);

		load_sub_block_single_read_pbat(curr_sub_block->read,
			sub_block_inner_size, &file_flag, &curr_sub_block->sub_block_read_number);

		pthread_mutex_unlock(&i_readinputMutex);

		if (file_flag == 0)
		{
			return 0;
		}

		///这个时候已经拿到了read，并且出了临界区
		///在这里要对read做后处理
		post_process_single_reads_pbat(curr_sub_block);

		return 1;
	}
	**/

}



