/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Auxiliary.h"
#include "Process_sam_out.h"
#include <pthread.h>

pthread_mutex_t o_readinputMutex;
pthread_mutex_t o_queueMutex;
pthread_mutex_t o_terminateMutex;
pthread_cond_t o_flushCond;
pthread_cond_t o_readinputflushCond;
pthread_cond_t o_stallCond;
pthread_cond_t o_readinputstallCond;
pthread_mutex_t o_doneMutex;



FILE *_out_fp;


FILE *_out_fp_methy1;
FILE *_out_fp_methy2;
FILE *_out_fp_methy3;




int sub_block_init_length = 10000;


Output_buffer buffer_out;
Output_buffer_sub_block tmp_buffer_sub_block;

Output_methy_buffer buffer_methy_out;
Methylation tmp_methy_sub_block;

Pair_Methylation tmp_methy_sub_block_PE;

Output_methy_buffer_pair buffer_methy_out_PE;


void output_single_methy_CpG
(bitmapper_bs_iter tmp_pos,
char* chrome_name,
int nmethyl,
int total)
{
	int nunmethyl = total - nmethyl;

	double precent;

	if (total == 0)
	{
		precent = 0;
	}
	else
	{
		precent = (double)nmethyl / (double)total;
	}


	fprintf(_out_fp_methy1, "%s\t%llu\t%llu\t%d\t%d\t%d\n",
		chrome_name, tmp_pos, tmp_pos + 1, (int)(precent * 100), nmethyl, nunmethyl);
}



void output_single_methy_CHG
(bitmapper_bs_iter tmp_pos,
char* chrome_name,
int nmethyl,
int total)
{
	int nunmethyl = total - nmethyl;

	double precent;

	if (total == 0)
	{
		precent = 0;
	}
	else
	{
		precent = (double)nmethyl / (double)total;
	}


	fprintf(_out_fp_methy2, "%s\t%llu\t%llu\t%d\t%d\t%d\n",
		chrome_name, tmp_pos, tmp_pos + 1, (int)(precent * 100), nmethyl, nunmethyl);
}


void output_single_methy_CHH
(bitmapper_bs_iter tmp_pos,
char* chrome_name,
int nmethyl,
int total)
{
	int nunmethyl = total - nmethyl;

	double precent;

	if (total == 0)
	{
		precent = 0;
	}
	else
	{
		precent = (double)nmethyl / (double)total;
	}


	fprintf(_out_fp_methy3, "%s\t%llu\t%llu\t%d\t%d\t%d\n",
		chrome_name, tmp_pos, tmp_pos + 1, (int)(precent * 100), nmethyl, nunmethyl);
}




int init_output_methy(char *fileName)
{
	char output_name[NAME_LENGTH];



	sprintf(output_name, "%s_CpG.bedGraph", fileName);

	fprintf(stderr, "%s\n", output_name);

	_out_fp_methy1 = fopen(output_name, "w");

	if (_out_fp_methy1 == NULL)
	{
		return 0;
	}





	sprintf(output_name, "%s_CHG.bedGraph", fileName);

	fprintf(stderr, "%s\n", output_name);

	_out_fp_methy2 = fopen(output_name, "w");

	if (_out_fp_methy2 == NULL)
	{
		return 0;
	}




	sprintf(output_name, "%s_CHH.bedGraph", fileName);

	fprintf(stderr, "%s\n", output_name);

	_out_fp_methy3 = fopen(output_name, "w");

	if (_out_fp_methy3 == NULL)
	{
		return 0;
	}





	return 1;
}






///void finalizeTXOutput()
void finalizeOutput()
{
  fclose(_out_fp);
}

void init_buffer_sub_block(Output_buffer_sub_block* sub_block)
{
	sub_block->length = 0;
	sub_block->size = sub_block_init_length;
	sub_block->buffer = (char*)malloc(sub_block->size);
}




///初始化缓冲区
void init_output_buffer(int thread_number)
{


	///ouput_buffer_size
	///size是空间大小
	buffer_out.sub_block_size = ouput_buffer_size * thread_number;
	buffer_out.sub_block_number = 0;

	buffer_out.sub_buffer = (Output_buffer_sub_block*)malloc(sizeof(Output_buffer_sub_block)*buffer_out.sub_block_size);



	for (int i = 0; i < buffer_out.sub_block_size; i++)
	{
		init_buffer_sub_block(&(buffer_out.sub_buffer[i]));
	}

	buffer_out.all_buffer_end = 0;
}


void init_output_methy_buffer(int thread_number)
{


	///ouput_buffer_size
	///size是空间大小
	buffer_methy_out.sub_block_size = methylation_buffer_times * thread_number;
	
	buffer_methy_out.sub_block_number = 0;

	buffer_methy_out.sub_buffer =
		(Methylation*)malloc(sizeof(Methylation)*buffer_methy_out.sub_block_size);


	
	for (int i = 0; i < buffer_methy_out.sub_block_size; i++)
	{
		init_methylation(&(buffer_methy_out.sub_buffer[i]), methylation_size);
	}

	buffer_methy_out.all_buffer_end = 0;
	
}



void init_output_methy_buffer_pair(int thread_number)
{


	///ouput_buffer_size
	///size是空间大小
	buffer_methy_out_PE.sub_block_size = methylation_buffer_times * thread_number;

	buffer_methy_out_PE.sub_block_number = 0;

	buffer_methy_out_PE.sub_buffer =
		(Pair_Methylation*)malloc(sizeof(Pair_Methylation)*buffer_methy_out_PE.sub_block_size);



	for (int i = 0; i < buffer_methy_out_PE.sub_block_size; i++)
	{
		///注意这里要除以2
		init_pair_methylation(&(buffer_methy_out_PE.sub_buffer[i]), methylation_size/2);
	}

	buffer_methy_out_PE.all_buffer_end = 0;

}


inline int if_empty_buffer()
{

	if (buffer_out.sub_block_number == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_empty_methy_buffer()
{

	if (buffer_methy_out.sub_block_number == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}



inline int if_empty_methy_buffer_pair()
{

	if (buffer_methy_out_PE.sub_block_number == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_full_buffer()
{

	if (buffer_out.sub_block_number >= buffer_out.sub_block_size)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_full_methy_buffer()
{

	if (buffer_methy_out.sub_block_number >= buffer_methy_out.sub_block_size)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


inline int if_full_methy_buffer_pair()
{

	if (buffer_methy_out_PE.sub_block_number >= buffer_methy_out_PE.sub_block_size)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


///这个是要取一个结果取到curr_sub_block中
inline void pop_single_buffer(Output_buffer_sub_block* curr_sub_block)
{
	buffer_out.sub_block_number--;
	///只交换指针
	///记得要交换，直接赋值会内存泄露，造成野指针
	char *k;
	k = buffer_out.sub_buffer[buffer_out.sub_block_number].buffer;
	buffer_out.sub_buffer[buffer_out.sub_block_number].buffer = curr_sub_block->buffer;
	curr_sub_block->buffer = k;

	
	long long tmp_size;
	tmp_size = curr_sub_block->size;
	curr_sub_block->size = buffer_out.sub_buffer[buffer_out.sub_block_number].size;
	buffer_out.sub_buffer[buffer_out.sub_block_number].size = tmp_size;


	curr_sub_block->length = buffer_out.sub_buffer[buffer_out.sub_block_number].length;
	buffer_out.sub_buffer[buffer_out.sub_block_number].length = 0;




}








///buffer_methy_out


inline void pop_single_methy_buffer(Methylation* curr_sub_block)
{
	buffer_methy_out.sub_block_number--;

	


	

	///这个空间应该是genome_cut+1
	bitmapper_bs_iter* cut_index;   ///这个要交换
	///这四个数组构成一个完整的元素
	bitmapper_bs_iter* sites;  ///这个要交换
	uint16_t* r_length;  ///这个要交换
	uint16_t* r_size;  ///这个要交换
	char** reads;  ///这个要交换


	cut_index = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].cut_index;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].cut_index = curr_sub_block->cut_index;
	curr_sub_block->cut_index = cut_index;

	sites = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].sites;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].sites = curr_sub_block->sites;
	curr_sub_block->sites = sites;

	r_length = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_length;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_length = curr_sub_block->r_length;
	curr_sub_block->r_length = r_length;

	r_size = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_size;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_size = curr_sub_block->r_size;
	curr_sub_block->r_size = r_size;


	reads = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].reads;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].reads = curr_sub_block->reads;
	curr_sub_block->reads = reads;

	
	curr_sub_block->current_size = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].current_size;
	///这个clear应该是可以不要的
	clear_methylation(&(buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number]));




}











inline void pop_single_methy_buffer_pair(Pair_Methylation* curr_sub_block)
{
	buffer_methy_out_PE.sub_block_number--;






	///这个空间应该是genome_cut+1
	bitmapper_bs_iter* cut_index;   ///这个要交换
	Single_Methylation tmp_R;///这个要交换


	cut_index = buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].cut_index;
	buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].cut_index = curr_sub_block->cut_index;
	curr_sub_block->cut_index = cut_index;

	///第一个read
	tmp_R = buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[0];
	buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[0] = curr_sub_block->R[0];
	curr_sub_block->R[0] = tmp_R;


	///第二个read
	tmp_R = buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[1];
	buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[1] = curr_sub_block->R[1];
	curr_sub_block->R[1] = tmp_R;


	curr_sub_block->current_size = buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].current_size;
	///这个clear应该是可以不要的
	clear_methylation_pair(&(buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number]));




}







///这个是要把curr_sub_block的结果放到队列中
inline void push_single_buffer(Output_buffer_sub_block* curr_sub_block)
{


	///只交换指针
	///记得要交换，直接赋值会内存泄露，造成野指针
	char *k;
	k = buffer_out.sub_buffer[buffer_out.sub_block_number].buffer;
	buffer_out.sub_buffer[buffer_out.sub_block_number].buffer = curr_sub_block->buffer;
	curr_sub_block->buffer = k;


	long long tmp_size;
	tmp_size = curr_sub_block->size;
	curr_sub_block->size = buffer_out.sub_buffer[buffer_out.sub_block_number].size;
	buffer_out.sub_buffer[buffer_out.sub_block_number].size = tmp_size;




	buffer_out.sub_buffer[buffer_out.sub_block_number].length = curr_sub_block->length;
	curr_sub_block->length = 0;



	buffer_out.sub_block_number++;
}




///这个是要把curr_sub_block的结果放到队列中
inline void push_single_methy_buffer(Methylation* curr_sub_block)
{



	///这个空间应该是genome_cut+1
	bitmapper_bs_iter* cut_index;   ///这个要交换
	///这四个数组构成一个完整的元素
	bitmapper_bs_iter* sites;  ///这个要交换
	uint16_t* r_length;  ///这个要交换
	uint16_t* r_size;  ///这个要交换
	char** reads;  ///这个要交换


	cut_index = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].cut_index;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].cut_index = curr_sub_block->cut_index;
	curr_sub_block->cut_index = cut_index;

	sites = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].sites;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].sites = curr_sub_block->sites;
	curr_sub_block->sites = sites;

	r_length = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_length;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_length = curr_sub_block->r_length;
	curr_sub_block->r_length = r_length;

	r_size = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_size;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].r_size = curr_sub_block->r_size;
	curr_sub_block->r_size = r_size;


	reads = buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].reads;
	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].reads = curr_sub_block->reads;
	curr_sub_block->reads = reads;





	buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number].current_size = curr_sub_block->current_size;
	///这个应该放在临界区之外做
	///clear_methylation(curr_sub_block);



	buffer_methy_out.sub_block_number++;
}




///这个是要把curr_sub_block的结果放到队列中
inline void push_single_methy_buffer_pair(Pair_Methylation* curr_sub_block)
{



	///这个空间应该是genome_cut+1
	bitmapper_bs_iter* cut_index;   ///这个要交换
	Single_Methylation tmp_R;  ///这个要交换



	cut_index = buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].cut_index;
	buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].cut_index = curr_sub_block->cut_index;
	curr_sub_block->cut_index = cut_index;


	///read1
	tmp_R = buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[0];
	buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[0] = curr_sub_block->R[0];
	curr_sub_block->R[0] = tmp_R;


	///read2
	tmp_R = buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[1];
	buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].R[1] = curr_sub_block->R[1];
	curr_sub_block->R[1] = tmp_R;



	buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number].current_size = curr_sub_block->current_size;
	///这个应该放在临界区之外做
	///clear_methylation(curr_sub_block);



	buffer_methy_out_PE.sub_block_number++;
}


void finish_output_buffer()
{

	

		if (output_methy == 0)
		{
			///这个也要加互斥锁
			pthread_mutex_lock(&o_doneMutex);
			buffer_out.all_buffer_end++;


			if (buffer_out.all_buffer_end == THREAD_COUNT)
				pthread_cond_signal(&o_flushCond);
			pthread_mutex_unlock(&o_doneMutex);
		}
		else
		{



			///这个也要加互斥锁
			pthread_mutex_lock(&o_doneMutex);
			buffer_methy_out.all_buffer_end++;


			if (buffer_methy_out.all_buffer_end == THREAD_COUNT)
				pthread_cond_signal(&o_flushCond);
			pthread_mutex_unlock(&o_doneMutex);
		}
	

	
}



void finish_output_buffer_pair()
{
	if (output_methy == 0)
	{
		///这个也要加互斥锁
		pthread_mutex_lock(&o_doneMutex);
		buffer_out.all_buffer_end++;


		if (buffer_out.all_buffer_end == THREAD_COUNT)
			pthread_cond_signal(&o_flushCond);
		pthread_mutex_unlock(&o_doneMutex);
	}
	else
	{



		///这个也要加互斥锁
		pthread_mutex_lock(&o_doneMutex);
		buffer_methy_out_PE.all_buffer_end++;


		if (buffer_methy_out_PE.all_buffer_end == THREAD_COUNT)
			pthread_cond_signal(&o_flushCond);
		pthread_mutex_unlock(&o_doneMutex);
	}


}



void push_results_to_buffer(Output_buffer_sub_block* sub_block)
{
	
	pthread_mutex_lock(&o_queueMutex);

	///if_full_buffer
	while (if_full_buffer())
	{
		///按道理这个信号量似乎不用发
		///因为队列不可能一边满一边空
		pthread_cond_signal(&o_flushCond);


		pthread_cond_wait(&o_stallCond, &o_queueMutex);
	}


	push_single_buffer(sub_block);



	pthread_cond_signal(&o_flushCond);
	pthread_mutex_unlock(&o_queueMutex);
	
}


void push_methy_to_buffer(Methylation* methy)
{

	pthread_mutex_lock(&o_queueMutex);

	///if_full_buffer
	while (if_full_methy_buffer())
	{
		///按道理这个信号量似乎不用发
		///因为队列不可能一边满一边空
		pthread_cond_signal(&o_flushCond);


		pthread_cond_wait(&o_stallCond, &o_queueMutex);
	}


	push_single_methy_buffer(methy);



	pthread_cond_signal(&o_flushCond);
	pthread_mutex_unlock(&o_queueMutex);

	clear_methylation(methy);

}



void push_methy_to_buffer_pair(Pair_Methylation* methy)
{

	pthread_mutex_lock(&o_queueMutex);

	///if_full_buffer
	while (if_full_methy_buffer_pair())
	{
		///按道理这个信号量似乎不用发
		///因为队列不可能一边满一边空
		pthread_cond_signal(&o_flushCond);


		pthread_cond_wait(&o_stallCond, &o_queueMutex);
	}


	push_single_methy_buffer_pair(methy);



	pthread_cond_signal(&o_flushCond);
	pthread_mutex_unlock(&o_queueMutex);

	clear_methylation_pair(methy);

}


void* pop_buffer(void*)
{
	
	FILE* output_file = get_Ouput_Dec();

	init_buffer_sub_block(&tmp_buffer_sub_block);

	

	while (buffer_out.all_buffer_end<THREAD_COUNT)
	{

		pthread_mutex_lock(&o_queueMutex);


		while (if_empty_buffer() && (buffer_out.all_buffer_end<THREAD_COUNT))
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&o_stallCond);


			pthread_cond_wait(&o_flushCond, &o_queueMutex);
		}


		if (!if_empty_buffer())
		{
			pop_single_buffer(&tmp_buffer_sub_block);
		}

		pthread_cond_signal(&o_stallCond);
		pthread_mutex_unlock(&o_queueMutex);

		if (tmp_buffer_sub_block.length != 0)
		{
			fprintf(output_file, "%s", tmp_buffer_sub_block.buffer);
		}

	}


	while (buffer_out.sub_block_number>0)
	{
		buffer_out.sub_block_number--;
		fprintf(output_file, "%s", buffer_out.sub_buffer[buffer_out.sub_block_number].buffer);

	}

}




void* pop_methy_buffer_pair(void*)
{
	///要求这个大小和每个线程内部的block大小一样
	int my_methylation_size = methylation_size;
	///注意这里要除以2
	init_pair_methylation(&tmp_methy_sub_block_PE, my_methylation_size/2);


	while (buffer_methy_out_PE.all_buffer_end<THREAD_COUNT)
	{

		pthread_mutex_lock(&o_queueMutex);


		while (if_empty_methy_buffer_pair() && (buffer_methy_out_PE.all_buffer_end<THREAD_COUNT))
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&o_stallCond);


			pthread_cond_wait(&o_flushCond, &o_queueMutex);
		}


		if (!if_empty_methy_buffer_pair())
		{
			pop_single_methy_buffer_pair(&tmp_methy_sub_block_PE);
		}

		pthread_cond_signal(&o_stallCond);
		pthread_mutex_unlock(&o_queueMutex);

		if (tmp_methy_sub_block_PE.current_size != 0)
		{
			///这个纯输出，不用划分，划分在计算线程里做
			///assign_cuts(&tmp_methy_sub_block);

			output_methylation_pair(&tmp_methy_sub_block_PE);

			///按道理也是不用清理的
			///clear_methylation(&tmp_methy_sub_block);
		}

	}


	while (buffer_methy_out_PE.sub_block_number>0)
	{
		buffer_methy_out_PE.sub_block_number--;

		///这个纯输出，不用划分，划分在计算线程里做
		///assign_cuts(&tmp_methy_sub_block);

		output_methylation_pair(&(buffer_methy_out_PE.sub_buffer[buffer_methy_out_PE.sub_block_number]));

		///按道理也是不用清理的
		///clear_methylation(&tmp_methy_sub_block);

	}

}



void* pop_methy_buffer(void*)
{
	///要求这个大小和每个线程内部的block大小一样
	int my_methylation_size = methylation_size;

	init_methylation(&tmp_methy_sub_block, my_methylation_size);


	while (buffer_methy_out.all_buffer_end<THREAD_COUNT)
	{

		pthread_mutex_lock(&o_queueMutex);


		while (if_empty_methy_buffer() && (buffer_methy_out.all_buffer_end<THREAD_COUNT))
		{
			///按道理这个信号量似乎不用发
			///因为队列不可能一边满一边空
			pthread_cond_signal(&o_stallCond);


			pthread_cond_wait(&o_flushCond, &o_queueMutex);
		}


		if (!if_empty_methy_buffer())
		{
			pop_single_methy_buffer(&tmp_methy_sub_block);
		}

		pthread_cond_signal(&o_stallCond);
		pthread_mutex_unlock(&o_queueMutex);

		if (tmp_methy_sub_block.current_size != 0)
		{
			///这个纯输出，不用划分，划分在计算线程里做
			///assign_cuts(&tmp_methy_sub_block);

			output_methylation(&tmp_methy_sub_block);

			///按道理也是不用清理的
			///clear_methylation(&tmp_methy_sub_block);
		}

	}


	while (buffer_methy_out.sub_block_number>0)
	{
		buffer_methy_out.sub_block_number--;
		
		///这个纯输出，不用划分，划分在计算线程里做
		///assign_cuts(&tmp_methy_sub_block);

		output_methylation(&(buffer_methy_out.sub_buffer[buffer_methy_out.sub_block_number]));

		///按道理也是不用清理的
		///clear_methylation(&tmp_methy_sub_block);

	}

}


void OutPutSAM_Nounheader(_rg_name_l  *_ih_refGenName, int refChromeCont, int argc, char *argv[])
{
	fprintf(_out_fp, "@HD\tVN:1.4\tSO:unsorted\n");
	int i = 0;
	for (i = 0; i<refChromeCont; i++)
	{
		fprintf(_out_fp, "@SQ\tSN:%s\tLN:%u\n", _ih_refGenName[i]._rg_chrome_name, _ih_refGenName[i]._rg_chrome_length);
	}

	fprintf(_out_fp, "@PG\tID:BitMapperBS\tVN:%s\tCL:", versionN);
	for (size_t i = 0; i < argc; i++)
	{
		fprintf(_out_fp, "%s ", argv[i]);
	}

	fprintf(_out_fp, "\n");
}

FILE* get_Ouput_Dec()
{
    return _out_fp;
}


void multi_outputQ(SAM map, int thread_id)
{
    if(map.FLAG<=16)
    {
        fprintf(_out_fp, "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
	  map.QNAME,
	  map.FLAG,
	  map.RNAME,
	  map.POS,
	  map.MAPQ,
	  map.CIGAR,
	  map.MRNAME,
	  map.MPOS,
	  map.ISIZE,
	  map.SEQ,
	  map.QUAL);
    }
    else
    {
        fprintf(_out_fp, "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%c\t%c",
	  map.QNAME,
	  map.FLAG,
	  map.RNAME,
	  map.POS,
	  map.MAPQ,
	  map.CIGAR,
	  map.MRNAME,
	  map.MPOS,
	  map.ISIZE,
	  '*',
	  '*');
    }
  int i;

  for ( i = 0; i < map.optSize; i++)
    {
      switch (map.optFields[i].type)
	{
	case 'A':
	  fprintf(_out_fp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
	  break;
	case 'i':
	  fprintf(_out_fp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
	  break;
	case 'f':
	  fprintf(_out_fp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
	  break;
	case 'Z':
	case 'H':
	  fprintf(_out_fp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
	  break;
	}
    }

  fprintf(_out_fp, "\n");
}






int Output_gene ( char *fileName)
{


    _out_fp = fopen (fileName, "w");
    char tmpFileName[NAME_LENGTH];

    if (_out_fp == NULL)
	{
	  return 0;
	}
    ///finalizeOutput = &finalizeTXOutput;

  return 1;
}

