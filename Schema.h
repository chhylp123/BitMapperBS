/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/


#ifndef __MR_FAST__
#define __MR_FAST__

#include "Process_Reads.h"
#include <ctype.h>
///#include "Levenshtein_Cal.h"
#include <sys/queue.h>
#include <pthread.h>
#include "bwt.h"


#define MAP_CHUNKS 15
#define MAX_CIGAR_SIZE 1000
#define INIT_READ_LENGTH 155



#define Check_Space(read, c_length, n_length) if(n_length > c_length) {read = (char*)realloc(read, n_length);c_length = n_length;}

#define Get_Count(arr, i) (arr[i+1]-arr[i])








typedef struct _rg_name_l
{
	char _rg_chrome_name[SEQ_MAX_LENGTH];
	bitmapper_bs_iter _rg_chrome_length;
	bitmapper_bs_iter start_location;
	bitmapper_bs_iter end_location;
}_rg_name_l;


typedef struct seed_votes
{
	bitmapper_bs_iter site;
	bitmapper_bs_iter vote;
	unsigned int err;
	///int score;
	bitmapper_bs_iter end_site;
}seed_votes;

typedef struct
{
    unsigned int hash_key;
    unsigned int site;
    unsigned int *locs;
    unsigned int Cnt;
}Candidate;

typedef struct
{
    unsigned int start_offset;
    unsigned int start_locs;
    int err;
}Result_Read1;



typedef struct
{
    unsigned short FLAG;
    int readNumber;
    unsigned short direction;
    unsigned int map_location;
    int chmore_number;
    char cigar[MAX_CIGAR_SIZE];
    unsigned short iVal;
    char sVal[MAX_CIGAR_SIZE];
} tmp_result;

typedef struct
{
	int flag;
	bitmapper_bs_iter origin_site;
	bitmapper_bs_iter site;
	bitmapper_bs_iter end_site;
	int chrome_id;
	char cigar[SEQ_MAX_LENGTH];
	unsigned int err;
} map_result;


typedef struct
{
	///这个空间应该是genome_cut+1
	bitmapper_bs_iter* cut_index;
	bitmapper_bs_iter genome_cut;
	bitmapper_bs_iter cut_length;
	bitmapper_bs_iter* tmp_index; ///这个就是移位用的

	int total_size;
	int current_size;

	///这四个数组构成一个完整的元素
	bitmapper_bs_iter* sites;
	uint16_t* r_length;
	uint16_t* r_size;
	char** reads;

	///下面这四个元素纯粹是用来做swap的暂存变量
	bitmapper_bs_iter k_sites;
	uint16_t k_r_length;
	uint16_t k_r_size;
	char* k_reads;


}Methylation;


typedef struct
{
	FILE** files;
	char filename[SEQ_MAX_LENGTH];
	bitmapper_bs_iter* cut_index;
	bitmapper_bs_iter genome_cut;
}Methylation_Output;



struct result_queue{
    TAILQ_ENTRY(result_queue) result_tailq;
    tmp_result single_result;
};

TAILQ_HEAD(Result_Queue,result_queue);

//struct Result_Queue queue_header;  ///queue_header是队列的头部

extern long long			mappingCnt[MAX_Thread];
extern long long			mappedSeqCnt[MAX_Thread];
extern long long			completedSeqCnt[MAX_Thread];

extern Methylation_Output methy_out;


void Prepare_alignment(char* outputFileName, char *genFileName, _rg_name_l *chhy_ih_refGenName, int chhy_refChromeCont, int i_read_format);
int Map_Single_Seq( int thread_id);
int Map_Single_Seq_pbat(int thread_id);
///void* Map_Single_Seq_split(int thread_id);
void* Map_Single_Seq_split(void* arg);
void* Map_Single_Seq_split_pbat(void* arg);
void* Map_Single_Seq_split_one_by_one(void* arg);
int Map_Single_Seq_muti_thread(int thread_id);
int Map_Single_Seq_pbat_muti_thread(int thread_id);


int Map_Pair_Seq(int thread_id);
void* Map_Pair_Seq_split(void* arg);
int Map_Pair_Seq_muti_thread(int thread_id);
void Prepare_methy(char *genFileName, _rg_name_l *chhy_ih_refGenName, int chhy_refChromeCont);


void get_mapping_informations(long long* number_of_read, long long* number_of_unique_mapped_read,
	long long* number_of_ambiguous_mapped_read, long long* number_of_unmapped_read);



///下面是我到函数



unsigned int get_total_reads_number();

void open_methy_output(Methylation_Output& methy_out);


inline void correct_read_using_cigar(char* cigar, char* input_read, char* output_read)
{
	///这下面是遍历cigar用的
	int i = 0;
	int cur_i = 0;

	///这下面是遍历read用的
	///int r_length = 0;
	int input_length = 0;
	int gap_size = 0;
	int convert_length;
	int output_length = 0;

	while (cigar[i] != '\0')
	{

		switch (cigar[i]){
		case 'M':
			cigar[i] = '\0';
			convert_length = atoi(cigar + cur_i);
			memcpy(output_read + output_length, input_read + input_length, convert_length);
			input_length = input_length + convert_length;
			output_length = output_length + convert_length;
			cur_i = i + 1;
			break; 
		case 'I':
			cigar[i] = '\0';
			convert_length = atoi(cigar + cur_i);
			input_length = input_length + convert_length;
			gap_size = gap_size - convert_length;
			cur_i = i + 1;
			break;
		case 'D':
			cigar[i] = '\0';
			convert_length = atoi(cigar + cur_i);
			memset(output_read + output_length, '_', convert_length);
			output_length = output_length + convert_length;
			gap_size = gap_size + convert_length;
			cur_i = i + 1;
			break;
		}

		i++;
	}

	output_read[output_length] = '\0';

	
}

inline void init_methylation(Methylation* methy, int size)
{
	(*methy).total_size = size;
	(*methy).current_size = 0;
	(*methy).r_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).r_size = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).sites = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*size);
	(*methy).reads = (char**)malloc(sizeof(char*)*size);
	(*methy).genome_cut = genome_cuts;
	(*methy).cut_length = cut_length;
	(*methy).cut_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));
	(*methy).tmp_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));
	int i;
	for (i = 0; i < size; i++)
	{
		(*methy).reads[i] = (char*)malloc(sizeof(char)*INIT_READ_LENGTH);
		(*methy).r_size[i] = INIT_READ_LENGTH;
		(*methy).r_length[i] = 0;
	}

	/**
	for (i = 0; i <= (*methy).genome_cut; i++)
	{
		(*methy).cut_index[i] = 0;
	}
	**/

	memset((*methy).cut_index, 0, sizeof(bitmapper_bs_iter)*((*methy).genome_cut + 1));
}

/**
inline void output_methylation(Methylation* methy, FILE *schema_out_fp)
{
	int i = 0;

	fprintf(schema_out_fp, "#%llu\t%llu\t%llu", 
		(*methy).current_size, (*methy).genome_cut, (*methy).cut_length);

	for (i = 0; i <= (*methy).genome_cut; i++)
	{

		fprintf(schema_out_fp, "\t%llu",
			(*methy).cut_index[i]);
	}
	fprintf(schema_out_fp, "\n");

	for (i = 0; i < (*methy).current_size; i++)
	{
		fprintf(schema_out_fp, "%llu\t%s\n", (*methy).sites[i], (*methy).reads[i]);
	}

}
**/

inline void output_methylation(Methylation* methy)
{
	bitmapper_bs_iter i;
	for (i = 0; i <= (*methy).genome_cut; i++)
	{
		methy_out.cut_index[i] += Get_Count((*methy).cut_index, i);
	}

	bitmapper_bs_iter index;
	bitmapper_bs_iter end_i;
	FILE* file;

	///外层循环是每个cut
	for (index = 0; index < (*methy).genome_cut; index++)
	{
		i = (*methy).cut_index[index];
		end_i = (*methy).cut_index[index + 1];
		file = methy_out.files[index];

		for (; i < end_i; i++)
		{
			fprintf(file, "%llu\t%llu\t%s\n", (*methy).r_length[i], (*methy).sites[i], (*methy).reads[i]);
		}
	}

}



inline int get_cut_id(Methylation* methy, bitmapper_bs_iter site)
{
	return site / (*methy).cut_length;
}

inline void swap(Methylation* methy, bitmapper_bs_iter i, bitmapper_bs_iter j)
{

	(*methy).k_sites = (*methy).sites[i];
	(*methy).k_r_length = (*methy).r_length[i];
	(*methy).k_r_size = (*methy).r_size[i];
	(*methy).k_reads = (*methy).reads[i];

	(*methy).sites[i] = (*methy).sites[j];
	(*methy).r_length[i] = (*methy).r_length[j];
	(*methy).r_size[i] = (*methy).r_size[j];
	(*methy).reads[i] = (*methy).reads[j];

	(*methy).sites[j] = (*methy).k_sites;
	(*methy).r_length[j] = (*methy).k_r_length;
	(*methy).r_size[j] = (*methy).k_r_size;
	(*methy).reads[j] = (*methy).k_reads;


}

inline void assign_cuts(Methylation* methy)
{
	bitmapper_bs_iter i = 0;


	(*methy).tmp_index[0] = 0;

	for (i = 0; i < (*methy).genome_cut; i++)
	{
		(*methy).tmp_index[i + 1] = (*methy).tmp_index[i] + (*methy).cut_index[i];
	}

	memcpy((*methy).cut_index, (*methy).tmp_index, sizeof(bitmapper_bs_iter)*((*methy).genome_cut + 1));


	bitmapper_bs_iter index;
	bitmapper_bs_iter end_i;
	bitmapper_bs_iter tmp_index;

	///外层循环是循环每个桶的
	for (index = 0; index < (*methy).genome_cut; index++)
	{
		i = (*methy).tmp_index[index];
		end_i = (*methy).cut_index[index + 1];
		///内层是循环每个桶
		while (i< end_i)
		{
			tmp_index = get_cut_id(methy, (*methy).sites[i]);

			if (tmp_index == index)
			{
				i++;
				(*methy).tmp_index[index]++; ///这个可以不加
			}
			else
			{
				swap(methy, i, (*methy).tmp_index[tmp_index]);
				(*methy).tmp_index[tmp_index]++;
			}
		}

	}



	/**
	for (index = 0; index < (*methy).genome_cut; index++)
	{
		i = (*methy).cut_index[index];
		end_i = (*methy).cut_index[index + 1];

		for (; i < end_i; i++)
		{
			tmp_index = get_cut_id(methy, (*methy).sites[i]);


			fprintf(stderr, "%llu(%llu)\t", (*methy).sites[i], tmp_index);

		}
		fprintf(stderr, "\n");
	}



	fprintf(stderr, "***********************\n\n\n");
	**/

	/**
	for (index = 0; index < (*methy).genome_cut; index++)
	{
		i = (*methy).cut_index[index];
		end_i = (*methy).cut_index[index + 1];

		if ((*methy).tmp_index[index] != (*methy).cut_index[index + 1])
		{
			fprintf(stderr, "index: %llu\n", index);
			fprintf(stderr, "tmp_index: %llu\n", (*methy).tmp_index[index]);
			fprintf(stderr, "cut_index: %llu\n", (*methy).cut_index[index + 1]);
		}

		for (; i < end_i; i++)
		{
			tmp_index = get_cut_id(methy, (*methy).sites[i]);

			if (tmp_index != index)
			{
				fprintf(stderr, "error\n");
			}
		}
	}
	
	**/

	
	/**
	for (i = 0; i < (*methy).genome_cut + 1; i++)
	{
		if ((*methy).tmp_index[i] != (*methy).cut_index[i])
		{
			fprintf(stderr, "error\n");
		}
		
	}
	**/
	

	/**
	for (i = 0; i < (*methy).genome_cut; i++)
	{
		if (Get_Count((*methy).tmp_index, i) != (*methy).cut_index[i])
			fprintf(stderr, "error\n");
	}
	**/
	

	/**
	fprintf(stderr, "*************************************\n");

	for (i = 0; i < (*methy).genome_cut; i++)
	{
		fprintf(stderr, "i: %d, %llu\n", i, (*methy).cut_index[i]);
	}

	fprintf(stderr, "\n");

	for (i = 0; i < (*methy).current_size; i++)
	{
		(*methy).cut_index[get_cut_id(methy, (*methy).sites[i])]--;
	}

	for (i = 0; i < (*methy).genome_cut; i++)
	{
		if ((*methy).cut_index[i]!=0)
		{
			fprintf(stderr, "error\n");
		}
		fprintf(stderr, "i: %d, %llu\n", i, (*methy).cut_index[i]);
	}

	fprintf(stderr, "\n");


	

	fprintf(stderr, "*************************************\n");
	**/
	
}



inline void clear_methylation(Methylation* methy)
{
	(*methy).current_size = 0;

	/**
	int i;
	for (i = 0; i < (*methy).genome_cut; i++)
	{
		(*methy).cut_index[i] = 0;
	}
	**/

	memset((*methy).cut_index, 0, sizeof(bitmapper_bs_iter)*((*methy).genome_cut + 1));
}



inline void C_to_T(char *Seq, char *bsSeq, int length, int* C_site)
{
	int i;
	(*C_site) = -1;
	for (i = 0; i<length; i++)
	{
		if (rc_table[Seq[i]] == 'C')
		///if (Seq[i] == 'C')
		///if (Seq[length - 1 - i] == 'C')
		{
			(*C_site) = i;
			bsSeq[i] = 'T';
		}
		else
		{
			bsSeq[i] = rc_table[Seq[i]];
			///bsSeq[i] = Seq[i];
		}
	}

	bsSeq[length] = '\0';
}


inline void C_to_T_forward(char *Seq, char *bsSeq, int length, int* C_site)
{
	int i;
	(*C_site) = -1;
	for (i = 0; i<length; i++)
	{
		bsSeq[i] = Seq[length - 1 - i];

		if (bsSeq[i] == 'C')
		{
			(*C_site) = i;
			bsSeq[i] = 'T';
		}

	}

	bsSeq[length] = '\0';
}


inline void C_to_T_forward_array(char *Seq, char *bsSeq, int length, int* C_site_number, int* C_site_array)
{
	int i;
	(*C_site_number) = 0;
	for (i = 0; i<length; i++)
	{
		bsSeq[i] = Seq[length - 1 - i];

		if (bsSeq[i] == 'C')
		{
			C_site_array[(*C_site_number)] = length - 1 - i;
			///C_site_array[(*C_site_number)] = i;
			(*C_site_number)++;
			bsSeq[i] = 'T';
		}

	}

	if ((*C_site_number) == 0)
	{
		C_site_array[0] = SEQ_MAX_LENGTH;
	}

	bsSeq[length] = '\0';
}


void methy_extract(int thread_id, char* file_name);


#endif
