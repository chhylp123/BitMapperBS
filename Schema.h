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
///这个是字符长度, 现在不用了
#define INIT_READ_LENGTH 155  
///这个是uint64_t的长度, 一个uint64_t可以存21个3-bit, 所以这个字符长度是INIT_READ_3_BIT_LENGTH*21
#define INIT_READ_3_BIT_LENGTH 8  


#define FORWARD_DUP_MASK 128
#define RC_DUP_MASK 64

#define Check_Space(read, c_length, n_length) if(n_length > c_length) {read = (char*)realloc(read, n_length);c_length = n_length;}

///注意c_length_3_bit是uint64_t的长度
///n_length是char的长度
///所以c_length_3_bit的char长度是c_length_3_bit*21
#define Check_Space_3_Bit(read, c_length_3_bit, n_length) if(n_length > c_length_3_bit * 21) \
										{c_length_3_bit = n_length / 21 + (n_length % 21 != 0); \
										read = (bitmapper_bs_iter*)realloc(read, c_length_3_bit * sizeof(bitmapper_bs_iter)); }


#define Get_Count(arr, i) (arr[i+1]-arr[i])

/**
A->0
C->1
G->2
T->3
N-4
**/
#define Get_BS_conversion(C_or_G)    if (C_or_G == 1)/*C methylation*/\
									{\
										match = 1;/* 'C';*/\
										convert = 3;/* 'T';*/\
										rc_match = 2;/* 'G';*/\
										rc_convert = 0;/* 'A';*/\
									}\
									else  /*G methylation*/\
									{\
										match = 2;/* 'G';*/\
										convert = 0;/* 'A';*/\
										rc_match = 1;/* 'C';*/\
										rc_convert = 3;/* 'T';*/\
									}


#define POS_MASK 15
#define CONTEXT_MASK 48

extern char char_to_3bit[128];



typedef struct deduplicate_PE
{
	uint8_t* forward_bits;
	uint8_t* rc_bits;
	long long each_site_bits;
	long long total_sites;
	long long each_site_uint8;
}deduplicate_PE;



inline void init_deduplicate_PE(deduplicate_PE* de, long long each_site_bits, long long total_sites)
{



	de->each_site_bits = each_site_bits;
	de->total_sites = total_sites;

	de->each_site_uint8 = each_site_bits / 8;

	if (each_site_bits % 8 != 0)
	{
		de->each_site_uint8++;
	}

	de->forward_bits = (uint8_t*)malloc(sizeof(uint8_t)*de->each_site_uint8*de->total_sites);
	de->rc_bits = (uint8_t*)malloc(sizeof(uint8_t)*de->each_site_uint8*de->total_sites);

	memset(de->forward_bits, 0, sizeof(uint8_t)*de->each_site_uint8*de->total_sites);
	memset(de->rc_bits, 0, sizeof(uint8_t)*de->each_site_uint8*de->total_sites);

}


inline void re_init_deduplicate_PE(deduplicate_PE* de)
{
	memset(de->forward_bits, 0, sizeof(uint8_t)*de->each_site_uint8*de->total_sites);
	memset(de->rc_bits, 0, sizeof(uint8_t)*de->each_site_uint8*de->total_sites);
}

///把i这个位置置1
inline void set_deduplicate_PE_mask(deduplicate_PE* de, uint8_t* read_bits, bitmapper_bs_iter site, bitmapper_bs_iter i)
{
	bitmapper_bs_iter offset = i & 7;
	uint8_t mask = 1;
	mask = mask << offset;

	offset = i >> 3;

	read_bits[offset + site*de->each_site_uint8] |= mask;
}


///把i这个位置置0
inline void unset_deduplicate_PE_mask(deduplicate_PE* de, uint8_t* read_bits, bitmapper_bs_iter site, bitmapper_bs_iter i)
{
	bitmapper_bs_iter offset = i & 7;
	uint8_t mask = 254;
	mask = mask << offset;

	offset = i >> 3;

	read_bits[offset + site*de->each_site_uint8] &= mask;
}


///把i这个位置的值返回回来
inline uint8_t get_deduplicate_PE_mask(deduplicate_PE* de, uint8_t* read_bits, bitmapper_bs_iter site, bitmapper_bs_iter i)
{
	bitmapper_bs_iter inner_offset = i & 7;
	bitmapper_bs_iter offset = i >> 3;


	uint8_t mask = 1;

	mask = (read_bits[offset + site*de->each_site_uint8] >> inner_offset) & mask;

	return mask;
}






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
	int score;
} map_result;



typedef struct
{
	///bitmapper_bs_iter minDistance_pair;
	bitmapper_bs_iter maxDistance_pair;
	///bitmapper_bs_iter length;
	bitmapper_bs_iter genome_cuts;
	bitmapper_bs_iter** count;
} pair_distance_count;



typedef struct
{
	bitmapper_bs_iter count;
	short index;
} pe_distance_info;


typedef struct
{
	bitmapper_bs_iter maxDistance_pair;
	bitmapper_bs_iter genome_cuts;
	pe_distance_info* count;
	
} pair_distance_result;


inline void init_pe_distance(pair_distance_count* distances, bitmapper_bs_iter minDistance_pair, bitmapper_bs_iter maxDistance_pair)
{
	distances->genome_cuts = genome_cuts;
	///distances->minDistance_pair = minDistance_pair;
	distances->maxDistance_pair = maxDistance_pair;
	///distances->length = distances->maxDistance_pair - distances->minDistance_pair + 1;
	distances->count = (bitmapper_bs_iter**)malloc(sizeof(bitmapper_bs_iter*)*distances->genome_cuts);
	int i;
	for (i = 0; i < distances->genome_cuts; i++)
	{
		///distances->count[i] = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*distances->length);
		///这个一定一定要+1
		distances->count[i] = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(distances->maxDistance_pair + 1));
		///memset(distances->count[i], 0, sizeof(bitmapper_bs_iter)*(distances->length));
		///这个一定一定要+1
		memset(distances->count[i], 0, sizeof(bitmapper_bs_iter)*(distances->maxDistance_pair + 1));
	}

}

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
	///uint16_t* r_size;
	//char** reads;

	uint16_t* r_real_length;
	uint16_t* r_size_3_bit;
	bitmapper_bs_iter** reads_3_bit;



	///下面这四个元素纯粹是用来做swap的暂存变量
	bitmapper_bs_iter k_sites;
	uint16_t k_r_length;
	//uint16_t k_r_size;
	//char* k_reads;

	uint16_t k_r_real_length;
	uint16_t k_r_size_3_bit;
	bitmapper_bs_iter* k_reads_3_bit;


}Methylation;


typedef struct
{
	///这四个数组构成一个完整的元素
	bitmapper_bs_iter* sites;
	uint16_t* r_length;
	///uint16_t* r_size;
	///char** reads;

	uint16_t* r_real_length;
	uint16_t* r_size_3_bit;
	bitmapper_bs_iter** reads_3_bit;


	///下面这四个元素纯粹是用来做swap的暂存变量
	bitmapper_bs_iter k_sites;
	uint16_t k_r_length;
	uint16_t k_r_real_length;
	///uint16_t k_r_size;
	///char* k_reads;
	uint16_t k_r_size_3_bit;
	bitmapper_bs_iter* k_reads_3_bit;

}Single_Methylation;

typedef struct
{
	int total_size;
	int current_size;
	Single_Methylation R[2];


	///这个空间应该是genome_cut+1
	bitmapper_bs_iter* cut_index;
	bitmapper_bs_iter genome_cut;
	bitmapper_bs_iter cut_length;
	bitmapper_bs_iter* tmp_index; ///这个就是移位用的

}Pair_Methylation;






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


void Prepare_alignment(char* outputFileName, char *genFileName, _rg_name_l *chhy_ih_refGenName, int chhy_refChromeCont, int i_read_format, int is_pairedEnd);
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
	long long* number_of_ambiguous_mapped_read, long long* number_of_unmapped_read, long long* number_of_mapped_bases, 
	long long* number_of_mapped_errors);



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



inline void correct_read_using_cigar_3_bit(char* cigar, char* input_read, 
	bitmapper_bs_iter* output_read_3_bit)
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

	int word64_index = 0;
	int word64_inner_bit_index = 0;
	int inner_i = 0;
	bitmapper_bs_iter tmp;
	

	while (cigar[i] != '\0')
	{

		switch (cigar[i]){
		case 'M':
			cigar[i] = '\0';
			convert_length = atoi(cigar + cur_i);

			/******************这里要改***********************/
			///memcpy(output_read + output_length, input_read + input_length, convert_length);
			for (inner_i = 0; inner_i < convert_length; inner_i++)
			{
				if (word64_inner_bit_index == 0)
				{
					output_read_3_bit[word64_index] = 0;
				}

				tmp = char_to_3bit[input_read[input_length + inner_i]];
				tmp = tmp << word64_inner_bit_index * 3;
		
				output_read_3_bit[word64_index] = output_read_3_bit[word64_index] | tmp;

				word64_inner_bit_index++;

				if (word64_inner_bit_index == 21)
				{
					word64_inner_bit_index = 0;
					word64_index++;
				}

				

			}
			/******************这里要改***********************/



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

			/******************这里要改***********************/
			///memset(output_read + output_length, '_', convert_length);
			for (inner_i = 0; inner_i < convert_length; inner_i++)
			{
				if (word64_inner_bit_index == 0)
				{
					output_read_3_bit[word64_index] = 0;
				}

				tmp = char_to_3bit['N'];
				tmp = tmp << word64_inner_bit_index * 3;

				output_read_3_bit[word64_index] = output_read_3_bit[word64_index] | tmp;

				word64_inner_bit_index++;

				if (word64_inner_bit_index == 21)
				{
					word64_inner_bit_index = 0;
					word64_index++;
				}



			}



			/******************这里要改***********************/








			output_length = output_length + convert_length;
			gap_size = gap_size + convert_length;
			cur_i = i + 1;
			break;
		}

		i++;
	}

	///output_read[output_length] = '\0';


}



inline void correct_read_using_cigar_3_bit_over_lap_length_r2(char* cigar, char* input_read,
	bitmapper_bs_iter* output_read_3_bit,
	long long min1, long long len1, long long min2, long long len2)
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

	int word64_index = 0;
	int word64_inner_bit_index = 0;
	int inner_i = 0;
	bitmapper_bs_iter tmp;

	int is_overlap = 0;






	long long min = min1;
	long long max = min1 + len1;

	if (min1 > min2)
		min = min2;

	if (min1 + len1 < min2 + len2)
		max = min2 + len2;

	long long overlap_length = max - min - len1 - len2;

	long long region_start2;
	int region_length2;
	///如果不重叠是正值或者0, 重叠是负值
	if (overlap_length < 0)
	{
		is_overlap = 1;

		overlap_length = overlap_length*-1;

		if (min1 <= min2)
		{
			region_start2 = 0;
			region_length2 = overlap_length;
		}
		else if (min1>min2)
		{
			region_start2 = min1 - min2;
			region_length2 = overlap_length;
		}

		///memset(read2 + region_start2, 'N', region_length2);

	}









	while (cigar[i] != '\0')
	{

		switch (cigar[i]){
		case 'M':
			cigar[i] = '\0';
			convert_length = atoi(cigar + cur_i);

			/******************这里要改***********************/
			///memcpy(output_read + output_length, input_read + input_length, convert_length);
			///input_length = input_length + convert_length;
			///output_length = output_length + convert_length;
			for (inner_i = 0; inner_i < convert_length; inner_i++, output_length++, input_length++)
			{
				if (word64_inner_bit_index == 0)
				{
					output_read_3_bit[word64_index] = 0;
				}

				if (is_overlap && output_length >= region_start2&&output_length<region_start2 + region_length2)
				{
					tmp = char_to_3bit['N'];
				}
				else
				{
					tmp = char_to_3bit[input_read[input_length]];
				}
				tmp = tmp << word64_inner_bit_index * 3;

				output_read_3_bit[word64_index] = output_read_3_bit[word64_index] | tmp;

				word64_inner_bit_index++;

				if (word64_inner_bit_index == 21)
				{
					word64_inner_bit_index = 0;
					word64_index++;
				}



			}
			/******************这里要改***********************/



			
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

			/******************这里要改***********************/
			///memset(output_read + output_length, '_', convert_length);
			for (inner_i = 0; inner_i < convert_length; inner_i++)
			{
				if (word64_inner_bit_index == 0)
				{
					output_read_3_bit[word64_index] = 0;
				}

				tmp = char_to_3bit['N'];
				tmp = tmp << word64_inner_bit_index * 3;

				output_read_3_bit[word64_index] = output_read_3_bit[word64_index] | tmp;

				word64_inner_bit_index++;

				if (word64_inner_bit_index == 21)
				{
					word64_inner_bit_index = 0;
					word64_index++;
				}



			}



			/******************这里要改***********************/








			output_length = output_length + convert_length;
			gap_size = gap_size + convert_length;
			cur_i = i + 1;
			break;
		}

		i++;
	}

	///output_read[output_length] = '\0';


}



inline void init_methylation(Methylation* methy, int size)
{
	(*methy).total_size = size;
	(*methy).current_size = 0;
	(*methy).r_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).sites = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*size);
	/******************这里改了********************/
	///(*methy).r_size = (uint16_t*)malloc(sizeof(uint16_t)*size);
	///(*methy).reads = (char**)malloc(sizeof(char*)*size);

	(*methy).r_real_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).r_size_3_bit = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).reads_3_bit = (bitmapper_bs_iter**)malloc(sizeof(bitmapper_bs_iter*)*size);
	/******************这里改了********************/
	
	
	
	(*methy).genome_cut = genome_cuts;
	(*methy).cut_length = cut_length;
	(*methy).cut_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));
	(*methy).tmp_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));
	int i;
	for (i = 0; i < size; i++)
	{
		/******************这里改了********************/
		///(*methy).reads[i] = (char*)malloc(sizeof(char)*INIT_READ_LENGTH);
		///(*methy).r_size[i] = INIT_READ_LENGTH;
		///(*methy).r_length[i] = 0;
		(*methy).reads_3_bit[i] = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*INIT_READ_3_BIT_LENGTH);
		(*methy).r_size_3_bit[i] = INIT_READ_3_BIT_LENGTH;
		(*methy).r_length[i] = 0;
		(*methy).r_real_length[i] = 0;
		/******************这里改了********************/
		
	}


	memset((*methy).cut_index, 0, sizeof(bitmapper_bs_iter)*((*methy).genome_cut + 1));
}



inline bitmapper_bs_iter get_key_pair_methylation(Pair_Methylation* methy, int i)
{
	/**
	if ((*methy).R[0].r_length[i] & 32)  ///比到正向链
	{
		return (*methy).R[0].sites[i];
	}
	else
	{
		return (*methy).R[1].sites[i];
	}

	**/

	
	if ((*methy).R[1].sites[i] < (*methy).R[0].sites[i])
	{
		return (*methy).R[1].sites[i];
	}
	else
	{
		return (*methy).R[0].sites[i];
	}

}

/**
inline bitmapper_bs_iter get_key_pair_methylation(Pair_Methylation* methy, int i)
{

	if ((*methy).R[1].sites[i] < (*methy).R[0].sites[i])
	{
		
		(*methy).R[0].k_sites = (*methy).R[0].sites[i];
		(*methy).R[0].k_r_length = (*methy).R[0].r_length[i];
		(*methy).R[0].k_r_size = (*methy).R[0].r_size[i];
		(*methy).R[0].k_reads = (*methy).R[0].reads[i];

		(*methy).R[0].sites[i] = (*methy).R[1].sites[i];
		(*methy).R[0].r_length[i] = (*methy).R[1].r_length[i];
		(*methy).R[0].r_size[i] = (*methy).R[1].r_size[i];
		(*methy).R[0].reads[i] = (*methy).R[1].reads[i];

		(*methy).R[1].sites[i] = (*methy).R[0].k_sites;
		(*methy).R[1].r_length[i] = (*methy).R[0].k_r_length;
		(*methy).R[1].r_size[i] = (*methy).R[0].k_r_size;
		(*methy).R[1].reads[i] = (*methy).R[0].k_reads;

	}




	return (*methy).R[0].sites[i];
	
}
**/


inline void init_pair_methylation(Pair_Methylation* methy, int size)
{
	methy->total_size = size;
	methy->current_size = 0;
	(*methy).genome_cut = genome_cuts;
	(*methy).cut_length = cut_length;
	(*methy).cut_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));
	(*methy).tmp_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));


	memset((*methy).cut_index, 0, sizeof(bitmapper_bs_iter)*((*methy).genome_cut + 1));





	/**read1**/
	(*methy).R[0].r_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].sites = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*size);
	/******************这里改了********************/
	/**
	(*methy).R[0].r_size = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].reads = (char**)malloc(sizeof(char*)*size);
	**/
	(*methy).R[0].r_real_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].r_size_3_bit = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].reads_3_bit = (bitmapper_bs_iter**)malloc(sizeof(bitmapper_bs_iter*)*size);
	/******************这里改了********************/
	int i;
	for (i = 0; i < size; i++)
	{
		/******************这里改了********************/
		/**
		(*methy).R[0].reads[i] = (char*)malloc(sizeof(char)*INIT_READ_LENGTH);
		(*methy).R[0].r_size[i] = INIT_READ_LENGTH;
		(*methy).R[0].r_length[i] = 0;
		**/
		(*methy).R[0].reads_3_bit[i] = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*INIT_READ_3_BIT_LENGTH);
		(*methy).R[0].r_size_3_bit[i] = INIT_READ_3_BIT_LENGTH;
		(*methy).R[0].r_length[i] = 0;
		(*methy).R[0].r_real_length[i] = 0;

		/******************这里改了********************/
	}
	/**read1**/





	/**read2**/
	(*methy).R[1].r_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[1].sites = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*size);
	/******************这里改了********************/
	/**
	(*methy).R[1].r_size = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[1].reads = (char**)malloc(sizeof(char*)*size);
	**/
	(*methy).R[1].r_real_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[1].r_size_3_bit = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[1].reads_3_bit = (bitmapper_bs_iter**)malloc(sizeof(bitmapper_bs_iter*)*size);
	/******************这里改了********************/
	for (i = 0; i < size; i++)
	{
		/******************这里改了********************/
		/**
		(*methy).R[1].reads[i] = (char*)malloc(sizeof(char)*INIT_READ_LENGTH);
		(*methy).R[1].r_size[i] = INIT_READ_LENGTH;
		(*methy).R[1].r_length[i] = 0;
		**/
		(*methy).R[1].reads_3_bit[i] = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*INIT_READ_3_BIT_LENGTH);
		(*methy).R[1].r_size_3_bit[i] = INIT_READ_3_BIT_LENGTH;
		(*methy).R[1].r_length[i] = 0;
		(*methy).R[1].r_real_length[i] = 0;
		/******************这里改了********************/
	}
	/**read2**/



}



inline void init_single_methylation(Pair_Methylation* methy, int size)
{
	methy->total_size = size;
	methy->current_size = 0;
	(*methy).genome_cut = genome_cuts;
	(*methy).cut_length = cut_length;
	(*methy).cut_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));
	(*methy).tmp_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*(genome_cuts + 1));


	memset((*methy).cut_index, 0, sizeof(bitmapper_bs_iter)*((*methy).genome_cut + 1));



	/**read1**/
	(*methy).R[0].r_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].sites = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*size);
	/******************这里改了********************/
	/**
	(*methy).R[0].r_size = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].reads = (char**)malloc(sizeof(char*)*size);
	**/
	(*methy).R[0].r_real_length = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].r_size_3_bit = (uint16_t*)malloc(sizeof(uint16_t)*size);
	(*methy).R[0].reads_3_bit = (bitmapper_bs_iter**)malloc(sizeof(bitmapper_bs_iter*)*size);
	/******************这里改了********************/
	int i;
	for (i = 0; i < size; i++)
	{
		/******************这里改了********************/
		/**
		(*methy).R[0].reads[i] = (char*)malloc(sizeof(char)*INIT_READ_LENGTH);
		(*methy).R[0].r_size[i] = INIT_READ_LENGTH;
		(*methy).R[0].r_length[i] = 0;
		**/
		(*methy).R[0].reads_3_bit[i] = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*INIT_READ_3_BIT_LENGTH);
		(*methy).R[0].r_size_3_bit[i] = INIT_READ_3_BIT_LENGTH;
		(*methy).R[0].r_length[i] = 0;
		(*methy).R[0].r_real_length[i] = 0;

		/******************这里改了********************/
	}
	/**read1**/

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

	int read_length;

	///外层循环是每个cut
	for (index = 0; index < (*methy).genome_cut; index++)
	{
		i = (*methy).cut_index[index];
		end_i = (*methy).cut_index[index + 1];
		file = methy_out.files[index];

		for (; i < end_i; i++)
		{
			///fprintf(file, "%llu\t%llu\t%s\n", (*methy).r_length[i], (*methy).sites[i], (*methy).reads[i]);
			fwrite(&((*methy).r_length[i]), sizeof((*methy).r_length[i]), 1, file);
			fwrite(&((*methy).sites[i]), sizeof((*methy).sites[i]), 1, file);
			/**********这里要改*********/
			///read_length = strlen((*methy).reads[i]);
			///fwrite(&read_length, sizeof(read_length), 1, file);
			///fwrite((*methy).reads[i], 1, read_length, file);

			read_length = (*methy).r_real_length[i];
			fwrite(&read_length, sizeof(read_length), 1, file);
			fwrite((*methy).reads_3_bit[i], sizeof(bitmapper_bs_iter),
				(read_length / 21 + (read_length % 21 != 0)), file);
			/**********这里要改*********/

		}
	}

}



inline void output_methylation_pair(Pair_Methylation* methy)
{
	bitmapper_bs_iter i;
	for (i = 0; i <= (*methy).genome_cut; i++)
	{
		methy_out.cut_index[i] += Get_Count((*methy).cut_index, i);
	}

	bitmapper_bs_iter index;
	bitmapper_bs_iter end_i;
	FILE* file;
	int read_length;

	///外层循环是每个cut
	for (index = 0; index < (*methy).genome_cut; index++)
	{
		i = (*methy).cut_index[index];
		end_i = (*methy).cut_index[index + 1];
		file = methy_out.files[index];

		for (; i < end_i; i++)
		{


			///read1的结果
			fwrite(&((*methy).R[0].r_length[i]), sizeof((*methy).R[0].r_length[i]), 1, file);
			fwrite(&((*methy).R[0].sites[i]), sizeof((*methy).R[0].sites[i]), 1, file);
			/**********这里要改*********/
			///read_length = strlen((*methy).R[0].reads[i]);
			///fwrite(&read_length, sizeof(read_length), 1, file);
			///fwrite((*methy).R[0].reads[i], 1, read_length, file);
			read_length = (*methy).R[0].r_real_length[i];
			fwrite(&read_length, sizeof(read_length), 1, file);
			fwrite((*methy).R[0].reads_3_bit[i], sizeof(bitmapper_bs_iter),
				(read_length / 21 + (read_length % 21 != 0)), file);
			/**********这里要改*********/

			///read2的结果
			fwrite(&((*methy).R[1].r_length[i]), sizeof((*methy).R[1].r_length[i]), 1, file);
			fwrite(&((*methy).R[1].sites[i]), sizeof((*methy).R[1].sites[i]), 1, file);
			/**********这里要改*********/
			///read_length = strlen((*methy).R[1].reads[i]);
			///fwrite(&read_length, sizeof(read_length), 1, file);
			///fwrite((*methy).R[1].reads[i], 1, read_length, file);
			read_length = (*methy).R[1].r_real_length[i];
			fwrite(&read_length, sizeof(read_length), 1, file);
			fwrite((*methy).R[1].reads_3_bit[i], sizeof(bitmapper_bs_iter),
				(read_length / 21 + (read_length % 21 != 0)), file);
			/**********这里要改*********/

		}
	}

}



inline int get_cut_id(Methylation* methy, bitmapper_bs_iter site)
{
	return site / (*methy).cut_length;
}


inline int get_cut_id_pair(Pair_Methylation* methy, bitmapper_bs_iter site)
{
	return site / (*methy).cut_length;
}


inline void swap(Methylation* methy, bitmapper_bs_iter i, bitmapper_bs_iter j)
{

	(*methy).k_sites = (*methy).sites[i];
	(*methy).k_r_length = (*methy).r_length[i];
	/************这个要改**************/
	///(*methy).k_r_size = (*methy).r_size[i];
	///(*methy).k_reads = (*methy).reads[i];
	(*methy).k_r_size_3_bit = (*methy).r_size_3_bit[i];
	(*methy).k_reads_3_bit = (*methy).reads_3_bit[i];
	(*methy).k_r_real_length = (*methy).r_real_length[i];
	/************这个要改**************/

	(*methy).sites[i] = (*methy).sites[j];
	(*methy).r_length[i] = (*methy).r_length[j];
	/************这个要改**************/
	///(*methy).r_size[i] = (*methy).r_size[j];
	///(*methy).reads[i] = (*methy).reads[j];
	(*methy).r_size_3_bit[i] = (*methy).r_size_3_bit[j];
	(*methy).reads_3_bit[i] = (*methy).reads_3_bit[j];
	(*methy).r_real_length[i] = (*methy).r_real_length[j];
	/************这个要改**************/

	(*methy).sites[j] = (*methy).k_sites;
	(*methy).r_length[j] = (*methy).k_r_length;
	/************这个要改**************/
	///(*methy).r_size[j] = (*methy).k_r_size;
	///(*methy).reads[j] = (*methy).k_reads;
	(*methy).r_size_3_bit[j] = (*methy).k_r_size_3_bit;
	(*methy).reads_3_bit[j] = (*methy).k_reads_3_bit;
	(*methy).r_real_length[j] = (*methy).k_r_real_length;
	/************这个要改**************/


}


inline void swap_pair(Pair_Methylation* methy, bitmapper_bs_iter i, bitmapper_bs_iter j)
{

	(*methy).R[0].k_sites = (*methy).R[0].sites[i];
	(*methy).R[0].k_r_length = (*methy).R[0].r_length[i];
	/************这个要改**************/
	///(*methy).R[0].k_r_size = (*methy).R[0].r_size[i];
	///(*methy).R[0].k_reads = (*methy).R[0].reads[i];
	(*methy).R[0].k_r_size_3_bit = (*methy).R[0].r_size_3_bit[i];
	(*methy).R[0].k_reads_3_bit = (*methy).R[0].reads_3_bit[i];
	(*methy).R[0].k_r_real_length = (*methy).R[0].r_real_length[i];
	/************这个要改**************/

	(*methy).R[0].sites[i] = (*methy).R[0].sites[j];
	(*methy).R[0].r_length[i] = (*methy).R[0].r_length[j];
	/************这个要改**************/
	//(*methy).R[0].r_size[i] = (*methy).R[0].r_size[j];
	//(*methy).R[0].reads[i] = (*methy).R[0].reads[j];
	(*methy).R[0].r_size_3_bit[i] = (*methy).R[0].r_size_3_bit[j];
	(*methy).R[0].reads_3_bit[i] = (*methy).R[0].reads_3_bit[j];
	(*methy).R[0].r_real_length[i] = (*methy).R[0].r_real_length[j];
	/************这个要改**************/

	(*methy).R[0].sites[j] = (*methy).R[0].k_sites;
	(*methy).R[0].r_length[j] = (*methy).R[0].k_r_length;
	/************这个要改**************/
	///(*methy).R[0].r_size[j] = (*methy).R[0].k_r_size;
	///(*methy).R[0].reads[j] = (*methy).R[0].k_reads;
	(*methy).R[0].r_size_3_bit[j] = (*methy).R[0].k_r_size_3_bit;
	(*methy).R[0].reads_3_bit[j] = (*methy).R[0].k_reads_3_bit;
	(*methy).R[0].r_real_length[j] = (*methy).R[0].k_r_real_length;
	/************这个要改**************/





	/************这个要改**************/
	(*methy).R[1].k_sites = (*methy).R[1].sites[i];
	(*methy).R[1].k_r_length = (*methy).R[1].r_length[i];
	///(*methy).R[1].k_r_size = (*methy).R[1].r_size[i];
	///(*methy).R[1].k_reads = (*methy).R[1].reads[i];
	(*methy).R[1].k_r_size_3_bit = (*methy).R[1].r_size_3_bit[i];
	(*methy).R[1].k_reads_3_bit = (*methy).R[1].reads_3_bit[i];
	(*methy).R[1].k_r_real_length = (*methy).R[1].r_real_length[i];
	/************这个要改**************/


	/************这个要改**************/
	(*methy).R[1].sites[i] = (*methy).R[1].sites[j];
	(*methy).R[1].r_length[i] = (*methy).R[1].r_length[j];
	///(*methy).R[1].r_size[i] = (*methy).R[1].r_size[j];
	///(*methy).R[1].reads[i] = (*methy).R[1].reads[j];
	(*methy).R[1].r_size_3_bit[i] = (*methy).R[1].r_size_3_bit[j];
	(*methy).R[1].reads_3_bit[i] = (*methy).R[1].reads_3_bit[j];
	(*methy).R[1].r_real_length[i] = (*methy).R[1].r_real_length[j];
	/************这个要改**************/


	/************这个要改**************/
	(*methy).R[1].sites[j] = (*methy).R[1].k_sites;
	(*methy).R[1].r_length[j] = (*methy).R[1].k_r_length;
	///(*methy).R[1].r_size[j] = (*methy).R[1].k_r_size;
	///(*methy).R[1].reads[j] = (*methy).R[1].k_reads;
	(*methy).R[1].r_size_3_bit[j] = (*methy).R[1].k_r_size_3_bit;
	(*methy).R[1].reads_3_bit[j] = (*methy).R[1].k_reads_3_bit;
	(*methy).R[1].r_real_length[j] = (*methy).R[1].k_r_real_length;
	/************这个要改**************/


}







inline void assign_cuts_pair(Pair_Methylation* methy)
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


			///tmp_index = get_cut_id_pair(methy, (*methy).sites[i]);
			tmp_index = get_cut_id_pair(methy, get_key_pair_methylation(methy, i));

			if (tmp_index == index)
			{
				i++;
				(*methy).tmp_index[index]++; ///这个可以不加
			}
			else
			{
				swap_pair(methy, i, (*methy).tmp_index[tmp_index]);
				(*methy).tmp_index[tmp_index]++;
			}
		}

	}

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



inline void clear_methylation_pair(Pair_Methylation* methy)
{
	(*methy).current_size = 0;


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


void methy_extract(int thread_id, char* file_name, int* need_context);

void methy_extract_PE(int thread_id, char* file_name, int PE_distance, int* need_context);


///min一定是个坐标, max一定是个坐标+长度
inline long long calculate_TLEN(long long min1, long long len1, long long min2, long long len2)
{

	long long min = min1;
	long long max = min1 + len1 - 1;

	if (min1 > min2)
		min = min2;

	if (max < min2 + len2 - 1)
		max = min2 + len2 - 1;

	return max - min + 1;
}


///min一定是个坐标, max一定是个坐标+长度
inline long long calculate_last_base_pos(long long min1, long long len1, long long min2, long long len2)
{

	long long min = min1;
	long long max = min1 + len1 - 1;

	if (min1 > min2)
		min = min2;

	if (max < min2 + len2 - 1)
		max = min2 + len2 - 1;

	///return max - min + 1;
	return max;
}


///min一定是个坐标, max一定是个坐标+长度
inline int over_lap_length(long long min1, long long len1, long long min2, long long len2,char* read2)
{

	long long min = min1;
	long long max = min1 + len1;

	if (min1 > min2)
		min = min2;

	if (min1 + len1 < min2 + len2)
		max = min2 + len2;

	long long overlap_length = max - min - len1 - len2;

	long long region_start2;
	int region_length2;
	///如果不重叠是正值或者0, 重叠是负值
	if (overlap_length < 0)
	{
		overlap_length = overlap_length*-1;

		if (min1<=min2)
		{
			region_start2 = 0;
			region_length2 = overlap_length;
		}
		else if (min1>min2)
		{
			region_start2 = min1 - min2;
			region_length2 = overlap_length;
		}

		memset(read2 + region_start2, 'N', region_length2);

		return 1;
	}

	return 0;

}


void out_paired_distance_statistic();

void get_genome_cuts(char* file_name);
















typedef struct
{
	///只负责尾部区域的锁
	pthread_mutex_t tail_Mutex;
	bitmapper_bs_iter head_zone_start;
	bitmapper_bs_iter head_zone_end;
	bitmapper_bs_iter tail_zone_start;
	bitmapper_bs_iter tail_zone_end;
}
Block_overlap_Mutex;

typedef struct
{
	char* buffer;
	long long size;
	long long length;

} Output_buffer_sub_block;



typedef struct
{
	int number_of_intervals;
	Pair_Methylation* intervals;
	bitmapper_bs_iter each_interval_length;
	bitmapper_bs_iter* each_buffer_interval_start;
	bitmapper_bs_iter* each_buffer_interval_end;
	Output_buffer_sub_block* CpG_buffer;
	Output_buffer_sub_block* CHG_buffer;
	Output_buffer_sub_block* CHH_buffer;




	uint16_t** M_methy;
	uint16_t** M_total;
	uint16_t** M_correct_methy;
	uint16_t** M_correct_total;
	char** M_ref_genome;
	long long* M_start_pos;
	long long* M_end_pos;
	int need_context[3];




	pthread_mutex_t* Mutex;
	pthread_cond_t* stallCond;
	pthread_cond_t* flushCond;

	pthread_mutex_t input_thread_Mutex;
	pthread_cond_t input_thread_flushCond;
	pthread_mutex_t abnormal_file_Mutex;
	pthread_mutex_t overlap_file_Mutex;

	pthread_mutex_t main_thread_Mutex;
	pthread_cond_t main_thread_flushCond;


	pthread_mutex_t* process_Mutex;
	pthread_cond_t* process_stallCond;

	int end;
	int all_end;

	int all_completed;

	pthread_mutex_t all_completed_Mutex;

	bitmapper_bs_iter total_start;
	bitmapper_bs_iter total_end;



	/*************methy_extract**************/
	pthread_mutex_t all_completed_Mutex_methy;
	pthread_mutex_t main_thread_Mutex_methy;
	pthread_cond_t main_thread_flushCond_methy;
	int all_completed_methy;

	pthread_mutex_t* process_Mutex_methy;
	pthread_cond_t* process_stallCond_methy;

	/*************methy_extract**************/


} Inpute_PE_methy_alignment_buffer;



typedef struct
{
	bitmapper_bs_iter* buffer;
	int size;

}buffer_3_bit;


typedef struct
{
	FILE* read_file;
	FILE* abnormal_file;
	FILE* overlap_file;
	char* ref_genome;
	uint16_t* methy;
	uint16_t* total;
	uint16_t* correct_methy;
	uint16_t* correct_total;
	bitmapper_bs_iter start_pos;
	bitmapper_bs_iter end_pos;
	bitmapper_bs_iter length;
	deduplicate_PE PE_de;
	int PE_distance;
	bitmapper_bs_iter extra_length;
	long long* total_read;
	long long* duplicate_read; 
	long long* unique_read;
	long long* abnormal_read;
}PE_methy_parameters;


inline void init_PE_methy_parameters(PE_methy_parameters *parameters, int PE_distance, bitmapper_bs_iter cut_length)
{
	parameters->PE_distance = PE_distance;
	parameters->extra_length = SEQ_MAX_LENGTH * 2 + parameters->PE_distance;


	parameters->methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->correct_methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->correct_total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->ref_genome = (char*)malloc(sizeof(char)*(cut_length + parameters->extra_length));
	///按道理来说, PE_de的长度不该有这个extra_length, 但是多加一点以免出错吧
	init_deduplicate_PE(&(parameters->PE_de), PE_distance, cut_length + parameters->extra_length);


		
	parameters->total_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);
	parameters->duplicate_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);
	parameters->unique_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);
	parameters->abnormal_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);
	

	int i = 0;

	for (i = 0; i < THREAD_COUNT; i++)
	{
		parameters->total_read[i] = 0;
		parameters->duplicate_read[i] = 0;
		parameters->unique_read[i] = 0;
		parameters->abnormal_read[i] = 0;
	}
}




inline void init_single_methy_parameters(PE_methy_parameters *parameters, bitmapper_bs_iter cut_length)
{
	parameters->extra_length = SEQ_MAX_LENGTH * 2;


	parameters->methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->correct_methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->correct_total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + parameters->extra_length));
	parameters->ref_genome = (char*)malloc(sizeof(char)*(cut_length + parameters->extra_length));



	parameters->total_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);
	parameters->duplicate_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);
	parameters->unique_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);
	parameters->abnormal_read = (long long*)malloc(sizeof(long long)*THREAD_COUNT);


	int i = 0;

	for (i = 0; i < THREAD_COUNT; i++)
	{
		parameters->total_read[i] = 0;
		parameters->duplicate_read[i] = 0;
		parameters->unique_read[i] = 0;
		parameters->abnormal_read[i] = 0;
	}
}


#define Init_buffer_3_bit(b) b.size = INIT_READ_3_BIT_LENGTH; b.buffer = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*b.size);


void init_methy_input_buffer(long long each_sub_buffer_size, int number_of_threads,
	bitmapper_bs_iter total_start, bitmapper_bs_iter total_end, int is_update, int is_paired);

void* input_methy_muti_threads(void* arg);
void* input_methy_muti_threads_single_end(void*);

void methy_extract_PE_mutiple_thread(int thread_id, char* file_name, int PE_distance, int* need_context);

void methy_extract_mutiple_thread(int thread_id, char* file_name, int* need_context);


#endif
