/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<stdint.h>
#include <math.h>
#include <dirent.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <sys/time.h>
#include <algorithm>
#include <immintrin.h>
#include <iostream>
#include <sstream>
#include <map>
#include <fstream>  
#include "Auxiliary.h"
#include "Process_Reads.h"
#include "Index.h"
#include "Process_sam_out.h"
#include "Levenshtein_Cal.h"
#include "SAM_queue.h"

#include "Ref_Genome.h"


_rg_name_l  *_ih_refGenName;
int refChromeCont;

char *versionN = "1.0.0.0";
long long mappingCnt[MAX_Thread];
unsigned int done;
long long mappedSeqCnt[MAX_Thread];
long long completedSeqCnt[MAX_Thread];
///char *Mapped_File;
char *_msf_refGen = NULL;
unsigned char *_BS_refGen = NULL;
bitmapper_bs_iter _msf_refGenLength = 0;
bitmapper_bs_iter _msf_refGenLength_2_bit = 0;
Commpress_HashTable  chhy_HashTable;
bitmapper_bs_iter total_SA_length;

long long unique_mapped_read[MAX_Thread];
long long ambiguous_mapped_read[MAX_Thread];
///long long unmapped_read[MAX_Thread];

int batch_read_size = 50000;

OPT_FIELDS **_chhy_optionalFields;

FILE *schema_out_fp;

pthread_mutex_t readinputMutex;
pthread_mutex_t queueMutex;
pthread_mutex_t terminateMutex;
pthread_cond_t flushCond;
pthread_cond_t readinputflushCond;
pthread_cond_t stallCond;
pthread_cond_t readinputstallCond;
///pthread_cond_t doneMutex;
pthread_mutex_t doneMutex;

int outputMaxOutstanding=100000;
int queue_length=0;

tmp_result* output_stack;

tmp_result_string* output_stack_string;

QUEUE reads_queue;

char four_ATGC[256][4];
char rc_four_ATGC[256][4];


SAM_Queue sam_output_queue;

int all_read_end=0;

long long enq_i=0;

int ma = 1;
int mp = -3;
int gap_open = -5;
int gap_extension = -3;
int is_gap_open = 1;

char map_cigar[9];
///double all_debug_time = 0;
///double total = 0;

int sub_block_inner_size = 500;
///实际队列长度都要乘以线程数目
int sub_block_number = 100;
int ouput_buffer_size = 100;










int bitmapper_bs_iter_compareEntrySize(const void *a, const void *b) {
	bitmapper_bs_iter a_list = *(bitmapper_bs_iter *)a;
	bitmapper_bs_iter b_list = *(bitmapper_bs_iter *)b;
	if (a_list>b_list)
		return 1;
	else if (a_list<b_list)
		return -1;
	else
		return 0;

}


inline void load_batch_paired_read(Read* read_batch1, Read* read_batch2, int batch_read_size,
	int* return_file_flag, int* obtain_reads_num)
{
	int inner_i = 0;
	int file_flag = 1;

	while (inner_i<batch_read_size)
	{
		///file_flag = inputReads_single_directly(&read_batch[inner_i]);

		file_flag = inputReads_paired_directly(&read_batch1[inner_i], &read_batch2[inner_i]);

		if (file_flag == 1)
		{
			inner_i++;
		}
		else if (file_flag == 0)
		{
			break;
		}
	}

	*return_file_flag = file_flag;
	*obtain_reads_num = inner_i;

}


void get_mapping_informations
(long long* number_of_read, long long* number_of_unique_mapped_read,
long long* number_of_ambiguous_mapped_read, long long* number_of_unmapped_read)
{
	*number_of_read = 0;
	*number_of_unique_mapped_read = 0;
	*number_of_ambiguous_mapped_read = 0;

	int i = 0;
	for (i = 0; i<THREAD_COUNT; i++)
	{
		*number_of_read += completedSeqCnt[i];
		*number_of_unique_mapped_read += unique_mapped_read[i];
		*number_of_ambiguous_mapped_read += ambiguous_mapped_read[i];
		///*number_of_unmapped_read += unmapped_read[i];

	}

	*number_of_unmapped_read = *number_of_read - *number_of_unique_mapped_read - *number_of_ambiguous_mapped_read;
}



bool compare_seed_votes(const seed_votes& s1, const seed_votes& s2)
{
	return s1.vote > s2.vote; //从大到小排序
}


bool compare_site_err(const seed_votes& s1, const seed_votes& s2)
{
	return s1.err < s2.err; //从大到小排序
}


int lists_compareEntrySize(const void *a, const void *b) {
    Candidate a_list=*(Candidate *)a;
    Candidate b_list=*(Candidate *)b;
    if(a_list.Cnt>b_list.Cnt)
        return 1;
    else if(a_list.Cnt<b_list.Cnt)
        return -1;
    else
        return 0;

}

int unsigned_int_compareEntrySize(const void *a, const void *b) {
	unsigned int a_list = *(unsigned int *)a;
	unsigned int b_list = *(unsigned int *)b;
	if (a_list>b_list)
		return 1;
	else if (a_list<b_list)
		return -1;
	else
		return 0;

}

unsigned int get_total_reads_number()
{
    return enq_i;
}



void Prepare_alignment(char *genFileName,_rg_name_l *chhy_ih_refGenName,int chhy_refChromeCont)
{


	int ijkijk = 0;

    done=0;
    queue_length=0;

    int offset=0;
	int sqrt_2=1;


	while(1)
	{
        if(SEQ_LENGTH<=(sqrt_2-1))
        {
            break;
        }
        offset++;
        sqrt_2=sqrt_2*2;

	}

	if(offset>=(2*thread_e+2))
    {
        each_length=offset;
    }
    else
    {
        each_length=2*thread_e+2;
    }

    int i;




    _ih_refGenName=chhy_ih_refGenName;
    refChromeCont=chhy_refChromeCont;


	_chhy_optionalFields = (OPT_FIELDS **)malloc(THREAD_COUNT * sizeof(OPT_FIELDS*));
    if (is_pairedEnd)
    {
        for(i=0;i<THREAD_COUNT;++i)
        {
			_chhy_optionalFields[i] = (OPT_FIELDS*)malloc(8 * sizeof(OPT_FIELDS));
        }
    }
    else
    {
        for(i=0;i<THREAD_COUNT;++i)
        {
			_chhy_optionalFields[i] = (OPT_FIELDS*)malloc(2 * sizeof(OPT_FIELDS));
        }
    }




	_BS_refGen = getRefGenome();
    _msf_refGenLength = getRefGenomeLength();
	_msf_refGenLength_2_bit = getRefGenomeLength_2bit();
	
	/**
	fprintf(stderr, "_msf_refGenLength: %llu\n", _msf_refGenLength);
	fprintf(stderr, "_msf_refGenLength_2_bit: %llu\n", _msf_refGenLength_2_bit);
	**/

    chhy_HashTable=getHashTable();



	char ATGC[4];

	ATGC[0] = 'A';
	ATGC[1] = 'C';
	ATGC[2] = 'G';
	ATGC[3] = 'T';

	unsigned int ijkijk_mask = (unsigned int)3;

	for (ijkijk = 0; ijkijk<256; ijkijk++)
	{

		four_ATGC[ijkijk][0] = ATGC[(ijkijk >> 6)];
		four_ATGC[ijkijk][1] = ATGC[(ijkijk >> 4)&3];
		four_ATGC[ijkijk][2] = ATGC[(ijkijk >> 2)&3];
		four_ATGC[ijkijk][3] = ATGC[(ijkijk)&3];
	}

	char rc_ATGC[256];


	rc_ATGC['A'] = 'T';
	rc_ATGC['C'] = 'G';
	rc_ATGC['G'] = 'C';
	rc_ATGC['T'] = 'A';

	for (ijkijk = 0; ijkijk<256; ijkijk++)
	{

		rc_four_ATGC[ijkijk][0] = rc_ATGC[four_ATGC[ijkijk][3]];
		rc_four_ATGC[ijkijk][1] = rc_ATGC[four_ATGC[ijkijk][2]];
		rc_four_ATGC[ijkijk][2] = rc_ATGC[four_ATGC[ijkijk][1]];
		rc_four_ATGC[ijkijk][3] = rc_ATGC[four_ATGC[ijkijk][0]];
	}




	/**
	for (ijkijk = 0; ijkijk < 256; ijkijk++)
	{
		std::cerr << "i    : " << ijkijk << std::endl;
		std::cerr << four_ATGC[ijkijk][0] << four_ATGC[ijkijk][1] << four_ATGC[ijkijk][2] << four_ATGC[ijkijk][3] << std::endl;

		std::cerr << "i(rc): " << ijkijk << std::endl;
		std::cerr << rc_four_ATGC[ijkijk][0] << rc_four_ATGC[ijkijk][1] << rc_four_ATGC[ijkijk][2] << rc_four_ATGC[ijkijk][3] << std::endl;

	}
	**/

	schema_out_fp = get_Ouput_Dec();


	map_cigar[0] = 'M';
	map_cigar[1] = 'M';
	map_cigar[2] = 'D';
	map_cigar[3] = 'I';

	map_cigar[5] = 'M';
	map_cigar[6] = 'M';
	map_cigar[7] = 'D';
	map_cigar[8] = 'I';


}







inline void reverse_and_adjust_site
(bitmapper_bs_iter* locates, bitmapper_bs_iter* locates_occ, bitmapper_bs_iter seed_length, bitmapper_bs_iter seed_offset)
{
	/**
	bitmapper_bs_iter i = 0, avaiable_i = 0;
	bitmapper_bs_iter tmp_site;
	**/
	bitmapper_bs_iter i = 0;

	for (i = 0; i < *locates_occ; i++)
	{
		///locates[i] = total_SA_length - 1 - (locates[i] + seed_length - 1) - seed_offset;
		locates[i] = total_SA_length - locates[i] - seed_length - seed_offset;
		/**
		tmp_site = locates[i] + seed_length + seed_offset;
		if (total_SA_length >= tmp_site)
		{
			locates[avaiable_i] = total_SA_length - tmp_site;
			avaiable_i++;
		}
		**/

	}

	///*locates_occ = avaiable_i;
	*locates_occ = i;
}



inline void generate_candidate_votes_shift(bitmapper_bs_iter* candidates, bitmapper_bs_iter candidates_length,
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length, bitmapper_bs_iter error_threshold)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter pre_candidate;
	bitmapper_bs_iter candidate_votes_i = 0;
	bitmapper_bs_iter vote;

	pre_candidate = candidates[0];

	i = 1;

	vote = 1;

	candidate_votes_i = 0;



	while (pre_candidate < error_threshold && i < candidates_length)
	{

		if (candidates[i] == pre_candidate)
		{
			i++;
			vote++;
		}
		else
		{
			candidate_votes[candidate_votes_i].site = 0;
			candidate_votes[candidate_votes_i].vote = vote;
			candidate_votes_i++;



			vote = 1;
			pre_candidate = candidates[i];

			i++;
		}

	}





	while (i < candidates_length)
	{

		if(candidates[i] == pre_candidate)
		{
			i++;
			vote++;
		}
		else
		{

			candidate_votes[candidate_votes_i].site = pre_candidate - error_threshold;
			candidate_votes[candidate_votes_i].vote = vote;
			candidate_votes_i++;



			vote = 1;
			pre_candidate = candidates[i];

			i++;
		}

	}


	///最后要单独处理
	candidate_votes[candidate_votes_i].site = 0;
	if (pre_candidate >= error_threshold)
	{
		candidate_votes[candidate_votes_i].site = pre_candidate - error_threshold;
	}
	candidate_votes[candidate_votes_i].vote = vote;
	candidate_votes_i++;



	(*candidate_votes_length) = candidate_votes_i;


}

inline int select_suit_candidates(bitmapper_bs_iter site, seed_votes* result2_array, int best_mapp_occ2,
	int inner_maxDistance_pair, int inner_minDistance_pair, int* next_start_index)
{
	long long distance;

	

	///for (int i = 0; i < best_mapp_occ2; i++)
	for (int i = *next_start_index; i < best_mapp_occ2; i++)
	{
		if (result2_array[i].site > site)
		{
			distance = result2_array[i].site - site;

			if (distance > inner_maxDistance_pair)
			{
				return 0;
			}


			if (distance <= inner_maxDistance_pair
				&&
				distance >= inner_minDistance_pair)
			{
				return 1;
			}


		}
		else
		{
			distance = site - result2_array[i].site;

			if (distance > inner_maxDistance_pair)
			{
				(*next_start_index) = i + 1;
			}
			else if(distance >= inner_minDistance_pair)
			{
				return 1;
			}

		}

		

	}

	return 0;
}



inline int debug_select_suit_candidates(bitmapper_bs_iter site, map_result* result2_array, int best_mapp_occ2,
	int inner_maxDistance_pair, int inner_minDistance_pair)
{
	long long distance;



	for (int i = 0; i < best_mapp_occ2; i++)
	{
		if (result2_array[i].origin_site > site)
		{
			distance = result2_array[i].origin_site - site;

			if (distance > inner_maxDistance_pair)
			{
				return 0;
			}


			if (distance <= inner_maxDistance_pair
				&&
				distance >= inner_minDistance_pair)
			{
				return 1;
			}


		}
		else
		{
			distance = site - result2_array[i].origin_site;

			if (distance > inner_maxDistance_pair)
			{
				;
				//(*next_start_index) = i + 1;
			}
			else if (distance >= inner_minDistance_pair)
			{
				return 1;
			}

		}



	}

	return 0;
}






inline void generate_candidate_votes_shift_filter(bitmapper_bs_iter* candidates, bitmapper_bs_iter candidates_length,
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length, bitmapper_bs_iter error_threshold, 
	seed_votes* result2_array, int best_mapp_occ2, int inner_maxDistance_pair, int inner_minDistance_pair)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter pre_candidate;
	bitmapper_bs_iter candidate_votes_i = 0;
	bitmapper_bs_iter vote;
	int next_start_index=0;

	pre_candidate = candidates[0];

	i = 1;

	vote = 1;

	candidate_votes_i = 0;



	while (pre_candidate < error_threshold && i < candidates_length)
	{

		if (candidates[i] == pre_candidate)
		{
			i++;
			vote++;
		}
		else
		{
			candidate_votes[candidate_votes_i].site = 0;
			if (select_suit_candidates(candidate_votes[candidate_votes_i].site, result2_array, best_mapp_occ2,
				inner_maxDistance_pair, inner_minDistance_pair, &next_start_index))
			{
				candidate_votes[candidate_votes_i].vote = vote;
				candidate_votes_i++;
			}



				

		
			vote = 1;
			pre_candidate = candidates[i];

			i++;
		}

	}





	while (i < candidates_length)
	{

		if (candidates[i] == pre_candidate)
		{
			i++;
			vote++;
		}
		else
		{

			candidate_votes[candidate_votes_i].site = pre_candidate - error_threshold;
			if (select_suit_candidates(candidate_votes[candidate_votes_i].site, result2_array, best_mapp_occ2,
				inner_maxDistance_pair, inner_minDistance_pair, &next_start_index))
			{
				candidate_votes[candidate_votes_i].vote = vote;
				candidate_votes_i++;
			}


			
			



			vote = 1;
			pre_candidate = candidates[i];

			i++;
		}

	}


	///最后要单独处理
	candidate_votes[candidate_votes_i].site = 0;
	if (pre_candidate >= error_threshold)
	{
		candidate_votes[candidate_votes_i].site = pre_candidate - error_threshold;
	}

	if (select_suit_candidates(candidate_votes[candidate_votes_i].site, result2_array, best_mapp_occ2,
		inner_maxDistance_pair, inner_minDistance_pair, &next_start_index))
	{
		candidate_votes[candidate_votes_i].vote = vote;
		candidate_votes_i++;
	}



	(*candidate_votes_length) = candidate_votes_i;


}





inline void get_actuall_genome(char* tmp_ref, bitmapper_bs_iter start_site, bitmapper_bs_iter length)
{
	///double debug_start_time = Get_T();

	bitmapper_bs_iter convert_i = 0;

	bitmapper_bs_iter index = (start_site >> 2);

	bitmapper_bs_iter last = (start_site & 3);





	///if (index + convert_t_length > _msf_refGenLength_2_bit)
	if (start_site + length > _msf_refGenLength)
	{
		memset(tmp_ref, 0, length);

		return;

	}



	///这个实际上是第一8-bit里面有多少个我们要的字符
	convert_i = 4 - last;

	///four_i++;

	memcpy(tmp_ref, four_ATGC[_BS_refGen[index]] + last, convert_i);

	index++;

	while (convert_i<length)
	{


		memcpy(tmp_ref + convert_i, four_ATGC[_BS_refGen[index]], 4);

		convert_i = convert_i + 4;

		index++;
		
	}



	tmp_ref[convert_i] = '\0';

	///all_debug_time = all_debug_time + Get_T() - debug_start_time;

}










inline void get_actuall_rc_genome(char* tmp_ref, bitmapper_bs_iter rc_start_site, bitmapper_bs_iter length)
{

	////double debug_start_time = Get_T();
	

	bitmapper_bs_iter end_site = _msf_refGenLength - rc_start_site - 1;

	bitmapper_bs_iter convert_i = 0;

	bitmapper_bs_iter index = (end_site >> 2);

	bitmapper_bs_iter last = (end_site & 3);


	if (end_site < length - 1
		||
		index >= _msf_refGenLength_2_bit)
	{
		memset(tmp_ref, 0, length);

		return;

	}


	///这个实际上是第一8-bit里面有多少个我们要的字符
	convert_i = last + 1;

	memcpy(tmp_ref, rc_four_ATGC[_BS_refGen[index]] + 4 - convert_i, convert_i);


	///注意这里是--
	index--;

	while (convert_i < length)
	{

		memcpy(tmp_ref + convert_i, rc_four_ATGC[_BS_refGen[index]], 4);

		convert_i = convert_i + 4;

		///注意这里是--
		index--;


	}



	tmp_ref[convert_i] = '\0';

	///all_debug_time = all_debug_time + Get_T() - debug_start_time;

}

















inline void map_candidate_votes_mutiple_cut_end_to_end_4(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[4];
	unsigned int return_sites_error[4];

	char* t[4];



	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;
	t[2] = t[1] + p_length + 32;
	t[3] = t[2] + p_length + 32;
	

	i = 0;

	///这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	(*min_err_index) = -1;
	bitmapper_bs_iter min_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;

	while (i + 4 <= *candidate_votes_length)
	{

		if (candidate_votes[i].site < _msf_refGenLength)
		{
			get_actuall_genome(t[0], candidate_votes[i].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[0], candidate_votes[i].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 1].site < _msf_refGenLength)
		{
			get_actuall_genome(t[1], candidate_votes[i + 1].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[1], candidate_votes[i + 1].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 2].site < _msf_refGenLength)
		{
			get_actuall_genome(t[2], candidate_votes[i + 2].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[2], candidate_votes[i + 2].site - _msf_refGenLength, p_length);
		}



		if (candidate_votes[i + 3].site < _msf_refGenLength)
		{
			get_actuall_genome(t[3], candidate_votes[i + 3].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[3], candidate_votes[i + 3].site - _msf_refGenLength, p_length);
		}





		BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i + 1].end_site = return_sites[1];
		candidate_votes[i + 2].end_site = return_sites[2];
		candidate_votes[i + 3].end_site = return_sites[3];

		candidate_votes[i].err = return_sites_error[0];
		candidate_votes[i + 1].err = return_sites_error[1];
		candidate_votes[i + 2].err = return_sites_error[2];
		candidate_votes[i + 3].err = return_sites_error[3];


		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[0] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*min_err) = return_sites_error[0];
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}


		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[1] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}



		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[2] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*min_err) = return_sites_error[2];
			(*min_err_index) = i + 2;
			min_err_site = tmp_min_err_site;
		}



		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[3] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*min_err) = return_sites_error[3];
			(*min_err_index) = i + 3;
			min_err_site = tmp_min_err_site;
		}




		i = i + 4;

		/**
		if ((*min_err) == 0)
		{

		///fprintf(stderr, "ERROR: (*min_err) cannot be zero\n");

		(*candidate_votes_length) = i;



		candidate_votes[0].err = candidate_votes[(*min_err_index)].err;
		candidate_votes[0].end_site = candidate_votes[(*min_err_index)].end_site;
		candidate_votes[0].site = candidate_votes[(*min_err_index)].site;


		return;
		///break;
		}
		**/






	}

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);

	int last = i;
	for (; last < *candidate_votes_length; last++)
	{


		if (candidate_votes[last].site < _msf_refGenLength)
		{
			get_actuall_genome(t[last - i], candidate_votes[last].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[last - i], candidate_votes[last].site - _msf_refGenLength, p_length);
		}



	}


	int last_length = *candidate_votes_length - i;

	if (last_length == 1)
	{
		candidate_votes[i].end_site =
			BS_Reserve_Banded_BPM(t[0], p_length, read, t_length, error_threshold, &(candidate_votes[i].err));



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (candidate_votes[i].err == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		int inner_i = 0;

		for (inner_i = 0; inner_i < last_length; inner_i++)
		{

			candidate_votes[i + inner_i].end_site = return_sites[inner_i];
			candidate_votes[i + inner_i].err = return_sites_error[inner_i];

			tmp_min_err_site = candidate_votes[i + inner_i].site + candidate_votes[i + inner_i].end_site;
			///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
			if (return_sites_error[inner_i] == (*min_err)
				&&
				min_err_site != tmp_min_err_site)
			{
				(*min_err_index) = -2;
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*min_err) = return_sites_error[inner_i];
				(*min_err_index) = i + inner_i;
				min_err_site = tmp_min_err_site;
			}

		}

	}


	if ((*min_err_index) >= 0)
	{

		candidate_votes[0].err = candidate_votes[(*min_err_index)].err;
		candidate_votes[0].end_site = candidate_votes[(*min_err_index)].end_site;
		candidate_votes[0].site = candidate_votes[(*min_err_index)].site;
	}


}







inline void map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* best_mapping_occ,
	char* name)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[4];
	unsigned int return_sites_error[4];

	char* t[4];



	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;
	t[2] = t[1] + p_length + 32;
	t[3] = t[2] + p_length + 32;


	i = 0;

	///这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	///(*min_err_index) = -1;
	(*best_mapping_occ) = 0;
	bitmapper_bs_iter pre_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;



	while (i + 4 <= *candidate_votes_length)
	{

		if (candidate_votes[i].site < _msf_refGenLength)
		{
			get_actuall_genome(t[0], candidate_votes[i].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[0], candidate_votes[i].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 1].site < _msf_refGenLength)
		{
			get_actuall_genome(t[1], candidate_votes[i + 1].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[1], candidate_votes[i + 1].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 2].site < _msf_refGenLength)
		{
			get_actuall_genome(t[2], candidate_votes[i + 2].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[2], candidate_votes[i + 2].site - _msf_refGenLength, p_length);
		}



		if (candidate_votes[i + 3].site < _msf_refGenLength)
		{
			get_actuall_genome(t[3], candidate_votes[i + 3].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[3], candidate_votes[i + 3].site - _msf_refGenLength, p_length);
		}





		BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i + 1].end_site = return_sites[1];
		candidate_votes[i + 2].end_site = return_sites[2];
		candidate_votes[i + 3].end_site = return_sites[3];

		candidate_votes[i].err = return_sites_error[0];
		candidate_votes[i + 1].err = return_sites_error[1];
		candidate_votes[i + 2].err = return_sites_error[2];
		candidate_votes[i + 3].err = return_sites_error[3];
		


		
		


		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		if (return_sites_error[0] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		if (return_sites_error[1] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 1].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 1].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 1].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		if (return_sites_error[2] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 2].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 2].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 2].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		if (return_sites_error[3] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 3].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 3].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 3].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;


		i = i + 4;

	}

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);

	int last = i;
	for (; last < *candidate_votes_length; last++)
	{


		if (candidate_votes[last].site < _msf_refGenLength)
		{
			get_actuall_genome(t[last - i], candidate_votes[last].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[last - i], candidate_votes[last].site - _msf_refGenLength, p_length);
		}



	}


	int last_length = *candidate_votes_length - i;

	if (last_length == 1)
	{
		candidate_votes[i].end_site =
			BS_Reserve_Banded_BPM(t[0], p_length, read, t_length, error_threshold, &(candidate_votes[i].err));



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		if (candidate_votes[i].err <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		int inner_i = 0;

		for (inner_i = 0; inner_i < last_length; inner_i++)
		{

			candidate_votes[i + inner_i].end_site = return_sites[inner_i];
			candidate_votes[i + inner_i].err = return_sites_error[inner_i];


			tmp_min_err_site = candidate_votes[i + inner_i].site + candidate_votes[i + inner_i].end_site;
			if (return_sites_error[inner_i] <= error_threshold
				&&
				pre_err_site != tmp_min_err_site)
			{
				candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + inner_i].site;
				candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + inner_i].err;
				candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + inner_i].end_site;
				(*best_mapping_occ)++;
			}
			pre_err_site = tmp_min_err_site;


		}

	}

}


///就是加了一个返回值，返回值记录best mapping的下标罢了
inline void map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* best_mapping_occ,
	char* name)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[8];
	unsigned int return_sites_error[8];
	/**
	int best_mapping_index = -1;
	int current_error;

	((unsigned int)-1) - 1
		**/

	char* t[8];



	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;
	t[2] = t[1] + p_length + 32;
	t[3] = t[2] + p_length + 32;
	t[4] = t[3] + p_length + 32;
	t[5] = t[4] + p_length + 32;
	t[6] = t[5] + p_length + 32;
	t[7] = t[6] + p_length + 32;


	i = 0;

	///这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	///(*min_err_index) = -1;
	(*best_mapping_occ) = 0;
	bitmapper_bs_iter pre_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;



	while (i + 8 <= *candidate_votes_length)
	{

		if (candidate_votes[i].site < _msf_refGenLength)
		{
			get_actuall_genome(t[0], candidate_votes[i].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[0], candidate_votes[i].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 1].site < _msf_refGenLength)
		{
			get_actuall_genome(t[1], candidate_votes[i + 1].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[1], candidate_votes[i + 1].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 2].site < _msf_refGenLength)
		{
			get_actuall_genome(t[2], candidate_votes[i + 2].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[2], candidate_votes[i + 2].site - _msf_refGenLength, p_length);
		}



		if (candidate_votes[i + 3].site < _msf_refGenLength)
		{
			get_actuall_genome(t[3], candidate_votes[i + 3].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[3], candidate_votes[i + 3].site - _msf_refGenLength, p_length);
		}

		if (candidate_votes[i + 4].site < _msf_refGenLength)
		{
			get_actuall_genome(t[4], candidate_votes[i + 4].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[4], candidate_votes[i + 4].site - _msf_refGenLength, p_length);
		}

		if (candidate_votes[i + 5].site < _msf_refGenLength)
		{
			get_actuall_genome(t[5], candidate_votes[i + 5].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[5], candidate_votes[i + 5].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 6].site < _msf_refGenLength)
		{
			get_actuall_genome(t[6], candidate_votes[i + 6].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[6], candidate_votes[i + 6].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 7].site < _msf_refGenLength)
		{
			get_actuall_genome(t[7], candidate_votes[i + 7].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[7], candidate_votes[i + 7].site - _msf_refGenLength, p_length);
		}




		BS_Reserve_Banded_BPM_8_SSE(
			t[0], t[1], t[2], t[3],
			t[4], t[5], t[6], t[7],
			p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i + 1].end_site = return_sites[1];
		candidate_votes[i + 2].end_site = return_sites[2];
		candidate_votes[i + 3].end_site = return_sites[3];
		candidate_votes[i + 4].end_site = return_sites[4];
		candidate_votes[i + 5].end_site = return_sites[5];
		candidate_votes[i + 6].end_site = return_sites[6];
		candidate_votes[i + 7].end_site = return_sites[7];

		candidate_votes[i].err = return_sites_error[0];
		candidate_votes[i + 1].err = return_sites_error[1];
		candidate_votes[i + 2].err = return_sites_error[2];
		candidate_votes[i + 3].err = return_sites_error[3];
		candidate_votes[i + 4].err = return_sites_error[4];
		candidate_votes[i + 5].err = return_sites_error[5];
		candidate_votes[i + 6].err = return_sites_error[6];
		candidate_votes[i + 7].err = return_sites_error[7];







		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		if (return_sites_error[0] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		if (return_sites_error[1] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 1].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 1].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 1].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		if (return_sites_error[2] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 2].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 2].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 2].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		if (return_sites_error[3] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 3].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 3].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 3].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;


		tmp_min_err_site = candidate_votes[i + 4].site + candidate_votes[i + 4].end_site;
		if (return_sites_error[4] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 4].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 4].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 4].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 5].site + candidate_votes[i + 5].end_site;
		if (return_sites_error[5] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 5].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 5].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 5].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;


		tmp_min_err_site = candidate_votes[i + 6].site + candidate_votes[i + 6].end_site;
		if (return_sites_error[6] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 6].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 6].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 6].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		tmp_min_err_site = candidate_votes[i + 7].site + candidate_votes[i + 7].end_site;
		if (return_sites_error[7] <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + 7].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + 7].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + 7].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;



		i = i + 8;

	}

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 8);

	int last = i;
	for (; last < *candidate_votes_length; last++)
	{


		if (candidate_votes[last].site < _msf_refGenLength)
		{
			get_actuall_genome(t[last - i], candidate_votes[last].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[last - i], candidate_votes[last].site - _msf_refGenLength, p_length);
		}



	}


	int last_length = *candidate_votes_length - i;

	if (last_length == 1)
	{
		candidate_votes[i].end_site =
			BS_Reserve_Banded_BPM(t[0], p_length, read, t_length, error_threshold, &(candidate_votes[i].err));



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		if (candidate_votes[i].err <= error_threshold
			&&
			pre_err_site != tmp_min_err_site)
		{
			candidate_votes[(*best_mapping_occ)].site = candidate_votes[i].site;
			candidate_votes[(*best_mapping_occ)].err = candidate_votes[i].err;
			candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i].end_site;
			(*best_mapping_occ)++;
		}
		pre_err_site = tmp_min_err_site;

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_8_SSE(
			t[0], t[1], t[2], t[3],
			t[4], t[5], t[6], t[7],
			p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		int inner_i = 0;

		for (inner_i = 0; inner_i < last_length; inner_i++)
		{

			candidate_votes[i + inner_i].end_site = return_sites[inner_i];
			candidate_votes[i + inner_i].err = return_sites_error[inner_i];


			tmp_min_err_site = candidate_votes[i + inner_i].site + candidate_votes[i + inner_i].end_site;
			if (return_sites_error[inner_i] <= error_threshold
				&&
				pre_err_site != tmp_min_err_site)
			{
				candidate_votes[(*best_mapping_occ)].site = candidate_votes[i + inner_i].site;
				candidate_votes[(*best_mapping_occ)].err = candidate_votes[i + inner_i].err;
				candidate_votes[(*best_mapping_occ)].end_site = candidate_votes[i + inner_i].end_site;
				(*best_mapping_occ)++;
			}
			pre_err_site = tmp_min_err_site;


		}

	}

}








inline void map_candidate_votes_mutiple_cut_end_to_end_8(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[8];
	unsigned int return_sites_error[8];

	char* t[8];

	int inner_i = 0;


	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;
	t[2] = t[1] + p_length + 32;
	t[3] = t[2] + p_length + 32;
	t[4] = t[3] + p_length + 32;
	t[5] = t[4] + p_length + 32;
	t[6] = t[5] + p_length + 32;
	t[7] = t[6] + p_length + 32;



	i = 0;

	///这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	(*min_err_index) = -1;
	bitmapper_bs_iter min_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;

	while (i + 8 <= *candidate_votes_length)
	{


		if (candidate_votes[i].site < _msf_refGenLength)
		{
			get_actuall_genome(t[0], candidate_votes[i].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[0], candidate_votes[i].site - _msf_refGenLength, p_length);
		}

		if (candidate_votes[i + 1].site < _msf_refGenLength)
		{
			get_actuall_genome(t[1], candidate_votes[i + 1].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[1], candidate_votes[i + 1].site - _msf_refGenLength, p_length);
		}

		if (candidate_votes[i + 2].site < _msf_refGenLength)
		{
			get_actuall_genome(t[2], candidate_votes[i + 2].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[2], candidate_votes[i + 2].site - _msf_refGenLength, p_length);
		}

		if (candidate_votes[i + 3].site < _msf_refGenLength)
		{
			get_actuall_genome(t[3], candidate_votes[i + 3].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[3], candidate_votes[i + 3].site - _msf_refGenLength, p_length);
		}

		if (candidate_votes[i + 4].site < _msf_refGenLength)
		{
			get_actuall_genome(t[4], candidate_votes[i + 4].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[4], candidate_votes[i + 4].site - _msf_refGenLength, p_length);
		}

		if (candidate_votes[i + 5].site < _msf_refGenLength)
		{
			get_actuall_genome(t[5], candidate_votes[i + 5].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[5], candidate_votes[i + 5].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 6].site < _msf_refGenLength)
		{
			get_actuall_genome(t[6], candidate_votes[i + 6].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[6], candidate_votes[i + 6].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 7].site < _msf_refGenLength)
		{
			get_actuall_genome(t[7], candidate_votes[i + 7].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[7], candidate_votes[i + 7].site - _msf_refGenLength, p_length);
		}

		





		BS_Reserve_Banded_BPM_8_SSE(
			t[0], t[1], t[2], t[3],
			t[4], t[5], t[6], t[7],
			p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);


		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i].err = return_sites_error[0];

		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[0] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*min_err) = return_sites_error[0];
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}


		candidate_votes[i + 1].end_site = return_sites[1];
		candidate_votes[i + 1].err = return_sites_error[1];

		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[1] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}


		candidate_votes[i + 2].end_site = return_sites[2];
		candidate_votes[i + 2].err = return_sites_error[2];

		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[2] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*min_err) = return_sites_error[2];
			(*min_err_index) = i + 2;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 3].end_site = return_sites[3];
		candidate_votes[i + 3].err = return_sites_error[3];

		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[3] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*min_err) = return_sites_error[3];
			(*min_err_index) = i + 3;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 4].end_site = return_sites[4];
		candidate_votes[i + 4].err = return_sites_error[4];

		tmp_min_err_site = candidate_votes[i + 4].site + candidate_votes[i + 4].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[4] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[4]<(*min_err))
		{
			(*min_err) = return_sites_error[4];
			(*min_err_index) = i + 4;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 5].end_site = return_sites[5];
		candidate_votes[i + 5].err = return_sites_error[5];

		tmp_min_err_site = candidate_votes[i + 5].site + candidate_votes[i + 5].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[5] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[5]<(*min_err))
		{
			(*min_err) = return_sites_error[5];
			(*min_err_index) = i + 5;
			min_err_site = tmp_min_err_site;
		}


		candidate_votes[i + 6].end_site = return_sites[6];
		candidate_votes[i + 6].err = return_sites_error[6];

		tmp_min_err_site = candidate_votes[i + 6].site + candidate_votes[i + 6].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[6] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[6]<(*min_err))
		{
			(*min_err) = return_sites_error[6];
			(*min_err_index) = i + 6;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 7].end_site = return_sites[7];
		candidate_votes[i + 7].err = return_sites_error[7];

		tmp_min_err_site = candidate_votes[i + 7].site + candidate_votes[i + 7].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[7] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (return_sites_error[7]<(*min_err))
		{
			(*min_err) = return_sites_error[7];
			(*min_err_index) = i + 7;
			min_err_site = tmp_min_err_site;
		}



		i = i + 8;
	}

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 8);

	int last = i;
	for (; last < *candidate_votes_length; last++)
	{


		if (candidate_votes[last].site < _msf_refGenLength)
		{
			get_actuall_genome(t[last - i], candidate_votes[last].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[last - i], candidate_votes[last].site - _msf_refGenLength, p_length);
		}



	}


	int last_length = *candidate_votes_length - i;

	if (last_length == 1)
	{
		candidate_votes[i].end_site =
			BS_Reserve_Banded_BPM(t[0], p_length, read, t_length, error_threshold, &(candidate_votes[i].err));



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (candidate_votes[i].err == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

	}
	else if (last_length != 0)
	{
		if (last_length <= 4)
		{
			BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
				return_sites, return_sites_error, error_threshold, Peq_SSE);
		}
		else
		{
			BS_Reserve_Banded_BPM_8_SSE(
				t[0], t[1], t[2], t[3],
				t[4], t[5], t[6], t[7],
				p_length, read, t_length,
				return_sites, return_sites_error, error_threshold, Peq_SSE);
		}






		for (inner_i = 0; inner_i < last_length; inner_i++)
		{

			candidate_votes[i + inner_i].end_site = return_sites[inner_i];
			candidate_votes[i + inner_i].err = return_sites_error[inner_i];

			tmp_min_err_site = candidate_votes[i + inner_i].site + candidate_votes[i + inner_i].end_site;
			///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
			if (return_sites_error[inner_i] == (*min_err)
				&&
				min_err_site != tmp_min_err_site)
			{
				(*min_err_index) = -2;
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*min_err) = return_sites_error[inner_i];
				(*min_err_index) = i + inner_i;
				min_err_site = tmp_min_err_site;
			}

		}

	}


	if ((*min_err_index) >= 0)
	{

		candidate_votes[0].err = candidate_votes[(*min_err_index)].err;
		candidate_votes[0].end_site = candidate_votes[(*min_err_index)].end_site;
		candidate_votes[0].site = candidate_votes[(*min_err_index)].site;
	}


}

















inline int map_candidate_votes_mutiple_end_to_end_8(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[8];
	unsigned int return_sites_error[8];

	char* t[8];

	int inner_i = 0;


	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;
	t[2] = t[1] + p_length + 32;
	t[3] = t[2] + p_length + 32;
	t[4] = t[3] + p_length + 32;
	t[5] = t[4] + p_length + 32;
	t[6] = t[5] + p_length + 32;
	t[7] = t[6] + p_length + 32;



	i = 0;

	//这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	(*min_err_index) = -1;

	bitmapper_bs_iter min_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;

	while (i + 8 <= *candidate_votes_length)
	{

		if (candidate_votes[i].site < _msf_refGenLength)
		{
			get_actuall_genome(t[0], candidate_votes[i].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[0], candidate_votes[i].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 1].site < _msf_refGenLength)
		{
			get_actuall_genome(t[1], candidate_votes[i + 1].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[1], candidate_votes[i + 1].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 2].site < _msf_refGenLength)
		{
			get_actuall_genome(t[2], candidate_votes[i + 2].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[2], candidate_votes[i + 2].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 3].site < _msf_refGenLength)
		{
			get_actuall_genome(t[3], candidate_votes[i + 3].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[3], candidate_votes[i + 3].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 4].site < _msf_refGenLength)
		{
			get_actuall_genome(t[4], candidate_votes[i + 4].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[4], candidate_votes[i + 4].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 5].site < _msf_refGenLength)
		{
			get_actuall_genome(t[5], candidate_votes[i + 5].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[5], candidate_votes[i + 5].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 6].site < _msf_refGenLength)
		{
			get_actuall_genome(t[6], candidate_votes[i + 6].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[6], candidate_votes[i + 6].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 7].site < _msf_refGenLength)
		{
			get_actuall_genome(t[7], candidate_votes[i + 7].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[7], candidate_votes[i + 7].site - _msf_refGenLength, p_length);
		}


		

		BS_Reserve_Banded_BPM_8_SSE(
			t[0], t[1], t[2], t[3],
			t[4], t[5], t[6], t[7],
			p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);




		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i].err = return_sites_error[0];


		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[0] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*min_err) = return_sites_error[0];
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}


		candidate_votes[i + 1].end_site = return_sites[1];
		candidate_votes[i + 1].err = return_sites_error[1];


		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[1] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}

		candidate_votes[i + 2].end_site = return_sites[2];
		candidate_votes[i + 2].err = return_sites_error[2];


		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[2] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*min_err) = return_sites_error[2];
			(*min_err_index) = i + 2;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 3].end_site = return_sites[3];
		candidate_votes[i + 3].err = return_sites_error[3];


		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[3] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*min_err) = return_sites_error[3];
			(*min_err_index) = i + 3;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 4].end_site = return_sites[4];
		candidate_votes[i + 4].err = return_sites_error[4];


		tmp_min_err_site = candidate_votes[i + 4].site + candidate_votes[i + 4].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[4] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[4]<(*min_err))
		{
			(*min_err) = return_sites_error[4];
			(*min_err_index) = i + 4;
			min_err_site = tmp_min_err_site;
		}


		candidate_votes[i + 5].end_site = return_sites[5];
		candidate_votes[i + 5].err = return_sites_error[5];


		tmp_min_err_site = candidate_votes[i + 5].site + candidate_votes[i + 5].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[5] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[5]<(*min_err))
		{
			(*min_err) = return_sites_error[5];
			(*min_err_index) = i + 5;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 6].end_site = return_sites[6];
		candidate_votes[i + 6].err = return_sites_error[6];


		tmp_min_err_site = candidate_votes[i + 6].site + candidate_votes[i + 6].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[6] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[6]<(*min_err))
		{
			(*min_err) = return_sites_error[6];
			(*min_err_index) = i + 6;
			min_err_site = tmp_min_err_site;
		}



		candidate_votes[i + 7].end_site = return_sites[7];
		candidate_votes[i + 7].err = return_sites_error[7];


		tmp_min_err_site = candidate_votes[i + 7].site + candidate_votes[i + 7].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[7] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[7]<(*min_err))
		{
			(*min_err) = return_sites_error[7];
			(*min_err_index) = i + 7;
			min_err_site = tmp_min_err_site;
		}

		i = i + 8;
	}


	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 8);


	int last = i;
	for (; last < *candidate_votes_length; last++)
	{
		if (candidate_votes[last].site < _msf_refGenLength)
		{
			get_actuall_genome(t[last - i], candidate_votes[last].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[last - i], candidate_votes[last].site - _msf_refGenLength, p_length);
		}
	}


	int last_length = *candidate_votes_length - i;

	if (last_length == 1)
	{
		candidate_votes[i].end_site =
			BS_Reserve_Banded_BPM(t[0], p_length, read, t_length, error_threshold, &(candidate_votes[i].err));



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (candidate_votes[i].err == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}



	}
	else if (last_length != 0)
	{
		if (last_length <= 4)
		{
			BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
				return_sites, return_sites_error, error_threshold, Peq_SSE);
		}
		else
		{
			BS_Reserve_Banded_BPM_8_SSE(
				t[0], t[1], t[2], t[3],
				t[4], t[5], t[6], t[7],
				p_length, read, t_length,
				return_sites, return_sites_error, error_threshold, Peq_SSE);
		}




		for (inner_i = 0; inner_i < last_length; inner_i++)
		{

			candidate_votes[i + inner_i].end_site = return_sites[inner_i];
			candidate_votes[i + inner_i].err = return_sites_error[inner_i];


			tmp_min_err_site = candidate_votes[i + inner_i].site + candidate_votes[i + inner_i].end_site;
			///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
			if (return_sites_error[inner_i] == (*min_err)
				&&
				min_err_site != tmp_min_err_site)
			{
				(*min_err_index) = -2;

				if ((*min_err) == 0)
				{
					return 0;
				}
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*min_err) = return_sites_error[inner_i];
				(*min_err_index) = i + inner_i;
				min_err_site = tmp_min_err_site;
			}


		}

	}



	if ((*min_err_index) >= 0)
	{

		candidate_votes[0].err = candidate_votes[(*min_err_index)].err;
		candidate_votes[0].end_site = candidate_votes[(*min_err_index)].end_site;
		candidate_votes[0].site = candidate_votes[(*min_err_index)].site;
	}





	return 1;

}










inline int map_candidate_votes_mutiple_end_to_end_4(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[4];
	unsigned int return_sites_error[4];
	int mumber_of_exact_matches = 0;

	char* t[4];

	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;
	t[2] = t[1] + p_length + 32;
	t[3] = t[2] + p_length + 32;


	i = 0;

	//这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	(*min_err_index) = -1;

	bitmapper_bs_iter min_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;

	while (i + 4 <= *candidate_votes_length)
	{

		if (candidate_votes[i].site < _msf_refGenLength)
		{
			get_actuall_genome(t[0], candidate_votes[i].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[0], candidate_votes[i].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 1].site < _msf_refGenLength)
		{
			get_actuall_genome(t[1], candidate_votes[i + 1].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[1], candidate_votes[i + 1].site - _msf_refGenLength, p_length);
		}


		if (candidate_votes[i + 2].site < _msf_refGenLength)
		{
			get_actuall_genome(t[2], candidate_votes[i + 2].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[2], candidate_votes[i + 2].site - _msf_refGenLength, p_length);
		}



		if (candidate_votes[i + 3].site < _msf_refGenLength)
		{
			get_actuall_genome(t[3], candidate_votes[i + 3].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[3], candidate_votes[i + 3].site - _msf_refGenLength, p_length);
		}



		BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);


		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i + 1].end_site = return_sites[1];
		candidate_votes[i + 2].end_site = return_sites[2];
		candidate_votes[i + 3].end_site = return_sites[3];

		candidate_votes[i].err = return_sites_error[0];
		candidate_votes[i + 1].err = return_sites_error[1];
		candidate_votes[i + 2].err = return_sites_error[2];
		candidate_votes[i + 3].err = return_sites_error[3];



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[0] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*min_err) = return_sites_error[0];
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

		/**
		if (return_sites_error[0] == 0)
		{
		mumber_of_exact_matches++;
		}
		**/





		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[1] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}
		/**
		if (return_sites_error[1] == 0)
		{
		mumber_of_exact_matches++;
		}
		**/





		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[2] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*min_err) = return_sites_error[2];
			(*min_err_index) = i + 2;
			min_err_site = tmp_min_err_site;
		}
		/**
		if (return_sites_error[2] == 0)
		{
		mumber_of_exact_matches++;
		}
		**/







		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[3] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*min_err) = return_sites_error[3];
			(*min_err_index) = i + 3;
			min_err_site = tmp_min_err_site;
		}
		/**
		if (return_sites_error[3] == 0)
		{
		mumber_of_exact_matches++;
		}
		**/




		i = i + 4;

		/**
		///这是不符合唯一匹配的要求
		if (mumber_of_exact_matches >= 2)
		{
		(*min_err_index) = -2;
		return 0;
		}
		**/


	}


	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	int last = i;
	for (; last < *candidate_votes_length; last++)
	{
		if (candidate_votes[last].site < _msf_refGenLength)
		{
			get_actuall_genome(t[last - i], candidate_votes[last].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(t[last - i], candidate_votes[last].site - _msf_refGenLength, p_length);
		}
	}


	int last_length = *candidate_votes_length - i;

	if (last_length == 1)
	{
		candidate_votes[i].end_site =
			BS_Reserve_Banded_BPM(t[0], p_length, read, t_length, error_threshold, &(candidate_votes[i].err));



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (candidate_votes[i].err == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*min_err_index) = -2;

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

		/**
		if (candidate_votes[i].err == 0)
		{
		mumber_of_exact_matches++;
		}
		**/

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_4_SSE(t[0], t[1], t[2], t[3], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		int inner_i = 0;

		for (inner_i = 0; inner_i < last_length; inner_i++)
		{

			candidate_votes[i + inner_i].end_site = return_sites[inner_i];
			candidate_votes[i + inner_i].err = return_sites_error[inner_i];


			tmp_min_err_site = candidate_votes[i + inner_i].site + candidate_votes[i + inner_i].end_site;
			///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
			if (return_sites_error[inner_i] == (*min_err)
				&&
				min_err_site != tmp_min_err_site)
			{
				(*min_err_index) = -2;

				if ((*min_err) == 0)
				{
					return 0;
				}
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*min_err) = return_sites_error[inner_i];
				(*min_err_index) = i + inner_i;
				min_err_site = tmp_min_err_site;
			}
			/**
			if (return_sites_error[inner_i] == 0)
			{
			mumber_of_exact_matches++;
			}
			**/

		}

	}



	if ((*min_err_index) >= 0)
	{

		candidate_votes[0].err = candidate_votes[(*min_err_index)].err;
		candidate_votes[0].end_site = candidate_votes[(*min_err_index)].end_site;
		candidate_votes[0].site = candidate_votes[(*min_err_index)].site;
	}


	/**
	if (mumber_of_exact_matches >= 2)
	{
	(*min_err_index) = -2;
	return 0;
	}
	**/




	return 1;

}


















inline void output_sam_end_to_end_return(
	bitmapper_bs_iter site, 
	bitmapper_bs_iter end_site, 
	bitmapper_bs_iter start_site,
	int read_length,
	map_result* result
	)
{




	bitmapper_bs_iter map_location = site;

	int flag;

	if (map_location >= _msf_refGenLength)
	{

		map_location = map_location + end_site;

		map_location = _msf_refGenLength * 2 - map_location - 1;

		flag = 16;

	}
	else
	{
		map_location = map_location + start_site;

		flag = 0;
	}


	int now_ref_name = 0;


	for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
	{
		if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
			break;
	}


	map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;

	
	result->flag = flag;
	result->chrome_id = now_ref_name;
	result->site = map_location;



}




inline void directly_output_read1_return_buffer(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length,
	Output_buffer_sub_block* sub_block)
{



	if (name[0] == '@')
	{
		output_to_buffer_char_no_length(sub_block, name + 1);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{
		output_to_buffer_char_no_length(sub_block, name);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}



	///代表read比对上正向链了
	if (result->flag == 0)
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ1;
	}
	else ///代表read比对上反向互补链
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ1;
	}


	output_to_buffer_int(sub_block, result->flag);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	output_to_buffer_char_no_length(sub_block, _ih_refGenName[result->chrome_id]._rg_chrome_name);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	output_to_buffer_int(sub_block, result->site);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '2';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;



	output_to_buffer_char_no_length(sub_block, result->cigar);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;



	

	if (another_result->site > result->site)
	{
		/**
		result_string << "=\t" << another_result->site << "\t"
			<< another_result->site - result->site + another_matched_length << "\t";
		**/

		sub_block->buffer[sub_block->length] = '=';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_int(sub_block, another_result->site);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;


		output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

	}
	else if (another_result->site < result->site)
	{
		/**
		result_string << "=\t" << another_result->site << "\t-"
			<< result->site - another_result->site + matched_length << "\t";
		**/

		sub_block->buffer[sub_block->length] = '=';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_int(sub_block, another_result->site);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '-';
		sub_block->length++;

		output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{
		if (matched_length >= another_matched_length)
		{
			/**
			result_string << "=\t" << another_result->site << "\t"
				<< matched_length << "\t";
			**/

			sub_block->buffer[sub_block->length] = '=';
			sub_block->length++;

			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

			output_to_buffer_int(sub_block, another_result->site);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

			output_to_buffer_int(sub_block, matched_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			/**
			result_string << "=\t" << another_result->site << "\t"
				<< another_matched_length << "\t";
			**/

			sub_block->buffer[sub_block->length] = '=';
			sub_block->length++;

			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

			output_to_buffer_int(sub_block, another_result->site);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;


			output_to_buffer_int(sub_block, another_matched_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

		}

	}





	if (result->flag & RC_MATE_MAPPED)
	{
		/**
		result_string << read << "\t";

		result_string << qulity << "\t";
		**/


		output_to_buffer_char_length(sub_block, read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{
		

		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}

		/**
		result_string << r_read << "\t";
		result_string << qulity << "\t";
		**/

		output_to_buffer_char_length(sub_block, r_read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

	}


	///result_string << "NM:i:" << result->err << "\n";

	sub_block->buffer[sub_block->length] = 'N';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'M';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'i';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	output_to_buffer_int(sub_block, result->err);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\n';
	sub_block->length++;

	/**
	fprintf(stderr, "********\nsub_block->length: %llu\n", sub_block->length);
	fflush(stderr);
	fprintf(stderr, "sub_block->size: %llu\n", sub_block->size);
	fflush(stderr);
	**/

}




inline void directly_output_read1(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length)
{


	if (name[0] == '@')
	{
		fprintf(schema_out_fp, "%s\t", name + 1);
	}
	else
	{
		fprintf(schema_out_fp, "%s\t", name);
	}


	///代表read比对上正向链了
	if (result->flag == 0)
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ1;
	}
	else ///代表read比对上反向互补链
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ1;
	}

	fprintf(schema_out_fp, "%d\t", result->flag);

	fprintf(schema_out_fp, "%s\t", _ih_refGenName[result->chrome_id]._rg_chrome_name);


	fprintf(schema_out_fp, "%llu\t", result->site);

	fprintf(schema_out_fp, "255\t");

	fprintf(schema_out_fp, "%s\t", result->cigar);

	if (another_result->site > result->site)
	{
		fprintf(schema_out_fp, "=\t%llu\t%llu\t", 
			another_result->site, another_result->site - result->site + another_matched_length);
	}
	else if (another_result->site < result->site)
	{
		fprintf(schema_out_fp, "=\t%llu\t-%llu\t", 
			another_result->site, result->site - another_result->site + matched_length);


	}
	else
	{
		if (matched_length >= another_matched_length)
		{
			fprintf(schema_out_fp, "=\t%llu\t%llu\t",
				another_result->site, matched_length);
		}
		else
		{
			fprintf(schema_out_fp, "=\t%llu\t%llu\t",
				another_result->site, another_matched_length);
		}
		
	}
	

	



	///if (result->flag == 0)
	if (result->flag & RC_MATE_MAPPED)
	{
		fprintf(schema_out_fp, "%s\t", read);

		fprintf(schema_out_fp, "%s\t", qulity);
	}
	else
	{
		fprintf(schema_out_fp, "%s\t", r_read);

		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}

		fprintf(schema_out_fp, "%s\t", qulity);

	}








	fprintf(schema_out_fp, "NM:i:%d\n", result->err);




}









inline void directly_output_read2_return_buffer(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length,
	Output_buffer_sub_block* sub_block)
{



	///这个和正向read是相反的
	if (result->flag == 0)
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ2;
	}
	else
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ2;
	}



	if (name[0] == '@')
	{
		output_to_buffer_char_no_length(sub_block, name + 1);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{
		output_to_buffer_char_no_length(sub_block, name);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}



	output_to_buffer_int(sub_block, result->flag);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;



	output_to_buffer_char_no_length(sub_block, _ih_refGenName[result->chrome_id]._rg_chrome_name);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	output_to_buffer_int(sub_block, result->site);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '2';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;



	output_to_buffer_char_no_length(sub_block, result->cigar);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;



	if (another_result->site > result->site)
	{
		/**
		result_string << "=\t" << another_result->site << "\t"
			<< another_result->site - result->site + another_matched_length << "\t";
		**/


		sub_block->buffer[sub_block->length] = '=';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_int(sub_block, another_result->site);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;


		output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;


	}
	else if (another_result->site < result->site)
	{
		/**
		result_string << "=\t" << another_result->site << "\t-"
			<< result->site - another_result->site + matched_length << "\t";
		**/


		sub_block->buffer[sub_block->length] = '=';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_int(sub_block, another_result->site);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '-';
		sub_block->length++;

		output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

	}
	else
	{
		if (matched_length >= another_matched_length)
		{
			/**
			result_string << "=\t" << another_result->site << "\t-"
				<< matched_length << "\t";
			**/


			sub_block->buffer[sub_block->length] = '=';
			sub_block->length++;

			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

			output_to_buffer_int(sub_block, another_result->site);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '-';
			sub_block->length++;

			output_to_buffer_int(sub_block, matched_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

		}
		else
		{
			/**
			result_string << "=\t" << another_result->site << "\t-"
				<< another_matched_length << "\t";
			**/


			sub_block->buffer[sub_block->length] = '=';
			sub_block->length++;

			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

			output_to_buffer_int(sub_block, another_result->site);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '-';
			sub_block->length++;


			output_to_buffer_int(sub_block, another_matched_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}

	}



	if (result->flag & RC_MAPPED)
	{

		

		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}
		
		/**
		result_string << read << "\t";
		result_string << qulity << "\t";
		**/

		output_to_buffer_char_length(sub_block, read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

	}
	else
	{
		/**
		result_string << r_read << "\t";
		result_string << qulity << "\t";
		**/

		output_to_buffer_char_length(sub_block, r_read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

	}




	///result_string << "NM:i:" << result->err << "\n";


	sub_block->buffer[sub_block->length] = 'N';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'M';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'i';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	output_to_buffer_int(sub_block, result->err);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\n';
	sub_block->length++;


}








inline void directly_output_read2(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length)
{



	///result->flag == 0代表read2的反向互补链比对到了参考组的正向上
	//所以在这种情况下，read2本身实际上比到了参考组反向互补链上
	if (result->flag == 0)
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ2;
	}
	else 
	{
		result->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ2;
	}
	///result->flag == 16代表read2的反向互补链比对到了参考组反向互补链上
	//所以在这种情况下，read2本身实际上比到了参考组正向链上


	if (name[0] == '@')
	{
		fprintf(schema_out_fp, "%s\t", name + 1);
	}
	else
	{
		fprintf(schema_out_fp, "%s\t", name);
	}


	fprintf(schema_out_fp, "%d\t", result->flag);

	fprintf(schema_out_fp, "%s\t", _ih_refGenName[result->chrome_id]._rg_chrome_name);


	fprintf(schema_out_fp, "%llu\t", result->site);

	fprintf(schema_out_fp, "255\t");

	fprintf(schema_out_fp, "%s\t", result->cigar);

	///fprintf(schema_out_fp, "*\t0\t0\t");
	///fprintf(schema_out_fp, "=\t%llu\t%lld\t", another_result->site, (long long)(another_result->site) - (long long)(result->site));




	if (another_result->site > result->site)
	{
		fprintf(schema_out_fp, "=\t%llu\t%llu\t",
			another_result->site, another_result->site - result->site + another_matched_length);
	}
	else if (another_result->site < result->site)
	{
		fprintf(schema_out_fp, "=\t%llu\t-%llu\t",
			another_result->site, result->site - another_result->site + matched_length);
	}
	else
	{
		if (matched_length >= another_matched_length)
		{
			fprintf(schema_out_fp, "=\t%llu\t-%llu\t",
				another_result->site, matched_length);
		}
		else
		{
			fprintf(schema_out_fp, "=\t%llu\t-%llu\t",
				another_result->site, another_matched_length);
		}

	}

	///这里好像read存的是反向互补链
	///相对应的，r_read存的是正向链
	///但是qulity存的是正向的质量值
	///注意read2里所有的东西都是反着来的
	///所以read比对到反向互补链的话，就输出反向互补序列
	if (result->flag & RC_MAPPED)
	{

		fprintf(schema_out_fp, "%s\t", read);

		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}

		fprintf(schema_out_fp, "%s\t", qulity);


	}
	else
	{
		fprintf(schema_out_fp, "%s\t", r_read);

		fprintf(schema_out_fp, "%s\t", qulity);
	}








	fprintf(schema_out_fp, "NM:i:%d\n", result->err);




}






inline void output_sam_end_to_end(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length)
{




	bitmapper_bs_iter map_location = site;

	int flag;

	if (map_location >= _msf_refGenLength)
	{

		map_location = map_location + end_site;

		map_location = _msf_refGenLength * 2 - map_location - 1;

		flag = 16;

		///fprintf(stderr, "-\n");

	}
	else
	{
		map_location = map_location + start_site;

		flag = 0;

		///fprintf(stderr, "+\n");
	}


	int now_ref_name = 0;


	for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
	{
		if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
			break;
	}


	map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;


	if (name[0] == '@')
	{
		fprintf(schema_out_fp, "%s\t", name + 1);
	}
	else
	{
		fprintf(schema_out_fp, "%s\t", name);
	}


	fprintf(schema_out_fp, "%d\t", flag);

	fprintf(schema_out_fp, "%s\t", _ih_refGenName[now_ref_name]._rg_chrome_name);


	fprintf(schema_out_fp, "%llu\t", map_location);

	fprintf(schema_out_fp, "255\t");

	fprintf(schema_out_fp, "%s\t", best_cigar);

	fprintf(schema_out_fp, "*\t0\t0\t");



	if (flag == 0)
	{
		fprintf(schema_out_fp, "%s\t", read);

		fprintf(schema_out_fp, "%s\t", qulity);
	}
	else
	{
		fprintf(schema_out_fp, "%s\t", r_read);

		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}

		fprintf(schema_out_fp, "%s\t", qulity);

	}





	///fprintf(schema_out_fp, "end_site: %lld\t", end_site);


	fprintf(schema_out_fp, "NM:i:%d\n", err);


	/**
	fprintf(stderr, "%d\n", flag);

	fprintf(stderr, "%s\n", _ih_refGenName[now_ref_name]._rg_chrome_name);

	fprintf(stderr, "%llu\n", map_location);

	fprintf(stderr, "err: %u\n", err);
	**/



}






inline void output_sam_end_to_end_output_buffer(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length, Output_buffer_sub_block* sub_block)
{



	bitmapper_bs_iter map_location = site;

	int flag;

	if (map_location >= _msf_refGenLength)
	{

		map_location = map_location + end_site;

		map_location = _msf_refGenLength * 2 - map_location - 1;

		flag = 16;

	}
	else
	{
		map_location = map_location + start_site;

		flag = 0;
	}


	int now_ref_name = 0;


	for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
	{
		if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
			break;
	}


	map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;

	/**
	if (name[0] == '@')
	{
		///fprintf(schema_out_fp, "%s\t", name + 1);
		result_string << name + 1 << "\t";
	}
	else
	{
		///fprintf(schema_out_fp, "%s\t", name);
		result_string << name << "\t";
	}
	**/


	if (name[0] == '@')
	{
		output_to_buffer_char_no_length(sub_block, name + 1);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{
		output_to_buffer_char_no_length(sub_block, name);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}






	///result_string << flag << "\t";
	output_to_buffer_int(sub_block, flag);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;

	///result_string << _ih_refGenName[now_ref_name]._rg_chrome_name << "\t";
	output_to_buffer_char_no_length(sub_block, _ih_refGenName[now_ref_name]._rg_chrome_name);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	///result_string << map_location << "\t";
	output_to_buffer_int(sub_block, map_location);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	///result_string << "255\t";
	sub_block->buffer[sub_block->length] = '2';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	///result_string << best_cigar << "\t";
	output_to_buffer_char_no_length(sub_block, best_cigar);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	///result_string << "*\t0\t0\t";
	sub_block->buffer[sub_block->length] = '*';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '0';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '0';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;



	if (flag == 0)
	{
		/**
		fprintf(schema_out_fp, "%s\t", read);
		fprintf(schema_out_fp, "%s\t", qulity);
		**/

		
		


		///result_string << read << "\t";
		output_to_buffer_char_length(sub_block, read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		///result_string << qulity << "\t";
		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{


		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}

		

		///result_string << r_read << "\t";
		output_to_buffer_char_length(sub_block, r_read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		///result_string << qulity << "\t";
		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

	}



	///fprintf(schema_out_fp, "NM:i:%d\n", err);
	///result_string << "NM:i:" << err << "\n";
	sub_block->buffer[sub_block->length] = 'N';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'M';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'i';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	output_to_buffer_int(sub_block, err);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\n';
	sub_block->length++;





	///fprintf(stderr, "%s", result_string.str().c_str());

}













inline void output_sam_end_to_end_pbat_output_buffer(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length, Output_buffer_sub_block* sub_block)
{



	bitmapper_bs_iter map_location = site;

	int flag;

	if (map_location >= _msf_refGenLength)
	{

		map_location = map_location + end_site;

		map_location = _msf_refGenLength * 2 - map_location - 1;

		flag = 0;

	}
	else
	{
		map_location = map_location + start_site;

		flag = 16;
	}


	int now_ref_name = 0;


	for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
	{
		if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
			break;
	}


	map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;

	/**
	if (name[0] == '@')
	{
		result_string << name + 1 << "\t";
	}
	else
	{
		result_string << name << "\t";
	}
	**/
	if (name[0] == '@')
	{
		output_to_buffer_char_no_length(sub_block, name + 1);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{
		output_to_buffer_char_no_length(sub_block, name);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}


	/**
	result_string << flag << "\t";
	result_string << _ih_refGenName[now_ref_name]._rg_chrome_name << "\t";
	result_string << map_location << "\t";
	result_string << "255\t";
	result_string << best_cigar << "\t";
	result_string << "*\t0\t0\t";
	**/


	output_to_buffer_int(sub_block, flag);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;

	output_to_buffer_char_no_length(sub_block, _ih_refGenName[now_ref_name]._rg_chrome_name);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	output_to_buffer_int(sub_block, map_location);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	sub_block->buffer[sub_block->length] = '2';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '5';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	output_to_buffer_char_no_length(sub_block, best_cigar);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;


	sub_block->buffer[sub_block->length] = '*';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '0';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '0';
	sub_block->length++;
	sub_block->buffer[sub_block->length] = '\t';
	sub_block->length++;




	if (flag == 16)
	{

		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}

		/**
		result_string << read << "\t";
		result_string << qulity << "\t";
		**/


		output_to_buffer_char_length(sub_block, read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;
	}
	else
	{
		/**
		result_string << r_read << "\t";
		result_string << qulity << "\t";
		**/


		///result_string << r_read << "\t";
		output_to_buffer_char_length(sub_block, r_read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		///result_string << qulity << "\t";
		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

	}


	///result_string << "NM:i:" << err << "\n";


	sub_block->buffer[sub_block->length] = 'N';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'M';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = 'i';
	sub_block->length++;

	sub_block->buffer[sub_block->length] = ':';
	sub_block->length++;

	output_to_buffer_int(sub_block, err);
	//上面扩容时有32的余量，所以下面\t不需要检查了
	sub_block->buffer[sub_block->length] = '\n';
	sub_block->length++;


}










inline void output_sam_end_to_end_pbat(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length)
{


	///fprintf(stderr, "******************\n");

	bitmapper_bs_iter map_location = site;

	int flag;

	if (map_location >= _msf_refGenLength)
	{

		map_location = map_location + end_site;

		map_location = _msf_refGenLength * 2 - map_location - 1;

		///flag = 16;
		flag = 0;

	}
	else
	{
		map_location = map_location + start_site;

		///flag = 0;
		flag = 16;
	}


	int now_ref_name = 0;


	for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
	{
		if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
			break;
	}


	map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;


	if (name[0] == '@')
	{
		fprintf(schema_out_fp, "%s\t", name + 1);
	}
	else
	{
		fprintf(schema_out_fp, "%s\t", name);
	}


	fprintf(schema_out_fp, "%d\t", flag);

	fprintf(schema_out_fp, "%s\t", _ih_refGenName[now_ref_name]._rg_chrome_name);


	fprintf(schema_out_fp, "%llu\t", map_location);

	fprintf(schema_out_fp, "255\t");

	fprintf(schema_out_fp, "%s\t", best_cigar);

	fprintf(schema_out_fp, "*\t0\t0\t");



	if (flag == 16)
	{
		fprintf(schema_out_fp, "%s\t", read);



		char tmp_c;
		int end_i = read_length / 2;

		for (int i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[read_length - i - 1];
			qulity[read_length - i - 1] = tmp_c;
		}

		fprintf(schema_out_fp, "%s\t", qulity);


		
	}
	else
	{
		fprintf(schema_out_fp, "%s\t", r_read);

		fprintf(schema_out_fp, "%s\t", qulity);

	}








	fprintf(schema_out_fp, "NM:i:%d\n", err);




}












inline void calculate_cigar_end_to_end(
	bitmapper_bs_iter site,
	char* path, int path_length,
	char* cigar,
	int t_length)
{

	int i = 0;
	char pre_ciga = 100;
	int pre_ciga_length = 0;




	cigar[0] = '\0';

	pre_ciga = 100;
	pre_ciga_length = 0;

	///反向互补链比对上
	if (site >= _msf_refGenLength)
	{

		for (i = 0; i < path_length; i++)
		{

			if (pre_ciga != map_cigar[path[i]])
			{
				if (pre_ciga_length != 0)
				{
					sprintf(cigar + strlen(cigar), "%d", pre_ciga_length);
					sprintf(cigar + strlen(cigar), "%c", pre_ciga);
				}

				pre_ciga = map_cigar[path[i]];
				pre_ciga_length = 1;

			}
			else
			{
				pre_ciga_length++;
			}
		}

		if (pre_ciga_length != 0)
		{
			sprintf(cigar + strlen(cigar), "%d", pre_ciga_length);
			sprintf(cigar + strlen(cigar), "%c", pre_ciga);
		}

	}
	else
	{
		for (i = path_length - 1; i >= 0; i--)
		{

			if (pre_ciga != map_cigar[path[i]])
			{
				if (pre_ciga_length != 0)
				{
					sprintf(cigar + strlen(cigar), "%d", pre_ciga_length);
					sprintf(cigar + strlen(cigar), "%c", pre_ciga);
				}

				pre_ciga = map_cigar[path[i]];
				pre_ciga_length = 1;

			}
			else
			{
				pre_ciga_length++;
			}
		}

		if (pre_ciga_length != 0)
		{
			sprintf(cigar + strlen(cigar), "%d", pre_ciga_length);
			sprintf(cigar + strlen(cigar), "%c", pre_ciga);
		}
	}


	


}







inline int calculate_best_map_cigar_end_to_end_return(bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar, map_result* result, int* matched_length)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;



	int start_site;
	int path_length;





	if (result->err != 0)
	{


		if (result->origin_site < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, result->origin_site, p_length);

		}
		else
		{
			get_actuall_rc_genome(tmp_ref, result->origin_site - _msf_refGenLength, p_length);

		}



		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(result->err),
			cigar, result->end_site, &start_site, path, &path_length, matrix_bit) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{



			///double start = Get_T();
			calculate_cigar_end_to_end(result->origin_site, path, path_length, best_cigar, t_length);
			///total = total + Get_T() - start;
		}


	}
	else
	{

		///start_site = 0;
		start_site = result->end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}



	///result->err = candidate_votes[0].err;
	output_sam_end_to_end_return(
		result->origin_site,
		result->end_site,
		start_site,
		read_length,
		result
		);



	(*matched_length) = result->end_site - start_site + 1;


	return 0;




}






inline int calculate_best_map_cigar_end_to_end_output_buffer(seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar, Output_buffer_sub_block* sub_block)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;



	int start_site;
	int path_length;





	if (candidate_votes[0].err != 0)
	{


		if (candidate_votes[0].site < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, candidate_votes[0].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(tmp_ref, candidate_votes[0].site - _msf_refGenLength, p_length);

		}






		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length);
		}

	}
	else
	{

		///start_site = 0;
		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}

	output_sam_end_to_end_output_buffer
		(name, read, r_read, qulity,
		candidate_votes[0].site,
		candidate_votes[0].end_site,
		start_site,
		candidate_votes[0].err, best_cigar,
		read_length, sub_block);
	/**
	output_sam_end_to_end_result_string(name, read, r_read, qulity,
		candidate_votes[0].site,
		candidate_votes[0].end_site,
		start_site,
		candidate_votes[0].err, best_cigar,
		read_length, result_string);
	**/



	return 0;




}







inline int calculate_best_map_cigar_end_to_end_pbat_output_buffer(seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar, Output_buffer_sub_block* sub_block)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;



	int start_site;
	int path_length;





	if (candidate_votes[0].err != 0)
	{


		if (candidate_votes[0].site < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, candidate_votes[0].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(tmp_ref, candidate_votes[0].site - _msf_refGenLength, p_length);

		}






		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length);
		}

	}
	else
	{

		///start_site = 0;
		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}



	output_sam_end_to_end_pbat_output_buffer(name, read, r_read, qulity,
		candidate_votes[0].site,
		candidate_votes[0].end_site,
		start_site,
		candidate_votes[0].err, best_cigar,
		read_length, sub_block);



	return 0;




}







inline int calculate_best_map_cigar_end_to_end(seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	


	int start_site;
	int path_length;





	if (candidate_votes[0].err != 0)
	{


		if (candidate_votes[0].site < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, candidate_votes[0].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(tmp_ref, candidate_votes[0].site - _msf_refGenLength, p_length);

		}






		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			///double start = Get_T();
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length);
			///total = total + Get_T() - start;
		}


	}
	else
	{

		///start_site = 0;
		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}


	
	output_sam_end_to_end(name, read, r_read, qulity,
		candidate_votes[0].site,
		candidate_votes[0].end_site,
		start_site,
		candidate_votes[0].err, best_cigar,
		read_length);
	


	return 0;




}






inline int calculate_best_map_cigar_end_to_end_pbat(seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;



	int start_site;
	int path_length;





	if (candidate_votes[0].err != 0)
	{


		if (candidate_votes[0].site < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, candidate_votes[0].site, p_length);

		}
		else
		{
			get_actuall_rc_genome(tmp_ref, candidate_votes[0].site - _msf_refGenLength, p_length);

		}




		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			///double start = Get_T();
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length);
			///total = total + Get_T() - start;
		}


		/**
		if (strcmp(name, "@SRR771401.1391") == 0)
		{
		fprintf(stderr, "read   : %s\n", read);
		fprintf(stderr, "tmp_ref: %s\n", tmp_ref);
		fprintf(stderr, "best_cigar: %s\n", best_cigar);
		fprintf(stderr, "start_site: %llu\n", start_site);
		fprintf(stderr, "end_site: %llu\n", candidate_votes[0].end_site);
		fprintf(stderr, "site: %llu\n", candidate_votes[0].site);
		fprintf(stderr, "read_length: %llu\n", read_length);
		fprintf(stderr, "err: %llu\n", candidate_votes[0].err);
		for (int iik = path_length - 1; iik >= 0; iik--)
		{
		fprintf(stderr, "%u", path[iik]);
		}

		fprintf(stderr, "\n");

		}
		**/




		/**
		if (strcmp(read, "AGGTTAGTGTTGAAGAGAGAAATAAAAAAAAAATAGTTTGTATAAAAATATTTGATGATATTTTGAGAAATTAAAATATGGTATAGT") == 0)
		{
		fprintf(stderr, "name: %s\n", name);
		fprintf(stderr, "read: %s\n", read);
		fprintf(stderr, "tmp_ref: %s\n", tmp_ref);
		fprintf(stderr, "best_cigar: %s\n", best_cigar);
		fprintf(stderr, "start_site: %d\n", start_site);
		fprintf(stderr, "end_site: %d\n", candidate_votes[0].end_site);

		for (int iik = path_length - 1; iik >= 0; iik--)
		{
		fprintf(stderr, "%u", path[iik]);
		}

		fprintf(stderr, "\n");

		}
		**/



	}
	else
	{

		///start_site = 0;
		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}



	output_sam_end_to_end_pbat(name, read, r_read, qulity,
		candidate_votes[0].site,
		candidate_votes[0].end_site,
		start_site,
		candidate_votes[0].err, best_cigar,
		read_length);



	return 0;




}








inline int try_process_unique_mismatch_end_to_end(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, int* get_error)
{
	bitmapper_bs_iter tmp_SA_length = 0;


	locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, *matched_length, 0);



	int i = 0;
	int error = 0;


	if (*matched_length > first_C_site)
	{
		*matched_length = first_C_site;
	}


	if (*matched_length != read_length)
	{
		int need_read_length = read_length - *matched_length;

		if (locates[0] < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, locates[0] + *matched_length, need_read_length);


		}
		else
		{

			get_actuall_rc_genome(tmp_ref, locates[0] + *matched_length - _msf_refGenLength,
				need_read_length);

		}

		int read_i = *matched_length;

		for (i = 0; i < need_read_length; i++)
		{
			if (read[read_i] != tmp_ref[i])
			{
				if (!(read[read_i] == 'T' && tmp_ref[i] == 'C'))
				{
					error++;
					///break;

					if (error == 1)
					{
						*matched_length = read_i;
					}
					else
					{
						break;
					}
					
				}
			}

			read_i++;

		}

		///*matched_length = read_i;


	}


	*get_error = error;


	///if (*matched_length == read_length)
	if (error == 0)
	{
		sprintf(cigar, "%dM", read_length);


		output_sam_end_to_end(name, read, r_read, qulity,
			locates[0],
			read_length - 1,
			0,
			0,
			cigar, read_length);

		return 1;

	}

	return 0;



}



inline int try_process_unique_mismatch_end_to_end_pbat(
		bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
		bitmapper_bs_iter* locates,
		bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, int* get_error)
	{
		bitmapper_bs_iter tmp_SA_length = 0;


		locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, *matched_length, 0);



		int i = 0;
		int error = 0;


		if (*matched_length > first_C_site)
		{
			*matched_length = first_C_site;
		}


		if (*matched_length != read_length)
		{
			int need_read_length = read_length - *matched_length;

			if (locates[0] < _msf_refGenLength)
			{
				get_actuall_genome(tmp_ref, locates[0] + *matched_length, need_read_length);


			}
			else
			{

				get_actuall_rc_genome(tmp_ref, locates[0] + *matched_length - _msf_refGenLength,
					need_read_length);

			}

			int read_i = *matched_length;

			for (i = 0; i < need_read_length; i++)
			{
				if (read[read_i] != tmp_ref[i])
				{
					if (!(read[read_i] == 'T' && tmp_ref[i] == 'C'))
					{
						error++;
						///break;

						if (error == 1)
						{
							*matched_length = read_i;
						}
						else
						{
							break;
						}

					}
				}

				read_i++;

			}

			///*matched_length = read_i;


		}


		*get_error = error;


		///if (*matched_length == read_length)
		if (error == 0)
		{
			sprintf(cigar, "%dM", read_length);

			output_sam_end_to_end_pbat(name, read, r_read, qulity,
				locates[0],
				read_length - 1,
				0,
				0,
				cigar, read_length);


			return 1;

		}

		return 0;



}






inline int try_process_unique_mismatch_end_to_end_output_buffer(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, 
	Output_buffer_sub_block* sub_block, int* get_error)
{
	bitmapper_bs_iter tmp_SA_length = 0;

	locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, *matched_length, 0);

	int i = 0;
	int error = 0;

	if (*matched_length > first_C_site)
	{
		*matched_length = first_C_site;
	}


	if (*matched_length != read_length)
	{
		int need_read_length = read_length - *matched_length;

		if (locates[0] < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, locates[0] + *matched_length, need_read_length);


		}
		else
		{
			get_actuall_rc_genome(tmp_ref, locates[0] + *matched_length - _msf_refGenLength,
				need_read_length);

		}





		int read_i = *matched_length;

		for (i = 0; i < need_read_length; i++)
		{
			if (read[read_i] != tmp_ref[i])
			{
				if (!(read[read_i] == 'T' && tmp_ref[i] == 'C'))
				{
					error++;
					///break;

					if (error == 1)
					{
						*matched_length = read_i;
					}
					else
					{
						break;
					}

				}
			}

			read_i++;

		}





	}


	*get_error = error;


	if (error == 0)
	{
		sprintf(cigar, "%dM", read_length);



		output_sam_end_to_end_output_buffer(name, read, r_read, qulity,
			locates[0],
			read_length - 1,
			0,
			0,
			cigar, read_length, sub_block);


		return 1;

	}

	return 0;



}













inline int try_process_unique_mismatch_end_to_end_pbat_get_output_buffer(
		bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
		bitmapper_bs_iter* locates,
		bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, 
		Output_buffer_sub_block* sub_block, int* get_error)
	{
		bitmapper_bs_iter tmp_SA_length = 0;

		locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, *matched_length, 0);

		int i = 0;
		int error = 0;

		if (*matched_length > first_C_site)
		{
			*matched_length = first_C_site;
		}


		if (*matched_length != read_length)
		{
			int need_read_length = read_length - *matched_length;

			if (locates[0] < _msf_refGenLength)
			{
				get_actuall_genome(tmp_ref, locates[0] + *matched_length, need_read_length);


			}
			else
			{
				get_actuall_rc_genome(tmp_ref, locates[0] + *matched_length - _msf_refGenLength,
					need_read_length);

			}





			int read_i = *matched_length;

			for (i = 0; i < need_read_length; i++)
			{
				if (read[read_i] != tmp_ref[i])
				{
					if (!(read[read_i] == 'T' && tmp_ref[i] == 'C'))
					{
						error++;
						///break;

						if (error == 1)
						{
							*matched_length = read_i;
						}
						else
						{
							break;
						}

					}
				}

				read_i++;

			}





		}


		*get_error = error;


		if (error == 0)
		{
			sprintf(cigar, "%dM", read_length);



			output_sam_end_to_end_pbat_output_buffer(name, read, r_read, qulity,
				locates[0],
				read_length - 1,
				0,
				0,
				cigar, read_length, sub_block);


			return 1;

		}

		return 0;



	}










inline int try_process_unique_mismatch_end_to_end_return_site(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site,
	bitmapper_bs_iter* return_site, unsigned int* return_err, int* get_error)
{
	






	bitmapper_bs_iter tmp_SA_length = 0;


	locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, *matched_length, 0);



	int i = 0;
	int error = 0;


	if (*matched_length > first_C_site)
	{
		*matched_length = first_C_site;
	}


	if (*matched_length != read_length)
	{
		int need_read_length = read_length - *matched_length;

		if (locates[0] < _msf_refGenLength)
		{
			get_actuall_genome(tmp_ref, locates[0] + *matched_length, need_read_length);


		}
		else
		{

			get_actuall_rc_genome(tmp_ref, locates[0] + *matched_length - _msf_refGenLength,
				need_read_length);

		}

		int read_i = *matched_length;

		for (i = 0; i < need_read_length; i++)
		{
			if (read[read_i] != tmp_ref[i])
			{
				if (!(read[read_i] == 'T' && tmp_ref[i] == 'C'))
				{
					error++;
					///break;

					if (error == 1)
					{
						*matched_length = read_i;
					}
					else
					{
						break;
					}

				}
			}

			read_i++;

		}

		///*matched_length = read_i;


	}


	*get_error = error;


	///if (*matched_length == read_length)
	if (error == 0)
	{
		(*return_site) = locates[0];
		(*return_err) = 0;

		return 1;

	}

	return 0;

}















inline int new_faster_verify_pairs(int best_mapp_occ1, int best_mapp_occ2, int error_threshold,
	seed_votes* result1_array, seed_votes* result2_array,
	long long* best_pair_1_index_return, long long* best_pair_2_index_return,
	int inner_maxDistance_pair,
	int inner_minDistance_pair
	)
{
	int mapping_pair = 0;
	int best_sum_err = 2 * error_threshold + 1;
	int inner_i, inner_j;
	long long distance;
	long long current_sum_err;
	long long best_pair_1_index, best_pair_2_index;
	long long first_index_i;
	long long second_index_i;

	if (best_mapp_occ1 > 0 && best_mapp_occ2 > 0)
	{


			first_index_i = 0;

			///外层循环是result1_array
			for (inner_i = 0; inner_i < best_mapp_occ1; inner_i++)
			{
				///内层循环是result2_array
				for (inner_j = first_index_i; inner_j < best_mapp_occ2; inner_j++)
				{

					///在内层循环里，result1_array[inner_i].site是不变的
					///变的是result2_array[inner_j].site
					if (result1_array[inner_i].site > result2_array[inner_j].site)
					{
						distance = result1_array[inner_i].site - result2_array[inner_j].site;

						if (distance>inner_maxDistance_pair)
						{
							first_index_i = inner_j + 1;
						}
						else if (distance >= inner_minDistance_pair)
						{

							current_sum_err = result1_array[inner_i].err + result2_array[inner_j].err;

							/**
							if (current_sum_err == 1 && best_mapp_occ1 == 3 && best_mapp_occ2 == 2)
							{
								fprintf(stderr, "result1_array[%lld].site: %lld\n", 
									inner_i, result1_array[inner_i].site);
								fprintf(stderr, "result2_array[%lld].site: %lld\n",
									inner_j, result2_array[inner_j].site);

								fprintf(stderr, "distance: %lld\n", distance);
							}
							**/
						
						
							if (distance <= maxDistance_pair + current_sum_err
								&&
								distance >= minDistance_pair - current_sum_err)
							{


								if (current_sum_err < best_sum_err)
								{
									best_sum_err = current_sum_err;
									best_pair_1_index = inner_i;
									best_pair_2_index = inner_j;
									mapping_pair = 1;
								}
								else if (current_sum_err == best_sum_err)
								{
									mapping_pair++;

									if (best_sum_err == 0)
									{
										return mapping_pair;
									}
								}
							}





						}




					}
					else
					{
						distance = result2_array[inner_j].site - result1_array[inner_i].site;

						if (distance > inner_maxDistance_pair)
						{
							break;
						}




						if (distance <= inner_maxDistance_pair
							&&
							distance >= inner_minDistance_pair)
						{

							current_sum_err = result1_array[inner_i].err + result2_array[inner_j].err;

							/**
							if (current_sum_err == 1 && best_mapp_occ1 == 3 && best_mapp_occ2 == 2)
							{
								fprintf(stderr, "result1_array[%lld].site: %lld\n",
									inner_i, result1_array[inner_i].site);
								fprintf(stderr, "result2_array[%lld].site: %lld\n",
									inner_j, result2_array[inner_j].site);

								fprintf(stderr, "distance: %lld\n", distance);
							}
							**/


							if (distance <= maxDistance_pair + current_sum_err
								&&
								distance >= minDistance_pair - current_sum_err)
							{

								if (current_sum_err < best_sum_err)
								{
									best_sum_err = current_sum_err;
									best_pair_1_index = inner_i;
									best_pair_2_index = inner_j;
									mapping_pair = 1;
								}
								else if (current_sum_err == best_sum_err)
								{
									mapping_pair++;

									if (best_sum_err == 0)
									{
										return mapping_pair;
									}

								}
							}
						}



					}








				}
			}

		




	}


	(*best_pair_1_index_return) = best_pair_1_index;
	(*best_pair_2_index_return) = best_pair_2_index;

	return mapping_pair;
}














inline void filter_pairs(
	int best_mapp_occ1, 
	int best_mapp_occ2, 
	seed_votes** result1_array, 
	seed_votes** result2_array,
	seed_votes** return1_array,
	seed_votes** return2_array,
	bitmapper_bs_iter* return_length1,
	bitmapper_bs_iter* return_length2,
	int inner_maxDistance_pair,
	int inner_minDistance_pair
	)
{
	int inner_i, inner_j;
	long long distance;
	long long first_index_i;
	long long second_index_i;
	long long length1, length2;
	length1 = 0;
	length2 = 0;

	
	///注意在这种情况下, result1_array和result2_array都是有序的
	if (best_mapp_occ1 > 0 && best_mapp_occ2 > 0)
	{


		first_index_i = 0;

		///外层循环是result1_array
		for (inner_i = 0; inner_i < best_mapp_occ1; inner_i++)
		{
			///内层循环是result2_array
			for (inner_j = first_index_i; inner_j < best_mapp_occ2; inner_j++)
			{

				///在内层循环里，result1_array[inner_i].site是不变的
				///变的是result2_array[inner_j].site
				///逻辑实际上是把result2_array中所有元素依次去和result1_array[inner_i].site去比
				if ((*result1_array)[inner_i].site >(*result2_array)[inner_j].site)
				{
					distance = (*result1_array)[inner_i].site - (*result2_array)[inner_j].site;

					if (distance>inner_maxDistance_pair)
					{
						first_index_i = inner_j + 1;
					}
					else if (distance >= inner_minDistance_pair)
					{
						///这样写不会溢出
						//因为如果length1 == 0, 就不会去判断第二个条件 
						if (length1 == 0 || (*result1_array)[inner_i].site > (*return1_array)[length1 - 1].site)
						{
							(*return1_array)[length1].site = (*result1_array)[inner_i].site;
							(*return1_array)[length1].err = (*result1_array)[inner_i].err;
							(*return1_array)[length1].end_site = (*result1_array)[inner_i].end_site;
							length1++;
						}

						if (length2 == 0 || (*result2_array)[inner_j].site > (*return2_array)[length2 - 1].site)
						{
							(*return2_array)[length2].site = (*result2_array)[inner_j].site;
							(*return2_array)[length2].err = (*result2_array)[inner_j].err;
							(*return2_array)[length2].end_site = (*result2_array)[inner_j].end_site;
							length2++;
						}
						
					}

				}
				else
				{
					distance = (*result2_array)[inner_j].site - (*result1_array)[inner_i].site;

					if (distance > inner_maxDistance_pair)
					{
						break;
					}




					if (distance <= inner_maxDistance_pair
						&&
						distance >= inner_minDistance_pair)
					{

						///这样写不会溢出
						//因为如果length1 == 0, 就不会去判断第二个条件 
						if (length1 == 0 || (*result1_array)[inner_i].site > (*return1_array)[length1 - 1].site)
						{
							(*return1_array)[length1].site = (*result1_array)[inner_i].site;
							(*return1_array)[length1].err = (*result1_array)[inner_i].err;
							(*return1_array)[length1].end_site = (*result1_array)[inner_i].end_site;
							length1++;
						}

						if (length2 == 0 || (*result2_array)[inner_j].site > (*return2_array)[length2 - 1].site)
						{
							(*return2_array)[length2].site = (*result2_array)[inner_j].site;
							(*return2_array)[length2].err = (*result2_array)[inner_j].err;
							(*return2_array)[length2].end_site = (*result2_array)[inner_j].end_site;
							length2++;
						}
					}

				}
			}
		}
	}
	


	(*return_length1) = length1;
	(*return_length2) = length2;

	seed_votes* k;

	k = *result1_array;
	*result1_array = *return1_array;
	*return1_array = k;


	k = *result2_array;
	*result2_array = *return2_array;
	*return2_array = k;

	
}




///用1的结果过滤2
inline void filter_pairs_single_side(
	int best_mapp_occ1,
	int best_mapp_occ2,
	seed_votes** result1_array,
	seed_votes** result2_array,
	bitmapper_bs_iter* return_length2,
	int inner_maxDistance_pair,
	int inner_minDistance_pair
	)
{
	int inner_i, inner_j;
	long long distance;
	long long first_index_i;
	long long second_index_i;
	long long length1, length2;
	length2 = 0;


	///注意在这种情况下, result1_array和result2_array都是有序的
	if (best_mapp_occ1 > 0 && best_mapp_occ2 > 0)
	{


		first_index_i = 0;

		///外层循环是result1_array
		for (inner_i = 0; inner_i < best_mapp_occ1; inner_i++)
		{
			///内层循环是result2_array
			for (inner_j = first_index_i; inner_j < best_mapp_occ2; inner_j++)
			{

				///在内层循环里，result1_array[inner_i].site是不变的
				///变的是result2_array[inner_j].site
				///逻辑实际上是把result2_array中所有元素依次去和result1_array[inner_i].site去比
				if ((*result1_array)[inner_i].site >(*result2_array)[inner_j].site)
				{
					distance = (*result1_array)[inner_i].site - (*result2_array)[inner_j].site;

					if (distance>inner_maxDistance_pair)
					{
						first_index_i = inner_j + 1;
					}
					else if (distance >= inner_minDistance_pair)
					{
						///这样写不会溢出
						//因为如果length2 == 0, 就不会去判断第二个条件 
						////if (length2 == 0 || (*result2_array)[inner_j].site > (*return2_array)[length2 - 1].site)
						{
							(*result2_array)[length2].site = (*result2_array)[inner_j].site;
							(*result2_array)[length2].err = (*result2_array)[inner_j].err;
							(*result2_array)[length2].end_site = (*result2_array)[inner_j].end_site;


							length2++;
						}

						first_index_i = inner_j + 1;

					}

				}
				else
				{
					distance = (*result2_array)[inner_j].site - (*result1_array)[inner_i].site;

					if (distance > inner_maxDistance_pair)
					{
						break;
					}




					if (distance <= inner_maxDistance_pair
						&&
						distance >= inner_minDistance_pair)
					{

						///这样写不会溢出
						//因为如果length2 == 0, 就不会去判断第二个条件 
						///if (length2 == 0 || (*result2_array)[inner_j].site > (*return2_array)[length2 - 1].site)
						{
							(*result2_array)[length2].site = (*result2_array)[inner_j].site;
							(*result2_array)[length2].err = (*result2_array)[inner_j].err;
							(*result2_array)[length2].end_site = (*result2_array)[inner_j].end_site;

							length2++;
						}

						first_index_i = inner_j + 1;
					}

				}
			}
		}
	}



	(*return_length2) = length2;

}









int inline process_rest_seed_filter_debug(
	bitmapper_bs_iter seed_id1,
	bitmapper_bs_iter max_seed_number1,
	bitmapper_bs_iter total_match_length1,
	Read* current_read,
	bitmapper_bs_iter current_seed_length1,
	bitmapper_bs_iter number_of_hits1,
	char* bsSeq1,
	bitmapper_bs_iter match_length1,
	bitmapper_bs_iter* locates1,
	bitmapper_bs_iter* candidates1,
	bitmapper_bs_iter candidate_length1,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	seed_votes* candidates_votes1,
	bitmapper_bs_iter candidates_votes_length1,
	bitmapper_bs_iter error_threshold,
	char* tmp_ref,
	__m256i* Peq_SSE,
	unsigned int min_err,
	seed_votes* candidates_votes2,
	int best_mapp_occ2,
	int inner_maxDistance_pair,
	int inner_minDistance_pair,
	int get_error1,
	long long first_seed_match_length1,
	int* read1_seed_start,
	int* read1_seed_length,
	bitmapper_bs_iter* full_seed_id1)
{
	int best_mapp_occ1;
	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter tmp_SA_length;

	long long second_seed_length1;
	int extra_seed_flag1;


	///这里要改
	extra_seed_flag1 = 1;

	/**
	if (strcmp("@SRR1532535.50167", (*current_read).name) == 0)
	{
		fprintf(stderr, "\n\n%s(2)\n", (*current_read).name);
		fprintf(stderr, "get_error1: %d\n", get_error1);
	}
	**/

	///这代表已经找到一个匹配位置了
	///所以不需要记录seed的分布
	if (get_error1 == 1)
	{

		second_seed_length1 = (*current_read).length - first_seed_match_length1;


		if (second_seed_length1 >= 17)
		{

			number_of_hits1 = count_hash_table(bsSeq1, second_seed_length1, &top, &bot, &pre_top, &pre_bot);

			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length,
					second_seed_length1, first_seed_match_length1);


				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;

				extra_seed_flag1 = 0;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq1, top, bot, pre_top, pre_bot, locates1, second_seed_length1, &tmp_SA_length);

					reverse_and_adjust_site(locates1, &tmp_SA_length, second_seed_length1, first_seed_match_length1);


					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

				extra_seed_flag1 = 0;

			}
			else if (number_of_hits1 > max_seed_matches)
			{
				extra_seed_flag1 = 1;
			}


		}
		else
		{
			extra_seed_flag1 = 1;
		}


	}




	///这里要改
	if (extra_seed_flag1 == 1)
	{


		///这是剩下的seed
		while (seed_id1 < max_seed_number1 && total_match_length1 < (*current_read).length)
		{



			current_seed_length1 = (*current_read).length - total_match_length1;


			number_of_hits1 =
				count_backward_as_much_1_terminate(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);


			read1_seed_start[*full_seed_id1] = total_match_length1;
			read1_seed_length[*full_seed_id1] = match_length1;
			(*full_seed_id1)++;


			/**
			if (strcmp("@SRR1532535.50167", (*current_read).name) == 0)
			{

				fprintf(stderr, "seed_id2: %lld\n", seed_id1);
				fprintf(stderr, "match_length2: %lld\n", match_length1);
				fprintf(stderr, "number_of_hits2: %lld\n", number_of_hits1);
				fprintf(stderr, "candidate_length2: %lld\n", candidate_length1);
				fprintf(stderr, "total_match_length2: %lld\n", total_match_length1);


				fprintf(stderr, "?????\n");
			}
			**/




			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
					reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

			}
			else
			{

				(*full_seed_id1)--;

				if (current_seed_length1 == match_length1)
				{
					break;
				}
			}









			if (match_length1 == 0)
			{
				total_match_length1 = total_match_length1 + match_step;
			}
			else
			{
				total_match_length1 = total_match_length1 + match_length1 / 2;
			}




			seed_id1++;





		}
	}


	/**
	if (strcmp("@SRR1532535.4275", (*current_read).name) == 0)
	{
	for (size_t haha_i = 0; haha_i < candidate_length1; haha_i++)
	{
	fprintf(stderr, "candidates1[%lld]: %lld\n", haha_i, candidates1[haha_i]);

	}

	fprintf(stderr, "####################\n");
	}
	**/

	///这里要改
	///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
	if (extra_seed_flag1 == 0
		&&
		(candidate_length1 == 1 ||
		(candidate_length1 == 2 && candidates1[0] == candidates1[1])
		))
	{

		candidates_votes1[0].site = candidates1[0];
		candidates_votes1[0].err = 1;
		candidates_votes1[0].end_site = (*current_read).length - 1;
		best_mapp_occ1 = 1;
	}
	else if (candidate_length1 != 0)
	{

		std::sort(candidates1, candidates1 + candidate_length1);

		generate_candidate_votes_shift_filter(candidates1, candidate_length1, candidates_votes1,
			&candidates_votes_length1, error_threshold, candidates_votes2, best_mapp_occ2, inner_maxDistance_pair, inner_minDistance_pair);


		if (error_threshold <= 15)
		{
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}
		else
		{
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}



	}
	else
	{
		///没有种子
		best_mapp_occ1 = 0;
		///return_flag1 = 0;
	}


	/**
	if (strcmp("@SRR1532535.50167", (*current_read).name) == 0)
	{
		///fprintf(stderr, "\n\n%s(2)\n", (*current_read).name);
		fprintf(stderr, "best_mapp_occ2: %d\n", best_mapp_occ1);
	}
	**/


	return best_mapp_occ1;







}







int inline select_best_seeds
(int* read1_seed_start, int* read1_seed_length,bitmapper_bs_iter full_seed_id1,
int* result_seed_start, int* result_seed_length, int max_seed_number1, int read_length)
{
	int result_seed_number = 0;


	if (full_seed_id1 >= 2)
	{
		result_seed_number = 2;

		result_seed_start[0] = read1_seed_start[0];
		result_seed_length[0] = read1_seed_start[1] - read1_seed_start[0];

		result_seed_start[1]
			= read1_seed_start[full_seed_id1 - 2] + read1_seed_length[full_seed_id1 - 2];
		result_seed_length[1] = read_length - result_seed_start[1];
	}
	else if (full_seed_id1 == 1)///if (full_seed_id1 == 1 && read1_seed_length[0] == read_length)
	{
		result_seed_number = 2;

		result_seed_start[0] = read1_seed_start[0];
		result_seed_length[0] = read_length / 2;

		result_seed_start[1] = result_seed_start[0] + result_seed_length[0];
		result_seed_length[1] = read_length - result_seed_start[1];
	}
	
	if (read1_seed_start[full_seed_id1 - 1] + read1_seed_length[full_seed_id1 - 1] < read_length)
	{
		result_seed_start[result_seed_number]
			= read1_seed_start[full_seed_id1 - 1] + read1_seed_length[full_seed_id1 - 1];
		result_seed_length[result_seed_number] = read_length - result_seed_start[result_seed_number];
		result_seed_number++;
	}
	



	return result_seed_number;
}






int inline reseed_filter(
	bitmapper_bs_iter seed_id1,
	bitmapper_bs_iter max_seed_number1,
	bitmapper_bs_iter total_match_length1,
	Read* current_read,
	bitmapper_bs_iter current_seed_length1,
	bitmapper_bs_iter number_of_hits1,
	char* bsSeq1,
	bitmapper_bs_iter match_length1,
	bitmapper_bs_iter* locates1,
	bitmapper_bs_iter* candidates1,
	bitmapper_bs_iter candidate_length1,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	seed_votes* candidates_votes1,
	bitmapper_bs_iter candidates_votes_length1,
	bitmapper_bs_iter error_threshold,
	char* tmp_ref,
	__m256i* Peq_SSE,
	unsigned int min_err,
	seed_votes* candidates_votes2,
	int best_mapp_occ2,
	int inner_maxDistance_pair,
	int inner_minDistance_pair,
	int get_error1,
	long long first_seed_match_length1,
	int* read1_seed_start,
	int* read1_seed_length,
	bitmapper_bs_iter full_seed_id1,
	int* result_seed_start,
	int* result_seed_length)
{
	int best_mapp_occ1;
	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter tmp_SA_length;

	char* tmp_bsSeq1;

	long long second_seed_length1;




	available_seed_length = 20;

	seed_id1 = 0;
	candidate_length1 = 0;
	locates1 = candidates1;
	
	int result_seed_number;
	
	result_seed_number = select_best_seeds(read1_seed_start, read1_seed_length, full_seed_id1,
		result_seed_start, result_seed_length, max_seed_number1, (*current_read).length);



	

	while (seed_id1 < result_seed_number)
	{
		total_match_length1 = result_seed_start[seed_id1];
		current_seed_length1 = (*current_read).length - total_match_length1;
		match_length1 = result_seed_length[seed_id1];

		//tmp_bsSeq1 = bsSeq1 + (*current_read).length - read1_seed_start[seed_id1] - read1_seed_length[seed_id1];

		tmp_bsSeq1 = bsSeq1 + current_seed_length1 - match_length1;

		number_of_hits1 = count_hash_table(tmp_bsSeq1, match_length1, &top, &bot, &pre_top, &pre_bot);

		///如果当前seed的只有一个匹配位置
		///则不管多长都是有效的
		if (number_of_hits1 == 1)
		{
			tmp_SA_length = 0;

			locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

			candidate_length1 = candidate_length1 + tmp_SA_length;

			locates1 = candidates1 + candidate_length1;


		}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
		else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
		{

			if (number_of_hits1 != 0)
			{

				tmp_SA_length = 0;

				///locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
				locate(tmp_bsSeq1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
				reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;
				locates1 = candidates1 + candidate_length1;

			}

		}
		else if (current_seed_length1 == match_length1)
		{
			break;
		}


		seed_id1++;

	}
	
	
	
	///available_seed_length = 20;

	if (full_seed_id1 > 1)
	{
		total_match_length1 = (read1_seed_start[0] + read1_seed_start[1])/2;
	}
	else
	{
		total_match_length1 = match_step / 2;
	}

	/**

	while (seed_id1 < max_seed_number1 && total_match_length1 < (*current_read).length)
	{



		current_seed_length1 = (*current_read).length - total_match_length1;
		match_length1 = available_seed_length;

		//tmp_bsSeq1 = bsSeq1 + (*current_read).length - read1_seed_start[seed_id1] - read1_seed_length[seed_id1];

		tmp_bsSeq1 = bsSeq1 + current_seed_length1 - match_length1;

		number_of_hits1 = count_hash_table(tmp_bsSeq1, match_length1, &top, &bot, &pre_top, &pre_bot);



		///如果当前seed的只有一个匹配位置
		///则不管多长都是有效的
		if (number_of_hits1 == 1)
		{
			tmp_SA_length = 0;

			locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

			candidate_length1 = candidate_length1 + tmp_SA_length;

			locates1 = candidates1 + candidate_length1;


		}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
		else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
		{

			if (number_of_hits1 != 0)
			{

				tmp_SA_length = 0;

				///locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
				locate(tmp_bsSeq1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
				reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;
				locates1 = candidates1 + candidate_length1;

			}

		}
		else if (current_seed_length1 == match_length1)
		{
			break;
		}





		total_match_length1 = total_match_length1 + match_step;

		seed_id1++;

	}

	**/




	
	while (seed_id1 < max_seed_number1 && total_match_length1 < (*current_read).length)
	{



		current_seed_length1 = (*current_read).length - total_match_length1;


		number_of_hits1 =
			count_backward_as_much_1_terminate(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);


	
		


		///如果当前seed的只有一个匹配位置
		///则不管多长都是有效的
		if (number_of_hits1 == 1)
		{
			tmp_SA_length = 0;

			locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

			candidate_length1 = candidate_length1 + tmp_SA_length;

			locates1 = candidates1 + candidate_length1;


		}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
		else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
		{

			if (number_of_hits1 != 0)
			{

				tmp_SA_length = 0;

				locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
				reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;
				locates1 = candidates1 + candidate_length1;

			}

		}
		else if (current_seed_length1 == match_length1)
		{
			break;
		}

		total_match_length1 = total_match_length1 + match_step;

		seed_id1++;

	}
	
	
	
	
	


	if (candidate_length1 != 0)
	{

		std::sort(candidates1, candidates1 + candidate_length1);

		generate_candidate_votes_shift_filter(candidates1, candidate_length1, candidates_votes1,
			&candidates_votes_length1, error_threshold, candidates_votes2, best_mapp_occ2, inner_maxDistance_pair, inner_minDistance_pair);


		if (error_threshold <= 15)
		{
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}
		else
		{
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}



	}
	else
	{
		///没有种子
		best_mapp_occ1 = 0;
		///return_flag1 = 0;
	}

	///fprintf(stderr, "best_mapp_occ1: %d\n", best_mapp_occ1);

	return best_mapp_occ1;
}








int inline reseed_filter_muti_thread(
	bitmapper_bs_iter seed_id1,
	bitmapper_bs_iter max_seed_number1,
	bitmapper_bs_iter total_match_length1,
	Read* current_read,
	bitmapper_bs_iter current_seed_length1,
	bitmapper_bs_iter number_of_hits1,
	char* bsSeq1,
	bitmapper_bs_iter match_length1,
	bitmapper_bs_iter* locates1,
	bitmapper_bs_iter* candidates1,
	bitmapper_bs_iter candidate_length1,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	seed_votes* candidates_votes1,
	bitmapper_bs_iter candidates_votes_length1,
	bitmapper_bs_iter error_threshold,
	char* tmp_ref,
	__m256i* Peq_SSE,
	unsigned int min_err,
	seed_votes* candidates_votes2,
	int best_mapp_occ2,
	int inner_maxDistance_pair,
	int inner_minDistance_pair,
	bwt_locate_queue* get_queue,
	int get_error1,
	long long first_seed_match_length1,
	int* read1_seed_start,
	int* read1_seed_length,
	bitmapper_bs_iter full_seed_id1,
	int* result_seed_start,
	int* result_seed_length)
{
	int best_mapp_occ1;
	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter tmp_SA_length;

	char* tmp_bsSeq1;

	long long second_seed_length1;




	available_seed_length = 20;

	seed_id1 = 0;
	candidate_length1 = 0;
	locates1 = candidates1;

	int result_seed_number;

	result_seed_number = select_best_seeds(read1_seed_start, read1_seed_length, full_seed_id1,
		result_seed_start, result_seed_length, max_seed_number1, (*current_read).length);





	while (seed_id1 < result_seed_number)
	{
		total_match_length1 = result_seed_start[seed_id1];
		current_seed_length1 = (*current_read).length - total_match_length1;
		match_length1 = result_seed_length[seed_id1];

		//tmp_bsSeq1 = bsSeq1 + (*current_read).length - read1_seed_start[seed_id1] - read1_seed_length[seed_id1];

		tmp_bsSeq1 = bsSeq1 + current_seed_length1 - match_length1;

		number_of_hits1 = count_hash_table(tmp_bsSeq1, match_length1, &top, &bot, &pre_top, &pre_bot);

		///如果当前seed的只有一个匹配位置
		///则不管多长都是有效的
		if (number_of_hits1 == 1)
		{
			tmp_SA_length = 0;

			locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

			candidate_length1 = candidate_length1 + tmp_SA_length;

			locates1 = candidates1 + candidate_length1;


		}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
		else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
		{

			if (number_of_hits1 != 0)
			{

				tmp_SA_length = 0;

				///locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
				locate_muti_thread(tmp_bsSeq1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length, get_queue);
				reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;
				locates1 = candidates1 + candidate_length1;

			}

		}
		else if (current_seed_length1 == match_length1)
		{
			break;
		}


		seed_id1++;

	}



	///available_seed_length = 20;

	if (full_seed_id1 > 1)
	{
		total_match_length1 = (read1_seed_start[0] + read1_seed_start[1]) / 2;
	}
	else
	{
		total_match_length1 = match_step / 2;
	}


	while (seed_id1 < max_seed_number1 && total_match_length1 < (*current_read).length)
	{



		current_seed_length1 = (*current_read).length - total_match_length1;


		number_of_hits1 =
			count_backward_as_much_1_terminate(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);






		///如果当前seed的只有一个匹配位置
		///则不管多长都是有效的
		if (number_of_hits1 == 1)
		{
			tmp_SA_length = 0;

			locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

			candidate_length1 = candidate_length1 + tmp_SA_length;

			locates1 = candidates1 + candidate_length1;


		}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
		else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
		{

			if (number_of_hits1 != 0)
			{

				tmp_SA_length = 0;

				locate_muti_thread(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length, get_queue);
				reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;
				locates1 = candidates1 + candidate_length1;

			}

		}
		else if (current_seed_length1 == match_length1)
		{
			break;
		}

		total_match_length1 = total_match_length1 + match_step;

		seed_id1++;

	}






	if (candidate_length1 != 0)
	{

		std::sort(candidates1, candidates1 + candidate_length1);

		generate_candidate_votes_shift_filter(candidates1, candidate_length1, candidates_votes1,
			&candidates_votes_length1, error_threshold, candidates_votes2, best_mapp_occ2, inner_maxDistance_pair, inner_minDistance_pair);


		if (error_threshold <= 15)
		{
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}
		else
		{
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}



	}
	else
	{
		///没有种子
		best_mapp_occ1 = 0;
		///return_flag1 = 0;
	}

	///fprintf(stderr, "best_mapp_occ1: %d\n", best_mapp_occ1);

	return best_mapp_occ1;
}




int inline process_rest_seed_filter_muti_thread(
	bitmapper_bs_iter seed_id1,
	bitmapper_bs_iter max_seed_number1,
	bitmapper_bs_iter total_match_length1,
	Read* current_read,
	bitmapper_bs_iter current_seed_length1,
	bitmapper_bs_iter number_of_hits1,
	char* bsSeq1,
	bitmapper_bs_iter match_length1,
	bitmapper_bs_iter* locates1,
	bitmapper_bs_iter* candidates1,
	bitmapper_bs_iter candidate_length1,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	seed_votes* candidates_votes1,
	bitmapper_bs_iter candidates_votes_length1,
	bitmapper_bs_iter error_threshold,
	char* tmp_ref,
	__m256i* Peq_SSE,
	unsigned int min_err,
	seed_votes* candidates_votes2,
	int best_mapp_occ2,
	int inner_maxDistance_pair,
	int inner_minDistance_pair,
	bwt_locate_queue* get_queue,
	int get_error1,
	long long first_seed_match_length1,
	int* read1_seed_start,
	int* read1_seed_length,
	bitmapper_bs_iter* full_seed_id1)
{
	int best_mapp_occ1;
	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter tmp_SA_length;

	long long second_seed_length1;
	int extra_seed_flag1;


	///这里要改
	extra_seed_flag1 = 1;

	/**
	if (strcmp("@SRR1532535.50167", (*current_read).name) == 0)
	{
	fprintf(stderr, "\n\n%s(2)\n", (*current_read).name);
	fprintf(stderr, "get_error1: %d\n", get_error1);
	}
	**/

	///这代表已经找到一个匹配位置了
	///所以不需要记录seed的分布
	if (get_error1 == 1)
	{

		second_seed_length1 = (*current_read).length - first_seed_match_length1;


		if (second_seed_length1 >= 17)
		{

			number_of_hits1 = count_hash_table(bsSeq1, second_seed_length1, &top, &bot, &pre_top, &pre_bot);

			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length,
					second_seed_length1, first_seed_match_length1);


				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;

				extra_seed_flag1 = 0;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate_muti_thread(bsSeq1, top, bot, pre_top, pre_bot, locates1, second_seed_length1, &tmp_SA_length, get_queue);

					reverse_and_adjust_site(locates1, &tmp_SA_length, second_seed_length1, first_seed_match_length1);


					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

				extra_seed_flag1 = 0;

			}
			else if (number_of_hits1 > max_seed_matches)
			{
				extra_seed_flag1 = 1;
			}


		}
		else
		{
			extra_seed_flag1 = 1;
		}


	}




	///这里要改
	if (extra_seed_flag1 == 1)
	{


		///这是剩下的seed
		while (seed_id1 < max_seed_number1 && total_match_length1 < (*current_read).length)
		{



			current_seed_length1 = (*current_read).length - total_match_length1;


			number_of_hits1 =
				count_backward_as_much_1_terminate(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);


			read1_seed_start[*full_seed_id1] = total_match_length1;
			read1_seed_length[*full_seed_id1] = match_length1;
			(*full_seed_id1)++;


			/**
			if (strcmp("@SRR1532535.50167", (*current_read).name) == 0)
			{

			fprintf(stderr, "seed_id2: %lld\n", seed_id1);
			fprintf(stderr, "match_length2: %lld\n", match_length1);
			fprintf(stderr, "number_of_hits2: %lld\n", number_of_hits1);
			fprintf(stderr, "candidate_length2: %lld\n", candidate_length1);
			fprintf(stderr, "total_match_length2: %lld\n", total_match_length1);


			fprintf(stderr, "?????\n");
			}
			**/




			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate_muti_thread(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length, get_queue);
					reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

			}
			else
			{

				(*full_seed_id1)--;

				if (current_seed_length1 == match_length1)
				{
					break;
				}
			}









			if (match_length1 == 0)
			{
				total_match_length1 = total_match_length1 + match_step;
			}
			else
			{
				total_match_length1 = total_match_length1 + match_length1 / 2;
			}




			seed_id1++;





		}
	}


	/**
	if (strcmp("@SRR1532535.4275", (*current_read).name) == 0)
	{
	for (size_t haha_i = 0; haha_i < candidate_length1; haha_i++)
	{
	fprintf(stderr, "candidates1[%lld]: %lld\n", haha_i, candidates1[haha_i]);

	}

	fprintf(stderr, "####################\n");
	}
	**/

	///这里要改
	///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
	if (extra_seed_flag1 == 0
		&&
		(candidate_length1 == 1 ||
		(candidate_length1 == 2 && candidates1[0] == candidates1[1])
		))
	{

		candidates_votes1[0].site = candidates1[0];
		candidates_votes1[0].err = 1;
		candidates_votes1[0].end_site = (*current_read).length - 1;
		best_mapp_occ1 = 1;
	}
	else if (candidate_length1 != 0)
	{

		std::sort(candidates1, candidates1 + candidate_length1);

		generate_candidate_votes_shift_filter(candidates1, candidate_length1, candidates_votes1,
			&candidates_votes_length1, error_threshold, candidates_votes2, best_mapp_occ2, inner_maxDistance_pair, inner_minDistance_pair);


		if (error_threshold <= 15)
		{
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}
		else
		{
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}



	}
	else
	{
		///没有种子
		best_mapp_occ1 = 0;
		///return_flag1 = 0;
	}


	/**
	if (strcmp("@SRR1532535.50167", (*current_read).name) == 0)
	{
	///fprintf(stderr, "\n\n%s(2)\n", (*current_read).name);
	fprintf(stderr, "best_mapp_occ2: %d\n", best_mapp_occ1);
	}
	**/


	return best_mapp_occ1;







}




int inline process_rest_seed_debug(
	bitmapper_bs_iter seed_id1,
	bitmapper_bs_iter max_seed_number1,
	bitmapper_bs_iter total_match_length1,
	Read* current_read,
	bitmapper_bs_iter current_seed_length1,
	bitmapper_bs_iter number_of_hits1,
	char* bsSeq1,
	bitmapper_bs_iter match_length1,
	bitmapper_bs_iter* locates1,
	bitmapper_bs_iter* candidates1,
	bitmapper_bs_iter candidate_length1,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	seed_votes* candidates_votes1,
	bitmapper_bs_iter candidates_votes_length1,
	bitmapper_bs_iter error_threshold,
	char* tmp_ref,
	__m256i* Peq_SSE,
	unsigned int min_err,
	int get_error1,
	long long first_seed_match_length1,
	int* read1_seed_start,
	int* read1_seed_length,
	bitmapper_bs_iter* full_seed_id1
	)
{
	int best_mapp_occ1;
	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter tmp_SA_length;

	long long second_seed_length1;
	int extra_seed_flag1;



	///这里要改
	extra_seed_flag1 = 1;
	///这代表已经找到一个匹配位置了
	///所以不需要记录seed的分布
	if (get_error1 == 1)
	{

		second_seed_length1 = (*current_read).length - first_seed_match_length1;


		if (second_seed_length1 >= 17)
		{

			number_of_hits1 = count_hash_table(bsSeq1, second_seed_length1, &top, &bot, &pre_top, &pre_bot);

			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length,
					second_seed_length1, first_seed_match_length1);


				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;

				extra_seed_flag1 = 0;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq1, top, bot, pre_top, pre_bot, locates1, second_seed_length1, &tmp_SA_length);

					reverse_and_adjust_site(locates1, &tmp_SA_length, second_seed_length1, first_seed_match_length1);


					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

				extra_seed_flag1 = 0;

			}
			else if (number_of_hits1 > max_seed_matches)
			{
				extra_seed_flag1 = 1;
			}


		}
		else
		{
			extra_seed_flag1 = 1;
		}


	}




	///这里要改
	if (extra_seed_flag1 == 1)
	{


		///这是剩下的seed
		while (seed_id1 < max_seed_number1 && total_match_length1 < (*current_read).length)
		{



			current_seed_length1 = (*current_read).length - total_match_length1;


			number_of_hits1 =
				count_backward_as_much_1_terminate(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);




			read1_seed_start[*full_seed_id1] = total_match_length1;
			read1_seed_length[*full_seed_id1] = match_length1;
			(*full_seed_id1)++;




			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
					reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

			}
			else
			{

				(*full_seed_id1)--;

				if (current_seed_length1 == match_length1)
				{
					break;
				}
			}
			






			if (match_length1 == 0)
			{
				total_match_length1 = total_match_length1 + match_step;
			}
			else
			{
				total_match_length1 = total_match_length1 + match_length1 / 2;
			}



			seed_id1++;


		}
	}


	///这里要改
	///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
	if (extra_seed_flag1 == 0
		&&
		(candidate_length1 == 1 ||
		(candidate_length1 == 2 && candidates1[0] == candidates1[1])
		))
	{

		candidates_votes1[0].site = candidates1[0];
		candidates_votes1[0].err = 1;
		candidates_votes1[0].end_site = (*current_read).length - 1;
		best_mapp_occ1 = 1;
	}
	else if (candidate_length1 != 0)
	{

		std::sort(candidates1, candidates1 + candidate_length1);

		generate_candidate_votes_shift(candidates1, candidate_length1, candidates_votes1, &candidates_votes_length1, error_threshold);


		if (error_threshold <= 15)
		{
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}
		else
		{
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}



	}
	else
	{
		///没有种子
		best_mapp_occ1 = 0;
		///return_flag1 = 0;
	}



	return best_mapp_occ1;







}






int inline process_rest_seed_muti_thread(
	bitmapper_bs_iter seed_id1,
	bitmapper_bs_iter max_seed_number1,
	bitmapper_bs_iter total_match_length1,
	Read* current_read,
	bitmapper_bs_iter current_seed_length1,
	bitmapper_bs_iter number_of_hits1,
	char* bsSeq1,
	bitmapper_bs_iter match_length1,
	bitmapper_bs_iter* locates1,
	bitmapper_bs_iter* candidates1,
	bitmapper_bs_iter candidate_length1,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	seed_votes* candidates_votes1,
	bitmapper_bs_iter candidates_votes_length1,
	bitmapper_bs_iter error_threshold,
	char* tmp_ref,
	__m256i* Peq_SSE,
	unsigned int min_err,
	bwt_locate_queue* get_queue,
	int get_error1,
	long long first_seed_match_length1,
	int* read1_seed_start,
	int* read1_seed_length,
	bitmapper_bs_iter* full_seed_id1
	)
{
	int best_mapp_occ1;
	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter tmp_SA_length;

	long long second_seed_length1;
	int extra_seed_flag1;



	///这里要改
	extra_seed_flag1 = 1;
	///这代表已经找到一个匹配位置了
	///所以不需要记录seed的分布
	if (get_error1 == 1)
	{

		second_seed_length1 = (*current_read).length - first_seed_match_length1;


		if (second_seed_length1 >= 17)
		{

			number_of_hits1 = count_hash_table(bsSeq1, second_seed_length1, &top, &bot, &pre_top, &pre_bot);

			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length,
					second_seed_length1, first_seed_match_length1);


				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;

				extra_seed_flag1 = 0;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate_muti_thread(bsSeq1, top, bot, pre_top, pre_bot, locates1, second_seed_length1, &tmp_SA_length, get_queue);

					reverse_and_adjust_site(locates1, &tmp_SA_length, second_seed_length1, first_seed_match_length1);


					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

				extra_seed_flag1 = 0;

			}
			else if (number_of_hits1 > max_seed_matches)
			{
				extra_seed_flag1 = 1;
			}


		}
		else
		{
			extra_seed_flag1 = 1;
		}


	}




	///这里要改
	if (extra_seed_flag1 == 1)
	{


		///这是剩下的seed
		while (seed_id1 < max_seed_number1 && total_match_length1 < (*current_read).length)
		{



			current_seed_length1 = (*current_read).length - total_match_length1;


			number_of_hits1 =
				count_backward_as_much_1_terminate(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);




			read1_seed_start[*full_seed_id1] = total_match_length1;
			read1_seed_length[*full_seed_id1] = match_length1;
			(*full_seed_id1)++;




			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			if (number_of_hits1 == 1)
			{
				tmp_SA_length = 0;

				locate_one_position_direct(locates1, top, &tmp_SA_length, total_SA_length, match_length1, total_match_length1);

				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate_muti_thread(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length, get_queue);
					reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;

				}

			}
			else
			{

				(*full_seed_id1)--;

				if (current_seed_length1 == match_length1)
				{
					break;
				}
			}







			if (match_length1 == 0)
			{
				total_match_length1 = total_match_length1 + match_step;
			}
			else
			{
				total_match_length1 = total_match_length1 + match_length1 / 2;
			}



			seed_id1++;


		}
	}


	///这里要改
	///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
	if (extra_seed_flag1 == 0
		&&
		(candidate_length1 == 1 ||
		(candidate_length1 == 2 && candidates1[0] == candidates1[1])
		))
	{

		candidates_votes1[0].site = candidates1[0];
		candidates_votes1[0].err = 1;
		candidates_votes1[0].end_site = (*current_read).length - 1;
		best_mapp_occ1 = 1;
	}
	else if (candidate_length1 != 0)
	{

		std::sort(candidates1, candidates1 + candidate_length1);

		generate_candidate_votes_shift(candidates1, candidate_length1, candidates_votes1, &candidates_votes_length1, error_threshold);


		if (error_threshold <= 15)
		{
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}
		else
		{
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
		}



	}
	else
	{
		///没有种子
		best_mapp_occ1 = 0;
		///return_flag1 = 0;
	}



	return best_mapp_occ1;







}





inline int verify_candidate_locations(
	bitmapper_bs_iter error_threshold,
	bitmapper_bs_iter candidates_votes_length,
	Read* current_read,
	char* bsSeq,
	char* cigar,
	char* tmp_ref,
	seed_votes *candidates_votes,
	__m256i* Peq_SSE)
{
	unsigned int min_err;
	int best_mapp_occ = 0;

	if (error_threshold <= 15)
	{
		map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes, &candidates_votes_length, (*current_read).length,
			error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ, (*current_read).name);
	}
	else
	{
		map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes, &candidates_votes_length, (*current_read).length,
			error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ, (*current_read).name);
	}

	return best_mapp_occ;
}


inline void get_candidates(
	int max_seed_number, 
	int max_candidates_occ,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	bitmapper_bs_iter error_threshold,
	bitmapper_bs_iter candidates_votes_length,
	Read* current_read,
	char* bsSeq,
	char* cigar,
	char* tmp_ref,
	bitmapper_bs_iter* candidates,
	seed_votes *candidates_votes,
	int* result_best_mapp_occ,
	bitmapper_bs_iter* result_candidates_votes_length)
{
	long long unique_best_read = -1;

	int direct_jump = 0;

	int best_mapp_occ = -1;

	int C_site;


	C_to_T_forward(current_read->seq, bsSeq, current_read->length, &C_site);

	bitmapper_bs_iter total_match_length = 0;

	bitmapper_bs_iter seed_id = 0;

	bitmapper_bs_iter* locates = candidates;

	bitmapper_bs_iter candidate_length = 0;

	int is_mutiple_map = 0;

	int read_jump = 0;

	///这里要改
	int get_error = -1;
	///这里要改
	int extra_seed_flag = 1;






	bitmapper_bs_iter current_seed_length;
	bitmapper_bs_iter number_of_hits;

	bitmapper_bs_iter top, bot;
	bitmapper_bs_iter pre_top, pre_bot;
	bitmapper_bs_iter match_length;
	bitmapper_bs_iter tmp_SA_length;
	long long first_seed_match_length;
	int return_flag;
	


	////5'端第一个seed单独处理
	if (seed_id < max_seed_number && total_match_length < current_read->length)
	{
		current_seed_length = current_read->length - total_match_length;


		number_of_hits =
			count_backward_as_much_1_terminate
			(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);

		///这里要改
		first_seed_match_length = match_length;


		if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end_return_site(
			current_read->length, current_read->seq, current_read->rseq, current_read->qual, cigar,
			current_read->name,
			locates,
			top, tmp_ref, &match_length, current_read->length - C_site - 1,
			&candidates_votes[0].site, &candidates_votes[0].err, &get_error))
		{
			candidates_votes[0].end_site = current_read->length - 1;
			best_mapp_occ = 1;
			///唯一匹配
			return_flag = 1;
			direct_jump = 1;
			read_jump = 1;
			goto read1_end;
		}


		///只有第一个seed才有可能出现这种情况
		if (match_length == current_read->length
			&&
			number_of_hits > 1
			&&
			number_of_hits <= max_candidates_occ)
		{
			is_mutiple_map = 1;

			///等于-1说明read里面没有C
			if (C_site == -1)
			{
				locate(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length);
				reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

				///排序是为了最后好合并双端的结果
				std::sort(locates, locates + tmp_SA_length);

				best_mapp_occ = tmp_SA_length;

				for (int inner_i = 0; inner_i < tmp_SA_length; inner_i++)
				{
					candidates_votes[inner_i].site = locates[inner_i];
					candidates_votes[inner_i].err = 0;
					candidates_votes[inner_i].end_site = current_read->length - 1;
				}




				///多匹配
				return_flag = 2;
				direct_jump = 1;
				read_jump = 1;
				goto read1_end;
			}
		}




		///如果当前seed的只有一个匹配位置
		///则不管多长都是有效的
		///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
		if (number_of_hits == 1)
		{

			tmp_SA_length = 1;

			candidate_length = candidate_length + tmp_SA_length;

			locates = candidates + candidate_length;


		}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
		else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
		{

			if (number_of_hits != 0)
			{

				tmp_SA_length = 0;

				locate(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length);
				reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

				candidate_length = candidate_length + tmp_SA_length;
				locates = candidates + candidate_length;
			}

		}



		if (match_length == 0)
		{
			total_match_length = total_match_length + match_step;
		}
		else
		{
			total_match_length = total_match_length + match_length / 2;
		}


		seed_id++;

	}
read1_end:



	///如果第一个种子没有找到匹配位置
	///那么就得接着取剩下的种子
	if (read_jump == 0)
	{
		long long second_seed_length;
		int extra_seed_flag = 1;

		///这代表第一个种子已经找到了一个误配数目为1的匹配
		if (get_error == 1)
		{

			second_seed_length = current_read->length - first_seed_match_length;


			if (second_seed_length >= 17)
			{

				number_of_hits = count_hash_table(bsSeq, second_seed_length, &top, &bot, &pre_top, &pre_bot);

				if (number_of_hits == 1)
				{
					tmp_SA_length = 0;

					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, second_seed_length, first_seed_match_length);

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;

					extra_seed_flag = 0;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (number_of_hits <= max_seed_matches)
				{

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate(bsSeq, top, bot, pre_top, pre_bot, locates, second_seed_length, &tmp_SA_length);

						reverse_and_adjust_site(locates, &tmp_SA_length, second_seed_length, first_seed_match_length);


						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}

					extra_seed_flag = 0;

				}
				else if (number_of_hits > max_seed_matches)
				{
					extra_seed_flag = 1;
				}


			}
			else
			{
				extra_seed_flag = 1;
			}

		}

		///为1有两种情况
		///1. 第一个种子没有找到一个误配为1的匹配
		///2. 尽管第一个种子找到了一个误配为1的匹配，但第二个种子频率太高或者太短了
		///这样的话还是要抽更多的种子出来
		if (extra_seed_flag == 1)
		{


			///这是剩下的seed
			while (seed_id < max_seed_number && total_match_length < current_read->length)
			{

				current_seed_length = current_read->length - total_match_length;


				number_of_hits =
					count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);

				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				if (number_of_hits == 1)
				{
					tmp_SA_length = 0;

					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, match_length, total_match_length);

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
				{

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length);
						reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}

				}
				else if (current_seed_length == match_length)
				{
					break;
				}






				if (match_length == 0)
				{
					total_match_length = total_match_length + match_step;
				}
				else
				{
					total_match_length = total_match_length + match_length / 2;
				}



				seed_id++;


			}
		}




		///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
		///extra_seed_flag == 0表明已经找到了一个1误配的匹配位置，并且执行了完整的鸽巢原理
		///在这种情况下，如果仍然只有一个匹配位置，那么就找到了best mapping
		if (extra_seed_flag == 0
			&&
			(candidate_length == 1 ||
			(candidate_length == 2 && candidates[0] == candidates[1])
			))
		{

			candidates_votes[0].site = candidates[0];
			candidates_votes[0].err = 1;
			candidates_votes[0].end_site = current_read->length - 1;
			best_mapp_occ = 1;
		}
		else if (candidate_length != 0)
		{

			std::sort(candidates, candidates + candidate_length);

			generate_candidate_votes_shift(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold);

		}
		else
		{
			///没有种子
			best_mapp_occ = 0;
			///return_flag1 = 0;
		}



	}


	/**
	当程序执行到这里时，有三种情况
	1.
	best_mapp_occ = -1，这说明还有候选位置没有验证，
	此时candidates_votes都是候选位置，candidates_votes_length是候选位置数目
	2.
	best_mapp_occ>0，此时说明已经找到了匹配位置
	3. 
	best_mapp_occ = 0，此时说明不匹配

	candidates_votes中的位置都是排过序的
	**/
	if (best_mapp_occ == 0)
	{
		candidates_votes_length = 0;
	}
	else if (best_mapp_occ > 0)
	{
		candidates_votes_length = best_mapp_occ;
	}


	(*result_best_mapp_occ) = best_mapp_occ;
	(*result_candidates_votes_length) = candidates_votes_length;

}




int Map_Pair_Seq_end_to_end_fast(int thread_id)
{
	long long enq_i = 0;
	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;


	char** ref = (char**)malloc(sizeof(char*)* 8);

	for (i = 0; i<8; i++)
	{
		ref[i] = (char*)malloc(sizeof(char)*(SEQ_LENGTH + 2 * thread_e));
	}


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;

	Read current_read1, current_read2;
	init_single_read(&current_read1);
	init_single_read(&current_read2);

	//bitmapper_bs_iter number_of_hits = 0;
	bitmapper_bs_iter number_of_hits1, number_of_hits2;

	unsigned int bsSeq_max_length = 1000;
	///char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	///bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter top, bot, match_step;
	bitmapper_bs_iter total_match_length1, total_match_length2;
	bitmapper_bs_iter match_length1, match_length2;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;
	///bitmapper_bs_iter available_seed_length = 30;

	///min_seed_length这个值不能改，改了就会报错
	///bitmapper_bs_iter min_seed_length = 17;
	///bitmapper_bs_iter min_seed_max_seed_matches = 100;

	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	bitmapper_bs_iter short_seed_max_seed_matches = 20;

	available_seed_length = over_all_seed_length;

	///fprintf(stderr, "over_all_seed_length: %d\n", over_all_seed_length);

	/**
	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	**/

	//bitmapper_bs_iter current_seed_length;
	bitmapper_bs_iter current_seed_length1, current_seed_length2;

	///bitmapper_bs_iter* locates;
	bitmapper_bs_iter *locates1, *locates2;
	///bitmapper_bs_iter* candidates;
	bitmapper_bs_iter *candidates1, *candidates2;
	///bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter candidate_length1, candidate_length2;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	///seed_votes* candidates_votes;
	seed_votes *candidates_votes1, *candidates_votes2;
	///bitmapper_bs_iter candidates_votes_length;
	bitmapper_bs_iter candidates_votes_length1, candidates_votes_length2;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	//bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter max_seed_number1, max_seed_number2;
	///bitmapper_bs_iter seed_id;
	bitmapper_bs_iter seed_id1, seed_id2;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter ambious_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	///bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;
	///int is_duplicate;
	///bitmapper_bs_iter min_index;


	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;
	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
	error_threshold = thread_e;
	}
	**/

	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	///max_seed_number = 9;
	max_seed_number1 = 25;
	max_seed_number2 = 25;
	int max_max_seed_number1 = max_seed_number1;
	int max_max_seed_number2 = max_seed_number2;



	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(large_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	double re_seed_threshold = 0.80;
	int reseed_length;
	bitmapper_bs_iter re_seed_error_threshold;

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	long long total_number_of_hit = 0;


	fprintf(stdout, "Welcome to BitMapperBS!\n");


	int is_mutiple_map1, is_mutiple_map2;

	long long distance;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	long long direct_cut = 0;



	int add_empty_read[4] = { 1, 0, 0, 0 };
	int add_matched_read[4] = { 0, 1, 1, 0 };
	int add_unique_matched_read[4] = { 0, 1, 0, 0 };
	int add_unmatched_read[4] = { 0, 0, 0, 1 };
	int best_mapp_occ1 = 0;
	int best_mapp_occ2 = 0;
	int inner_i, inner_j;
	long long mapping_pair;
	long long best_pair_1_index;
	long long best_pair_2_index;
	///long long best_sum_err = 2 * error_threshold + 1;
	long long best_sum_err;
	long long current_sum_err;

	long long unique_best_read1;
	long long unique_best_read2;
	map_result result1, result2;


	int max_candidates_occ = max_seed_number1 * max_seed_matches;


	candidates1 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates2 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates_votes1 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));
	candidates_votes2 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));

	seed_votes *tmp_votes1, *tmp_votes2;
	tmp_votes1 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));
	tmp_votes2 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));



	///这是个特殊情况，仅在很偶尔的地方才有用
	///就是一边read有超过max_candidates_occ个误配的时候，就不要了
	if (max_candidates_occ > 10000)
	{
		max_candidates_occ = 10000;
	}


	int matched_length1, matched_length2;



	char* bsSeq1 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	char* bsSeq2 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	int C_site1;
	int C_site2;
	int direct_jump1 = 0;
	int direct_jump2 = 0;
	long long current_best_err;

	int first_seed_read1_occ, first_seed_read2_occ;
	int first_seed_read1_length, first_seed_read2_length;



	long long debug_1 = 0, debug_2 = 0;


	fprintf(stdout, "max_candidates_occ: %lld\n", max_candidates_occ);


	double startTime = Get_T();

	int file_flag;
	int return_flag1, return_flag2;
	int read_1_jump, read_2_jump;

	int inner_maxDistance_pair;
	int inner_minDistance_pair;



	int get_error1, get_error2;
	long long first_seed_match_length1, first_seed_match_length2;
	long long second_seed_length1, second_seed_length2;
	int extra_seed_flag1, extra_seed_flag2;
	int paired_end_distance;

	//正向模式
	i = 0;
	while (1)
	{


		file_flag = inputReads_paired_directly(&current_read1, &current_read2);



		///如果等于0, 估计说明read读完了吧
		if (file_flag == 0)
		{
			break;
		}


		enq_i++;


		///等于3应该是有过多的N吧
		if (file_flag == 3)
		{
			continue;
		}




		error_threshold1 = thread_e_f * current_read1.length;
		if (error_threshold1 >= max_error_threshold)
		{
			error_threshold1 = max_error_threshold;
		}


		error_threshold2 = thread_e_f * current_read2.length;
		if (error_threshold2 >= max_error_threshold)
		{
			error_threshold2 = max_error_threshold;
		}

		if (error_threshold1 > error_threshold2)
		{
			large_error_threshold = error_threshold1;
		}
		else
		{
			large_error_threshold = error_threshold2;
		}

		best_sum_err = 2 * large_error_threshold + 1;



		inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
		inner_minDistance_pair = minDistance_pair - large_error_threshold * 2;



		read_1_jump = 0;
		read_2_jump = 0;


		//*************************************第一个read的第一个种子********************************************

		//*******控制种子数量************************
		max_seed_number1 = current_read1.length / 10 - 1;
		if (max_seed_number1 > max_max_seed_number1)
		{
			max_seed_number1 = max_max_seed_number1;
		}
		//*******控制种子数量************************

		/*************************************第二个read的第一个种子********************************************/

		/*******控制种子数量************************/
		max_seed_number2 = current_read2.length / 10 - 1;
		if (max_seed_number2 > max_max_seed_number2)
		{
			max_seed_number2 = max_max_seed_number2;
		}
		/*******控制种子数量************************/






		get_candidates(
			max_seed_number1, 
			max_candidates_occ,
			available_seed_length,
			max_seed_matches,
			match_step,
			error_threshold1,
			candidates_votes_length1,
			&current_read1,
			bsSeq1,
			cigar,
			tmp_ref,
			candidates1,
			candidates_votes1,
			&best_mapp_occ1,
			&candidates_votes_length1);


		get_candidates(
			max_seed_number2,
			max_candidates_occ,
			available_seed_length,
			max_seed_matches,
			match_step,
			error_threshold2,
			candidates_votes_length2,
			&current_read2,
			bsSeq2,
			cigar,
			tmp_ref,
			candidates2,
			candidates_votes2,
			&best_mapp_occ2,
			&candidates_votes_length2);













		if (best_mapp_occ1 > 0 && best_mapp_occ2 > 0)   ///说明两个都不要验证
		{
			best_mapp_occ1 = candidates_votes_length1;
			best_mapp_occ2 = candidates_votes_length2;
		}
		else
		{

			if (best_mapp_occ1 == 0 || best_mapp_occ2 == 0)
			{
				i++;
				continue;
			}


			filter_pairs(
				candidates_votes_length1,
				candidates_votes_length2,
				&candidates_votes1,
				&candidates_votes2,
				&tmp_votes1,
				&tmp_votes2,
				&candidates_votes_length1,
				&candidates_votes_length2,
				inner_maxDistance_pair,
				inner_minDistance_pair
				);

			if (candidates_votes_length1 == 0 || candidates_votes_length2 == 0)
			{
				i++;
				continue;
			}

			/**
			当程序执行到这里时，有三种情况
			1.
			best_mapp_occ = -1，这说明还有候选位置没有验证，
			此时candidates_votes都是候选位置，candidates_votes_length是候选位置数目
			2.
			best_mapp_occ>0，此时说明已经找到了匹配位置
			3.
			best_mapp_occ = 0，此时说明不匹配

			candidates_votes中的位置都是排过序的
			**/
			if (best_mapp_occ1 == -1 && best_mapp_occ2 == -1)///这是两个read都需要验证的情况
			{



				if (candidates_votes_length1 <= candidates_votes_length2)
				{
					best_mapp_occ1 = verify_candidate_locations(
						error_threshold1,
						candidates_votes_length1,
						&current_read1,
						bsSeq1,
						cigar,
						tmp_ref,
						candidates_votes1,
						Peq_SSE);


					if (best_mapp_occ1 == 0)
					{
						i++;
						continue;
					}

					///用1的结果过滤2
					///所以1用的是best_mapp_occ1, 而2用的是candidates_votes_length2
					filter_pairs_single_side(
						best_mapp_occ1,
						candidates_votes_length2,
						&candidates_votes1,
						&candidates_votes2,
						&candidates_votes_length2,
						inner_maxDistance_pair,
						inner_minDistance_pair
						);

					best_mapp_occ2 = verify_candidate_locations(
						error_threshold2,
						candidates_votes_length2,
						&current_read2,
						bsSeq2,
						cigar,
						tmp_ref,
						candidates_votes2,
						Peq_SSE);
				}
				else
				{

					best_mapp_occ2 = verify_candidate_locations(
						error_threshold2,
						candidates_votes_length2,
						&current_read2,
						bsSeq2,
						cigar,
						tmp_ref,
						candidates_votes2,
						Peq_SSE);


					if (best_mapp_occ2 == 0)
					{
						i++;
						continue;
					}



					///用2的结果过滤1
					///所以2用的是best_mapp_occ2, 而1用的是candidates_votes_length1
					filter_pairs_single_side(
						best_mapp_occ2,
						candidates_votes_length1,
						&candidates_votes2,
						&candidates_votes1,
						&candidates_votes_length1,
						inner_maxDistance_pair,
						inner_minDistance_pair
						);


					best_mapp_occ1 = verify_candidate_locations(
						error_threshold1,
						candidates_votes_length1,
						&current_read1,
						bsSeq1,
						cigar,
						tmp_ref,
						candidates_votes1,
						Peq_SSE);
				}

			}
			else if (best_mapp_occ1 != -1)
			{

				if (candidates_votes_length1 < best_mapp_occ1)
				{
					best_mapp_occ1 = candidates_votes_length1;
				}

				if (best_mapp_occ2 == -1)
				{
					best_mapp_occ2 = verify_candidate_locations(
						error_threshold2,
						candidates_votes_length2,
						&current_read2,
						bsSeq2,
						cigar,
						tmp_ref,
						candidates_votes2,
						Peq_SSE);
				}
			}
			else if (best_mapp_occ2 != -1)
			{
				if (candidates_votes_length2 < best_mapp_occ2)
				{
					best_mapp_occ2 = candidates_votes_length2;
				}



				if (best_mapp_occ1 == -1)
				{
					best_mapp_occ1 = verify_candidate_locations(
						error_threshold1,
						candidates_votes_length1,
						&current_read1,
						bsSeq1,
						cigar,
						tmp_ref,
						candidates_votes1,
						Peq_SSE);
				}
			}

		}





















		mapping_pair = 0;


		mapping_pair = new_faster_verify_pairs(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
				candidates_votes1, candidates_votes2,
				&best_pair_1_index, &best_pair_2_index,
				inner_maxDistance_pair, inner_minDistance_pair);
		


	end_i:
		if (mapping_pair == 1)
		{


			//****************第一个read的后处理
			result1.err = candidates_votes1[best_pair_1_index].err;
			result1.origin_site = candidates_votes1[best_pair_1_index].site;
			result1.end_site = candidates_votes1[best_pair_1_index].end_site;

			if (result1.err != 0)
			{
				calculate_best_map_cigar_end_to_end_return
					(&min_candidates_votes_length,
					current_read1.length, error_threshold1, tmp_ref,
					current_read1.seq, current_read1.rseq, current_read1.qual, cigar, path, matrix, matrix_bit,
					current_read1.name,
					result1.cigar,
					&result1, &matched_length1);



			}
			else
			{

				output_sam_end_to_end_return(
					result1.origin_site,
					result1.end_site,
					result1.end_site + 1 - current_read1.length,
					current_read1.length,
					&result1
					);

				sprintf(result1.cigar, "%dM", current_read1.length);
				matched_length1 = current_read1.length;
			}
			//****************第一个read的后处理




			//****************第二个read的后处理
			result2.err = candidates_votes2[best_pair_2_index].err;
			result2.origin_site = candidates_votes2[best_pair_2_index].site;
			result2.end_site = candidates_votes2[best_pair_2_index].end_site;

			if (result2.err != 0)
			{
				calculate_best_map_cigar_end_to_end_return
					(&min_candidates_votes_length,
					current_read2.length, error_threshold2, tmp_ref,
					current_read2.seq, current_read2.rseq, current_read2.qual, cigar, path, matrix, matrix_bit,
					current_read2.name,
					result2.cigar,
					&result2, &matched_length2);
			}
			else
			{

				output_sam_end_to_end_return(
					result2.origin_site,
					result2.end_site,
					result2.end_site + 1 - current_read2.length,
					current_read2.length,
					&result2
					);

				sprintf(result2.cigar, "%dM", current_read2.length);
				matched_length2 = current_read2.length;
			}
			//****************第二个read的后处理



			if (result2.site >= result1.site)
			{
				paired_end_distance = result2.site - result1.site + matched_length2;

			}
			else
			{
				paired_end_distance = result1.site - result2.site + matched_length1;
			}

			if (paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
			{
				unique_matched_read++;


				directly_output_read1(current_read1.name, current_read1.seq, current_read1.rseq, current_read1.qual,
					&result1, &result2, current_read1.length,
					matched_length1, matched_length2);

				directly_output_read2(current_read2.name, current_read2.seq, current_read2.rseq, current_read2.qual,
					&result2, &result1, current_read2.length,
					matched_length2, matched_length1);
			}





		}
		else if (mapping_pair > 1)
		{
			ambious_matched_read++;
		}

	}


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;


	fprintf(stderr, "debug_1: %lld\n", debug_1);
	fprintf(stderr, "debug_2: %lld\n", debug_2);


	return 1;
}



inline void get_candidates_muti_thread(
	int max_seed_number,
	int max_candidates_occ,
	bitmapper_bs_iter available_seed_length,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter match_step,
	bitmapper_bs_iter error_threshold,
	bitmapper_bs_iter candidates_votes_length,
	Read* current_read,
	char* bsSeq,
	char* cigar,
	char* tmp_ref,
	bitmapper_bs_iter* candidates,
	seed_votes *candidates_votes,
	int* result_best_mapp_occ,
	bitmapper_bs_iter* result_candidates_votes_length,
	bwt_locate_queue* get_queue)
{
	long long unique_best_read = -1;

	int direct_jump = 0;

	int best_mapp_occ = -1;

	int C_site;


	C_to_T_forward(current_read->seq, bsSeq, current_read->length, &C_site);

	bitmapper_bs_iter total_match_length = 0;

	bitmapper_bs_iter seed_id = 0;

	bitmapper_bs_iter* locates = candidates;

	bitmapper_bs_iter candidate_length = 0;

	int is_mutiple_map = 0;

	int read_jump = 0;

	///这里要改
	int get_error = -1;
	///这里要改
	int extra_seed_flag = 1;






	bitmapper_bs_iter current_seed_length;
	bitmapper_bs_iter number_of_hits;

	bitmapper_bs_iter top, bot;
	bitmapper_bs_iter pre_top, pre_bot;
	bitmapper_bs_iter match_length;
	bitmapper_bs_iter tmp_SA_length;
	long long first_seed_match_length;
	int return_flag;



	////5'端第一个seed单独处理
	if (seed_id < max_seed_number && total_match_length < current_read->length)
	{
		current_seed_length = current_read->length - total_match_length;


		number_of_hits =
			count_backward_as_much_1_terminate
			(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);

		///这里要改
		first_seed_match_length = match_length;


		if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end_return_site(
			current_read->length, current_read->seq, current_read->rseq, current_read->qual, cigar,
			current_read->name,
			locates,
			top, tmp_ref, &match_length, current_read->length - C_site - 1,
			&candidates_votes[0].site, &candidates_votes[0].err, &get_error))
		{
			candidates_votes[0].end_site = current_read->length - 1;
			best_mapp_occ = 1;
			///唯一匹配
			return_flag = 1;
			direct_jump = 1;
			read_jump = 1;
			goto read1_end;
		}


		///只有第一个seed才有可能出现这种情况
		if (match_length == current_read->length
			&&
			number_of_hits > 1
			&&
			number_of_hits <= max_candidates_occ)
		{
			is_mutiple_map = 1;

			///等于-1说明read里面没有C
			if (C_site == -1)
			{
				locate_muti_thread(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length, get_queue);
				reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

				///排序是为了最后好合并双端的结果
				std::sort(locates, locates + tmp_SA_length);

				best_mapp_occ = tmp_SA_length;

				for (int inner_i = 0; inner_i < tmp_SA_length; inner_i++)
				{
					candidates_votes[inner_i].site = locates[inner_i];
					candidates_votes[inner_i].err = 0;
					candidates_votes[inner_i].end_site = current_read->length - 1;
				}




				///多匹配
				return_flag = 2;
				direct_jump = 1;
				read_jump = 1;
				goto read1_end;
			}
		}




		///如果当前seed的只有一个匹配位置
		///则不管多长都是有效的
		///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
		if (number_of_hits == 1)
		{

			tmp_SA_length = 1;

			candidate_length = candidate_length + tmp_SA_length;

			locates = candidates + candidate_length;


		}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
		else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
		{

			if (number_of_hits != 0)
			{

				tmp_SA_length = 0;

				locate_muti_thread(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length, get_queue);
				reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

				candidate_length = candidate_length + tmp_SA_length;
				locates = candidates + candidate_length;
			}

		}



		if (match_length == 0)
		{
			total_match_length = total_match_length + match_step;
		}
		else
		{
			total_match_length = total_match_length + match_length / 2;
		}


		seed_id++;

	}
read1_end:



	///如果第一个种子没有找到匹配位置
	///那么就得接着取剩下的种子
	if (read_jump == 0)
	{
		long long second_seed_length;
		int extra_seed_flag = 1;

		///这代表第一个种子已经找到了一个误配数目为1的匹配
		if (get_error == 1)
		{

			second_seed_length = current_read->length - first_seed_match_length;


			if (second_seed_length >= 17)
			{

				number_of_hits = count_hash_table(bsSeq, second_seed_length, &top, &bot, &pre_top, &pre_bot);

				if (number_of_hits == 1)
				{
					tmp_SA_length = 0;

					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, second_seed_length, first_seed_match_length);

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;

					extra_seed_flag = 0;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (number_of_hits <= max_seed_matches)
				{

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate_muti_thread(bsSeq, top, bot, pre_top, pre_bot, locates, second_seed_length, &tmp_SA_length, get_queue);

						reverse_and_adjust_site(locates, &tmp_SA_length, second_seed_length, first_seed_match_length);


						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}

					extra_seed_flag = 0;

				}
				else if (number_of_hits > max_seed_matches)
				{
					extra_seed_flag = 1;
				}


			}
			else
			{
				extra_seed_flag = 1;
			}

		}

		///为1有两种情况
		///1. 第一个种子没有找到一个误配为1的匹配
		///2. 尽管第一个种子找到了一个误配为1的匹配，但第二个种子频率太高或者太短了
		///这样的话还是要抽更多的种子出来
		if (extra_seed_flag == 1)
		{


			///这是剩下的seed
			while (seed_id < max_seed_number && total_match_length < current_read->length)
			{

				current_seed_length = current_read->length - total_match_length;


				number_of_hits =
					count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);

				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				if (number_of_hits == 1)
				{
					tmp_SA_length = 0;

					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, match_length, total_match_length);

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
				{

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate_muti_thread(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length, get_queue);
						reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}

				}
				else if (current_seed_length == match_length)
				{
					break;
				}






				if (match_length == 0)
				{
					total_match_length = total_match_length + match_step;
				}
				else
				{
					total_match_length = total_match_length + match_length / 2;
				}



				seed_id++;


			}
		}




		///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
		///extra_seed_flag == 0表明已经找到了一个1误配的匹配位置，并且执行了完整的鸽巢原理
		///在这种情况下，如果仍然只有一个匹配位置，那么就找到了best mapping
		if (extra_seed_flag == 0
			&&
			(candidate_length == 1 ||
			(candidate_length == 2 && candidates[0] == candidates[1])
			))
		{

			candidates_votes[0].site = candidates[0];
			candidates_votes[0].err = 1;
			candidates_votes[0].end_site = current_read->length - 1;
			best_mapp_occ = 1;
		}
		else if (candidate_length != 0)
		{

			std::sort(candidates, candidates + candidate_length);

			generate_candidate_votes_shift(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold);

		}
		else
		{
			///没有种子
			best_mapp_occ = 0;
			///return_flag1 = 0;
		}



	}


	/**
	当程序执行到这里时，有三种情况
	1.
	best_mapp_occ = -1，这说明还有候选位置没有验证，
	此时candidates_votes都是候选位置，candidates_votes_length是候选位置数目
	2.
	best_mapp_occ>0，此时说明已经找到了匹配位置
	3.
	best_mapp_occ = 0，此时说明不匹配

	candidates_votes中的位置都是排过序的
	**/
	if (best_mapp_occ == 0)
	{
		candidates_votes_length = 0;
	}
	else if (best_mapp_occ > 0)
	{
		candidates_votes_length = best_mapp_occ;
	}


	(*result_best_mapp_occ) = best_mapp_occ;
	(*result_candidates_votes_length) = candidates_votes_length;

}



///输出成对匹配
///和一端unique一端不匹配的结果
int Map_Pair_Seq_end_to_end(int thread_id)
{
	long long enq_i = 0;
	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;


	char** ref = (char**)malloc(sizeof(char*)* 8);

	for (i = 0; i<8; i++)
	{
		ref[i] = (char*)malloc(sizeof(char)*(SEQ_LENGTH + 2 * thread_e));
	}


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;

	Read current_read1, current_read2;
	init_single_read(&current_read1);
	init_single_read(&current_read2);

	//bitmapper_bs_iter number_of_hits = 0;
	bitmapper_bs_iter number_of_hits1, number_of_hits2;

	unsigned int bsSeq_max_length = 1000;
	///char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	///bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter top, bot, match_step;
	bitmapper_bs_iter total_match_length1, total_match_length2;
	bitmapper_bs_iter match_length1, match_length2;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;
	///bitmapper_bs_iter available_seed_length = 30;

	///min_seed_length这个值不能改，改了就会报错
	///bitmapper_bs_iter min_seed_length = 17;
	///bitmapper_bs_iter min_seed_max_seed_matches = 100;

	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	bitmapper_bs_iter short_seed_max_seed_matches = 20;


	available_seed_length = over_all_seed_length;

	////fprintf(stderr, "over_all_seed_length: %d\n", over_all_seed_length);

	/**
	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	**/

	//bitmapper_bs_iter current_seed_length;
	bitmapper_bs_iter current_seed_length1, current_seed_length2;

	///bitmapper_bs_iter* locates;
	bitmapper_bs_iter *locates1, *locates2;
	///bitmapper_bs_iter* candidates;
	bitmapper_bs_iter *candidates1, *candidates2;
	///bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter candidate_length1, candidate_length2;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	///seed_votes* candidates_votes;
	seed_votes *candidates_votes1, *candidates_votes2;
	///bitmapper_bs_iter candidates_votes_length;
	bitmapper_bs_iter candidates_votes_length1, candidates_votes_length2;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	//bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter max_seed_number1, max_seed_number2;
	///bitmapper_bs_iter seed_id;
	bitmapper_bs_iter seed_id1, seed_id2;
	bitmapper_bs_iter full_seed_id1, full_seed_id2;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter ambious_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	///bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;
	///int is_duplicate;
	///bitmapper_bs_iter min_index;



	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;
	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
		error_threshold = thread_e;
	}
	**/


	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	///max_seed_number = 9;
	max_seed_number1 = 25;
	max_seed_number2 = 25;
	int max_max_seed_number1 = max_seed_number1;
	int max_max_seed_number2 = max_seed_number2;



	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(max_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	double re_seed_threshold = 0.80;
	int reseed_length;
	bitmapper_bs_iter re_seed_error_threshold;

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	long long total_number_of_hit = 0;


	fprintf(stdout, "Welcome to BitMapperBS!\n");


	int is_mutiple_map1, is_mutiple_map2;

	long long distance;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	long long direct_cut = 0;



	int add_empty_read[4] = { 1, 0, 0, 0 };
	int add_matched_read[4] = { 0, 1, 1, 0 };
	int add_unique_matched_read[4] = { 0, 1, 0, 0 };
	int add_unmatched_read[4] = { 0, 0, 0, 1 };
	int best_mapp_occ1 = 0;
	int best_mapp_occ2 = 0;
	int inner_i, inner_j;
	long long mapping_pair;
	long long best_pair_1_index;
	long long best_pair_2_index;
	///long long best_sum_err = 2 * error_threshold + 1;
	long long best_sum_err;
	long long current_sum_err;

	long long unique_best_read1;
	long long unique_best_read2;
	map_result result1, result2;


	int max_candidates_occ = max_seed_number1 * max_seed_matches * 2.5;


	candidates1 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates2 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates_votes1 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));
	candidates_votes2 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));

	///这是个特殊情况，仅在很偶尔的地方才有用
	///就是一边read有超过max_candidates_occ个误配的时候，就不要了
	if (max_candidates_occ > 10000)
	{
		max_candidates_occ = 10000;
	}


	int matched_length1, matched_length2;



	char* bsSeq1 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	char* bsSeq2 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	int C_site1;
	int C_site2;
	int direct_jump1 = 0;
	int direct_jump2 = 0;
	long long current_best_err;

	int first_seed_read1_occ, first_seed_read2_occ;
	int first_seed_read1_length, first_seed_read2_length;



	long long debug_1 = 0, debug_2 = 0;


	fprintf(stdout, "max_candidates_occ: %lld\n", max_candidates_occ);


	double startTime = Get_T();

	int file_flag;
	int return_flag1, return_flag2;
	int read_1_jump, read_2_jump;

	int inner_maxDistance_pair;
	int inner_minDistance_pair;



	int get_error1, get_error2;
	long long first_seed_match_length1, first_seed_match_length2;
	long long second_seed_length1, second_seed_length2;
	int extra_seed_flag1, extra_seed_flag2;
	int paired_end_distance;
	int read1_max_last_zone_start = -1;
	int read1_max_last_zone_length = -1;
	int read2_max_last_zone_start = -1;
	int read2_max_last_zone_length = -1;



	int* read1_seed_start = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* read1_seed_length = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* read2_seed_start = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* read2_seed_length = (int*)malloc(sizeof(int)* max_max_seed_number1);

	int* result_seed_start = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* result_seed_length = (int*)malloc(sizeof(int)* max_max_seed_number1);

	//正向模式
	i = 0;
	while (1)
	{


		file_flag = inputReads_paired_directly(&current_read1, &current_read2);



		///如果等于0, 估计说明read读完了吧
		if (file_flag == 0)
		{
			break;
		}


		enq_i++;


		///等于3应该是有过多的N吧
		if (file_flag == 3)
		{
			continue;
		}




		error_threshold1 = thread_e_f * current_read1.length;
		if (error_threshold1 >= max_error_threshold)
		{
			error_threshold1 = max_error_threshold;
		}


		error_threshold2 = thread_e_f * current_read2.length;
		if (error_threshold2 >= max_error_threshold)
		{
			error_threshold2 = max_error_threshold;
		}

		if (error_threshold1 > error_threshold2)
		{
			large_error_threshold = error_threshold1;
		}
		else
		{
			large_error_threshold = error_threshold2;
		}

		best_sum_err = 2 * large_error_threshold + 1;



		inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
		inner_minDistance_pair = minDistance_pair - large_error_threshold * 2;




		/*************************************处理第一个种子********************************************/

		read_1_jump = 0;
		read_2_jump = 0;


		//*************************************第一个read的第一个种子********************************************




		
		//*******控制种子数量************************
		max_seed_number1 = current_read1.length / 10 - 1;
		if (max_seed_number1 > max_max_seed_number1)
		{
			max_seed_number1 = max_max_seed_number1;
		}
		//*******控制种子数量************************




		unique_best_read1 = -1;

		direct_jump1 = 0;

		best_mapp_occ1 = 0;


		C_to_T_forward(current_read1.seq, bsSeq1, current_read1.length, &C_site1);

		total_match_length1 = 0;

		seed_id1 = 0;

		full_seed_id1 = 0;


		locates1 = candidates1;

		candidate_length1 = 0;

		is_mutiple_map1 = 0;


		///这里要改
		get_error1 = -1;
		///这里要改
		extra_seed_flag1 = 1;


		////5'端第一个seed单独处理
		if (seed_id1 < max_seed_number1 && total_match_length1 < current_read1.length)
		{
			current_seed_length1 = current_read1.length - total_match_length1;


			number_of_hits1 =
				count_backward_as_much_1_terminate
				(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);

			///这里要改
			first_seed_match_length1 = match_length1;

			/**
			if (strcmp("@SRR1532535.1001096", current_read1.name) == 0)
			{
				fprintf(stderr, "\n\n%s(1)\n", current_read1.name);
				fprintf(stderr, "number_of_hits1: %d\n", number_of_hits1);
				fprintf(stderr, "match_length1: %d\n", match_length1);
				fprintf(stderr, "current_read1.length: %d\n", current_read1.length);

				
			}
			**/
			

			/**
			有两种情况能直接跳转
			1. 找到一个精确的唯一匹配
			2. 找到一堆精确的唯一匹配
			这两种情况best_mapp_occ1都不是0
			所以reseed不用担心这两种情况
			**/
			////新加的
			read1_seed_start[full_seed_id1] = total_match_length1;
			read1_seed_length[full_seed_id1] = match_length1;
			full_seed_id1++;
			////新加的




			///满足这个条件就直接匹配
			if (number_of_hits1 == 1 && try_process_unique_mismatch_end_to_end_return_site(
				current_read1.length, current_read1.seq, current_read1.rseq, current_read1.qual, cigar,
				current_read1.name,
				locates1,
				top, tmp_ref, &match_length1, current_read1.length - C_site1 - 1,
				&candidates_votes1[0].site, &candidates_votes1[0].err, &get_error1))
			{

				candidates_votes1[0].end_site = current_read1.length - 1;
				best_mapp_occ1 = 1;
				///唯一匹配
				return_flag1 = 1;
				direct_jump1 = 1;
				read_1_jump = 1;
				goto read1_end;
			}

			/**
			if (strcmp("@SRR1532535.1001096", current_read1.name) == 0)
			{
				fprintf(stderr, "\n\n%s(1)\n", current_read1.name);
				fprintf(stderr, "number_of_hits1: %d\n", number_of_hits1);
				fprintf(stderr, "match_length1: %d\n", match_length1);
				fprintf(stderr, "current_read1.length: %d\n", current_read1.length);


			}
			**/





			///只有第一个seed才有可能出现这种情况
			if (match_length1 == current_read1.length
				&&
				number_of_hits1 > 1
				&&
				number_of_hits1 <= max_candidates_occ)
			{
				is_mutiple_map1 = 1;

				///等于-1说明read里面没有C
				if (C_site1 == -1)
				{
					locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
					reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

					///排序是为了最后好合并双端的结果
					std::sort(locates1, locates1 + tmp_SA_length);

					best_mapp_occ1 = tmp_SA_length;

					for (inner_i = 0; inner_i < tmp_SA_length; inner_i++)
					{
						candidates_votes1[inner_i].site = locates1[inner_i];
						candidates_votes1[inner_i].err = 0;
						candidates_votes1[inner_i].end_site = current_read1.length - 1;
					}




					///多匹配
					return_flag1 = 2;
					direct_jump1 = 1;
					read_1_jump = 1;
					goto read1_end;
				}
			}




			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
			if (number_of_hits1 == 1)
			{

				tmp_SA_length = 1;

				candidate_length1 = candidate_length1 + tmp_SA_length;

				locates1 = candidates1 + candidate_length1;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
			{

				if (number_of_hits1 != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length);
					reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

					candidate_length1 = candidate_length1 + tmp_SA_length;
					locates1 = candidates1 + candidate_length1;
				}

			}
			else
			{
				////新加的
				full_seed_id1--;
				////新加的
			}





			

			



			if (match_length1 == 0)
			{
				total_match_length1 = total_match_length1 + match_step;
			}
			else
			{
				total_match_length1 = total_match_length1 + match_length1 / 2;
			}


			seed_id1++;

		}
	read1_end:
		


		/*************************************第二个read的第一个种子********************************************/
		
		/*******控制种子数量************************/
		max_seed_number2 = current_read2.length / 10 - 1;
		if (max_seed_number2 > max_max_seed_number2)
		{
			max_seed_number2 = max_max_seed_number2;
		}
		/*******控制种子数量************************/

		unique_best_read2 = -1;

		direct_jump2 = 0;

		best_mapp_occ2 = 0;


		C_to_T_forward(current_read2.seq, bsSeq2, current_read2.length, &C_site2);

		total_match_length2 = 0;

		seed_id2 = 0;

		full_seed_id2 = 0;


		locates2 = candidates2;

		candidate_length2 = 0;

		is_mutiple_map2 = 0;

		///这里要改
		get_error2 = -1;
		///这里要改
		extra_seed_flag2 = 1;


		////5'端第一个seed单独处理
		if (seed_id2 < max_seed_number2 && total_match_length2 < current_read2.length)
		{
			current_seed_length2 = current_read2.length - total_match_length2;


			number_of_hits2 =
				count_backward_as_much_1_terminate
				(bsSeq2, current_seed_length2, &top, &bot, &pre_top, &pre_bot, &match_length2);


			///这里要改
			first_seed_match_length2 = match_length2;

			/**
			if (strcmp("@SRR1532535.50167", current_read2.name) == 0)
			{
				fprintf(stderr, "\n\n%s(2)\n", current_read2.name);
				fprintf(stderr, "number_of_hits2: %d\n", number_of_hits2);
				fprintf(stderr, "match_length2: %d\n", match_length2);
			}
			**/


			
			/**
			有两种情况能直接跳转
			1. 找到一个精确的唯一匹配
			2. 找到一堆精确的唯一匹配
			这两种情况best_mapp_occ2都不是0
			所以reseed不用担心这两种情况
			**/
			////新加的
			read2_seed_start[full_seed_id2] = total_match_length2;
			read2_seed_length[full_seed_id2] = match_length2;
			full_seed_id2++;
			////新加的


			if (number_of_hits2 == 1 && try_process_unique_mismatch_end_to_end_return_site(
				current_read2.length, current_read2.seq, current_read2.rseq, current_read2.qual, cigar,
				current_read2.name,
				locates2,
				top, tmp_ref, &match_length2, current_read2.length - C_site2 - 1,
				&candidates_votes2[0].site, &candidates_votes2[0].err, &get_error2))
			{
				candidates_votes2[0].end_site = current_read2.length - 1;
				best_mapp_occ2 = 1;
				///唯一匹配
				return_flag2 = 1;
				direct_jump2 = 1;
				read_2_jump = 1;
				goto read2_end;
			}

			/**
			if (strcmp("@SRR1532535.50167", current_read2.name) == 0)
			{
				///fprintf(stderr, "%s\n", current_read2.name);
				///fprintf(stderr, "number_of_hits2: %d\n", number_of_hits2);
				fprintf(stderr, "match_length2: %d\n", match_length2);
			}
			**/


			///只有第一个seed才有可能出现这种情况
			if (match_length2 == current_read2.length
				&&
				number_of_hits2 > 1
				&&
				number_of_hits2 <= max_candidates_occ)
			{
				is_mutiple_map2 = 1;

				///等于-1说明read里面没有C
				if (C_site2 == -1)
				{
					locate(bsSeq2 + current_seed_length2 - match_length2, top, bot, pre_top, pre_bot, locates2, match_length2, &tmp_SA_length);
					reverse_and_adjust_site(locates2, &tmp_SA_length, match_length2, total_match_length2);

					///排序是为了最后好合并双端的结果
					std::sort(locates2, locates2 + tmp_SA_length);

					best_mapp_occ2 = tmp_SA_length;

					for (inner_i = 0; inner_i < tmp_SA_length; inner_i++)
					{
						candidates_votes2[inner_i].site = locates2[inner_i];
						candidates_votes2[inner_i].err = 0;
						candidates_votes2[inner_i].end_site = current_read2.length - 1;
					}




					///多匹配
					return_flag2 = 2;
					direct_jump2 = 1;
					read_2_jump = 1;
					goto read2_end;
				}
			}




			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
			if (number_of_hits2 == 1)
			{

				tmp_SA_length = 1;

				candidate_length2 = candidate_length2 + tmp_SA_length;

				locates2 = candidates2 + candidate_length2;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length2 >= available_seed_length && number_of_hits2 <= max_seed_matches)
			{

				if (number_of_hits2 != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq2 + current_seed_length2 - match_length2, top, bot, pre_top, pre_bot, locates2, match_length2, &tmp_SA_length);
					reverse_and_adjust_site(locates2, &tmp_SA_length, match_length2, total_match_length2);

					candidate_length2 = candidate_length2 + tmp_SA_length;
					locates2 = candidates2 + candidate_length2;
				}

			}
			else
			{
				////新加的
				full_seed_id2--;
				////新加的
			}







			if (match_length2 == 0)
			{
				total_match_length2 = total_match_length2 + match_step;
			}
			else
			{
				total_match_length2 = total_match_length2 + match_length2 / 2;
			}


			seed_id2++;

		}

		read2_end:

		/*************************************处理第一个种子********************************************/
		
		/**
		if (strcmp("@SRR1532535.50167", current_read2.name) == 0)
		{
			///fprintf(stderr, "%s\n", current_read2.name);
			///fprintf(stderr, "number_of_hits2: %d\n", number_of_hits2);
			fprintf(stderr, "candidate_length1: %d\n", candidate_length1);
			fprintf(stderr, "candidate_length2: %d\n", candidate_length2);
		}
		**/

		
		///if (number_of_hits1 <= number_of_hits2 && candidate_length1 <= candidate_length2)
		if (candidate_length1 <= candidate_length2)
		{
			/*************************************第一个read********************************************/
			if (read_1_jump == 0)
			{

				best_mapp_occ1 = process_rest_seed_debug(
					seed_id1,
					max_seed_number1,
					total_match_length1,
					&current_read1,
					current_seed_length1,
					number_of_hits1,
					bsSeq1,
					match_length1,
					locates1,
					candidates1,
					candidate_length1,
					available_seed_length,
					max_seed_matches,
					match_step,
					candidates_votes1,
					candidates_votes_length1,
					error_threshold1,
					tmp_ref,
					Peq_SSE,
					min_err,
					get_error1,
					first_seed_match_length1,
					read1_seed_start,
					read1_seed_length,
					&full_seed_id1);


			}


			if (best_mapp_occ1 == 0)
			{
				i++;
				continue;
			}



			//*************************************第一个read********************************************



			


			/*************************************第二个read********************************************/

			if (read_2_jump == 0)
			{

				best_mapp_occ2 = process_rest_seed_filter_debug(
					seed_id2,
					max_seed_number2,
					total_match_length2,
					&current_read2,
					current_seed_length2,
					number_of_hits2,
					bsSeq2,
					match_length2,
					locates2,
					candidates2,
					candidate_length2,
					available_seed_length,
					max_seed_matches,
					match_step,
					candidates_votes2,
					candidates_votes_length2,
					error_threshold2,
					tmp_ref,
					Peq_SSE,
					min_err,
					candidates_votes1,
					best_mapp_occ1,
					inner_maxDistance_pair,
					inner_minDistance_pair,
					get_error2,
					first_seed_match_length2,
					read2_seed_start,
					read2_seed_length,
					&full_seed_id2);




			}








			/*************************************第二个read********************************************/
		}
		else
		{


			/*************************************第二个read********************************************/

			if (read_2_jump == 0)
			{
				
				best_mapp_occ2 = process_rest_seed_debug(
				seed_id2,
				max_seed_number2,
				total_match_length2,
				&current_read2,
				current_seed_length2,
				number_of_hits2,
				bsSeq2,
				match_length2,
				locates2,
				candidates2,
				candidate_length2,
				available_seed_length,
				max_seed_matches,
				match_step,
				candidates_votes2,
				candidates_votes_length2,
				error_threshold2,
				tmp_ref,
				Peq_SSE,
				min_err,
				get_error2,
				first_seed_match_length2,
				read2_seed_start,
				read2_seed_length,
				&full_seed_id2);


				
			}


			if (best_mapp_occ2 == 0)
			{
				i++;
				continue;
			}


			/*************************************第二个read********************************************/






			





			/*************************************第一个read********************************************/
			if (read_1_jump == 0)
			{

				best_mapp_occ1 = process_rest_seed_filter_debug(
					seed_id1,
					max_seed_number1,
					total_match_length1,
					&current_read1,
					current_seed_length1,
					number_of_hits1,
					bsSeq1,
					match_length1,
					locates1,
					candidates1,
					candidate_length1,
					available_seed_length,
					max_seed_matches,
					match_step,
					candidates_votes1,
					candidates_votes_length1,
					error_threshold1,
					tmp_ref,
					Peq_SSE,
					min_err,
					candidates_votes2,
					best_mapp_occ2,
					inner_maxDistance_pair,
					inner_minDistance_pair,
					get_error1,
					first_seed_match_length1,
					read1_seed_start,
					read1_seed_length,
					&full_seed_id1);


			}




			//*************************************第一个read********************************************

			
			
			


		}


	

		///if (best_mapp_occ1 != 0 && best_mapp_occ2 == 0)
		if (best_mapp_occ2 == 0)
		{
			////fprintf(stderr, "\nread2: %s\n", current_read2.name);

			best_mapp_occ2 = reseed_filter(
				seed_id2,
				max_seed_number2,
				total_match_length2,
				&current_read2,
				current_seed_length2,
				number_of_hits2,
				bsSeq2,
				match_length2,
				locates2,
				candidates2,
				candidate_length2,
				available_seed_length,
				max_seed_matches,
				match_step,
				candidates_votes2,
				candidates_votes_length2,
				error_threshold2,
				tmp_ref,
				Peq_SSE,
				min_err,
				candidates_votes1,
				best_mapp_occ1,
				inner_maxDistance_pair,
				inner_minDistance_pair,
				get_error2,
				first_seed_match_length2,
				read2_seed_start,
				read2_seed_length,
				full_seed_id2,
				result_seed_start,
				result_seed_length);
		}
		///else if (best_mapp_occ1 == 0 && best_mapp_occ2 != 0)
		else if (best_mapp_occ1 == 0)
		{

			///fprintf(stderr, "\nread1: %s\n", current_read1.name);

			best_mapp_occ1 = reseed_filter(
				seed_id1,
				max_seed_number1,
				total_match_length1,
				&current_read1,
				current_seed_length1,
				number_of_hits1,
				bsSeq1,
				match_length1,
				locates1,
				candidates1,
				candidate_length1,
				available_seed_length,
				max_seed_matches,
				match_step,
				candidates_votes1,
				candidates_votes_length1,
				error_threshold1,
				tmp_ref,
				Peq_SSE,
				min_err,
				candidates_votes2,
				best_mapp_occ2,
				inner_maxDistance_pair,
				inner_minDistance_pair,
				get_error1,
				first_seed_match_length1,
				read1_seed_start,
				read1_seed_length,
				full_seed_id1,
				result_seed_start,
				result_seed_length);
		}


		
	
		
		
		mapping_pair = 0;

		if (unique_best_read1 != -1 && unique_best_read2 != -1)
		{

			if (candidates_votes1[unique_best_read1].site > candidates_votes2[unique_best_read2].site)
			{
				distance = candidates_votes1[unique_best_read1].site - candidates_votes2[unique_best_read2].site;
			}
			else
			{
				distance = candidates_votes2[unique_best_read2].site - candidates_votes1[unique_best_read1].site;
			}

			if (distance <= inner_maxDistance_pair
				&&
				distance >= inner_minDistance_pair)
			{
				best_pair_1_index = unique_best_read1;
				best_pair_2_index = unique_best_read2;
				mapping_pair = 1;
			}
		}




		if (mapping_pair == 0)
		{
			mapping_pair = new_faster_verify_pairs(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
				candidates_votes1, candidates_votes2,
				&best_pair_1_index, &best_pair_2_index,
				inner_maxDistance_pair, inner_minDistance_pair);
		}
		

	end_i:
		if (mapping_pair == 1)
		{
			

			//****************第一个read的后处理
			result1.err = candidates_votes1[best_pair_1_index].err;
			result1.origin_site = candidates_votes1[best_pair_1_index].site;
			result1.end_site = candidates_votes1[best_pair_1_index].end_site;

			if (result1.err!=0)
			{
				calculate_best_map_cigar_end_to_end_return
					(&min_candidates_votes_length,
					current_read1.length, error_threshold1, tmp_ref,
					current_read1.seq, current_read1.rseq, current_read1.qual, cigar, path, matrix, matrix_bit,
					current_read1.name,
					result1.cigar,
					&result1, &matched_length1);

				

			}
			else
			{

				output_sam_end_to_end_return(
					result1.origin_site,
					result1.end_site,
					result1.end_site + 1 - current_read1.length,
					current_read1.length,
					&result1
					);

				sprintf(result1.cigar, "%dM", current_read1.length);
				matched_length1 = current_read1.length;
			}
			//****************第一个read的后处理




			//****************第二个read的后处理
			result2.err = candidates_votes2[best_pair_2_index].err;
			result2.origin_site = candidates_votes2[best_pair_2_index].site;
			result2.end_site = candidates_votes2[best_pair_2_index].end_site;

			if (result2.err != 0)
			{
				calculate_best_map_cigar_end_to_end_return
					(&min_candidates_votes_length,
					current_read2.length, error_threshold2, tmp_ref,
					current_read2.seq, current_read2.rseq, current_read2.qual, cigar, path, matrix, matrix_bit,
					current_read2.name,
					result2.cigar,
					&result2, &matched_length2);
			}
			else
			{

				output_sam_end_to_end_return(
					result2.origin_site,
					result2.end_site,
					result2.end_site + 1 - current_read2.length,
					current_read2.length,
					&result2
					);

				sprintf(result2.cigar, "%dM", current_read2.length);
				matched_length2 = current_read2.length;
			}
			//****************第二个read的后处理

			

			if (result2.site >= result1.site)
			{
				paired_end_distance = result2.site - result1.site + matched_length2;

			}
			else
			{
				paired_end_distance = result1.site - result2.site + matched_length1;
			}

			if (paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
			{
				unique_matched_read++;


				directly_output_read1(current_read1.name, current_read1.seq, current_read1.rseq, current_read1.qual,
					&result1, &result2, current_read1.length,
					matched_length1, matched_length2);

				directly_output_read2(current_read2.name, current_read2.seq, current_read2.rseq, current_read2.qual,
					&result2, &result1, current_read2.length,
					matched_length2, matched_length1);
			}

			
			


		}
		else if (mapping_pair > 1)
		{
			ambious_matched_read++;
		}

		/**
		if (mapping_pair < 1 && (best_mapp_occ1 == 0 || best_mapp_occ2 == 0))
		{
			fprintf(stderr, "%s\n", current_read1.name);
			fprintf(stderr, "best_mapp_occ1: %llu\n", best_mapp_occ1);
			fprintf(stderr, "best_mapp_occ2: %llu\n", best_mapp_occ2);
			fprintf(stderr, "*************\n");
		}
		**/



	
		/**
		fprintf(stderr, "%s\n", current_read1.name);
		fprintf(stderr, "%llu\n", best_mapp_occ1);
		fprintf(stderr, "%llu\n", best_mapp_occ2);
		fprintf(stderr, "*************\n");
		
		if (best_mapp_occ1&&best_mapp_occ2)
		{
			total_1++;
		}
		**/
		


		i++;
	}


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;


	fprintf(stderr, "debug_1: %lld\n", debug_1);
	fprintf(stderr, "debug_2: %lld\n", debug_2);
	

	return 1;
}






void* Map_Pair_Seq_split_fast(void* arg)
{

	int thread_id = *((int *)arg);

	FILE* output_file = get_Ouput_Dec();

	bwt_locate_queue get_queue;

	init_locate_queue_muti_thread(&get_queue);



	long long enq_i = 0;
	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;


	char** ref = (char**)malloc(sizeof(char*)* 8);

	for (i = 0; i<8; i++)
	{
		ref[i] = (char*)malloc(sizeof(char)*(SEQ_LENGTH + 2 * thread_e));
	}


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;

	///Read current_read1, current_read2;

	//bitmapper_bs_iter number_of_hits = 0;
	bitmapper_bs_iter number_of_hits1, number_of_hits2;

	unsigned int bsSeq_max_length = 1000;
	///char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	///bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter top, bot, match_step;
	bitmapper_bs_iter total_match_length1, total_match_length2;
	bitmapper_bs_iter match_length1, match_length2;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;
	///bitmapper_bs_iter available_seed_length = 30;

	///min_seed_length这个值不能改，改了就会报错
	///bitmapper_bs_iter min_seed_length = 17;
	///bitmapper_bs_iter min_seed_max_seed_matches = 100;

	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	bitmapper_bs_iter short_seed_max_seed_matches = 20;


	available_seed_length = over_all_seed_length;

	/**
	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	**/

	//bitmapper_bs_iter current_seed_length;
	bitmapper_bs_iter current_seed_length1, current_seed_length2;

	///bitmapper_bs_iter* locates;
	bitmapper_bs_iter *locates1, *locates2;
	///bitmapper_bs_iter* candidates;
	bitmapper_bs_iter *candidates1, *candidates2;
	///bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter candidate_length1, candidate_length2;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	///seed_votes* candidates_votes;
	seed_votes *candidates_votes1, *candidates_votes2;
	///bitmapper_bs_iter candidates_votes_length;
	bitmapper_bs_iter candidates_votes_length1, candidates_votes_length2;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	//bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter max_seed_number1, max_seed_number2;
	///bitmapper_bs_iter seed_id;
	bitmapper_bs_iter seed_id1, seed_id2;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter ambious_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	///bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;
	///int is_duplicate;
	///bitmapper_bs_iter min_index;

	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
		error_threshold = thread_e;
	}
	**/

	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;


	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	///max_seed_number = 9;
	max_seed_number1 = 25;
	max_seed_number2 = 25;
	int max_max_seed_number1 = max_seed_number1;
	int max_max_seed_number2 = max_seed_number2;



	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(max_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	double re_seed_threshold = 0.80;
	int reseed_length;
	bitmapper_bs_iter re_seed_error_threshold;

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	long long total_number_of_hit = 0;


	int is_mutiple_map1, is_mutiple_map2;

	long long distance;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	long long direct_cut = 0;



	int add_empty_read[4] = { 1, 0, 0, 0 };
	int add_matched_read[4] = { 0, 1, 1, 0 };
	int add_unique_matched_read[4] = { 0, 1, 0, 0 };
	int add_unmatched_read[4] = { 0, 0, 0, 1 };
	int best_mapp_occ1 = 0;
	int best_mapp_occ2 = 0;
	int inner_i, inner_j;
	long long mapping_pair;
	long long best_pair_1_index;
	long long best_pair_2_index;
	///long long best_sum_err = 2 * error_threshold + 1;
	long long best_sum_err;
	long long current_sum_err;

	long long unique_best_read1;
	long long unique_best_read2;
	map_result result1, result2;


	int max_candidates_occ = max_seed_number1 * max_seed_matches;


	candidates1 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates2 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates_votes1 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));
	candidates_votes2 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));

	seed_votes *tmp_votes1, *tmp_votes2;
	tmp_votes1 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));
	tmp_votes2 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));


	///这是个特殊情况，仅在很偶尔的地方才有用
	///就是一边read有超过max_candidates_occ个误配的时候，就不要了
	if (max_candidates_occ > 10000)
	{
		max_candidates_occ = 10000;
	}


	int matched_length1, matched_length2;



	char* bsSeq1 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	char* bsSeq2 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	int C_site1;
	int C_site2;
	int direct_jump1 = 0;
	int direct_jump2 = 0;
	long long current_best_err;

	int first_seed_read1_occ, first_seed_read2_occ;
	int first_seed_read1_length, first_seed_read2_length;



	Read* read_batch1;
	Read* read_batch2;

	Read_buffer_pe_sub_block curr_sub_block;

	init_single_sub_block_pe(&curr_sub_block);

	Output_buffer_sub_block current_sub_buffer;
	init_buffer_sub_block(&current_sub_buffer);

	double startTime = Get_T();

	int file_flag;
	int return_flag1, return_flag2;
	int read_1_jump, read_2_jump;

	int inner_maxDistance_pair;
	int inner_minDistance_pair;

	int obtain_reads_num;
	file_flag = 1;


	int get_error1, get_error2;
	long long first_seed_match_length1, first_seed_match_length2;
	long long second_seed_length1, second_seed_length2;
	int extra_seed_flag1, extra_seed_flag2;
	int paired_end_distance;

	//正向模式
	i = 0;
	while (file_flag != 0)
	{


		file_flag = get_pe_reads_mul_thread(&curr_sub_block);
		obtain_reads_num = curr_sub_block.sub_block_read_number;
		read_batch1 = curr_sub_block.read1;
		read_batch2 = curr_sub_block.read2;




		enq_i = enq_i + obtain_reads_num;

		i = 0;


		while (i < obtain_reads_num)
		{


			error_threshold1 = thread_e_f * read_batch1[i].length;
			if (error_threshold1 >= max_error_threshold)
			{
				error_threshold1 = max_error_threshold;
			}


			error_threshold2 = thread_e_f * read_batch2[i].length;
			if (error_threshold2 >= max_error_threshold)
			{
				error_threshold2 = max_error_threshold;
			}

			if (error_threshold1 > error_threshold2)
			{
				large_error_threshold = error_threshold1;
			}
			else
			{
				large_error_threshold = error_threshold2;
			}

			best_sum_err = 2 * large_error_threshold + 1;


			inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
			inner_minDistance_pair = minDistance_pair - large_error_threshold * 2;


			read_1_jump = 0;
			read_2_jump = 0;



			//*******控制种子数量************************
			max_seed_number1 = read_batch1[i].length / 10 - 1;
			if (max_seed_number1 > max_max_seed_number1)
			{
				max_seed_number1 = max_max_seed_number1;
			}
			//*******控制种子数量************************




			/*******控制种子数量************************/
			max_seed_number2 = read_batch2[i].length / 10 - 1;
			if (max_seed_number2 > max_max_seed_number2)
			{
				max_seed_number2 = max_max_seed_number2;
			}
			/*******控制种子数量************************/

			

			get_candidates_muti_thread(
				max_seed_number1,
				max_candidates_occ,
				available_seed_length,
				max_seed_matches,
				match_step,
				error_threshold1,
				candidates_votes_length1,
				&read_batch1[i],
				bsSeq1,
				cigar,
				tmp_ref,
				candidates1,
				candidates_votes1,
				&best_mapp_occ1,
				&candidates_votes_length1,
				&get_queue);


			get_candidates_muti_thread(
				max_seed_number2,
				max_candidates_occ,
				available_seed_length,
				max_seed_matches,
				match_step,
				error_threshold2,
				candidates_votes_length2,
				&read_batch2[i],
				bsSeq2,
				cigar,
				tmp_ref,
				candidates2,
				candidates_votes2,
				&best_mapp_occ2,
				&candidates_votes_length2,
				&get_queue);


			if (best_mapp_occ1 > 0 && best_mapp_occ2 > 0)   ///说明两个都不要验证
			{
				best_mapp_occ1 = candidates_votes_length1;
				best_mapp_occ2 = candidates_votes_length2;
			}
			else
			{

				if (best_mapp_occ1 == 0 || best_mapp_occ2 == 0)
				{
					i++;
					continue;
				}
			

				filter_pairs(
					candidates_votes_length1,
					candidates_votes_length2,
					&candidates_votes1,
					&candidates_votes2,
					&tmp_votes1,
					&tmp_votes2,
					&candidates_votes_length1,
					&candidates_votes_length2,
					inner_maxDistance_pair,
					inner_minDistance_pair
					);

				if (candidates_votes_length1 == 0 || candidates_votes_length2 == 0)
				{
					i++;
					continue;
				}

				/**
				当程序执行到这里时，有三种情况
				1.
				best_mapp_occ = -1，这说明还有候选位置没有验证，
				此时candidates_votes都是候选位置，candidates_votes_length是候选位置数目
				2.
				best_mapp_occ>0，此时说明已经找到了匹配位置
				3.
				best_mapp_occ = 0，此时说明不匹配

				candidates_votes中的位置都是排过序的
				**/
				if (best_mapp_occ1 == -1 && best_mapp_occ2 == -1)///这是两个read都需要验证的情况
				{

				

					if (candidates_votes_length1 <= candidates_votes_length2)
					{
						best_mapp_occ1 = verify_candidate_locations(
							error_threshold1,
							candidates_votes_length1,
							&read_batch1[i],
							bsSeq1,
							cigar,
							tmp_ref,
							candidates_votes1,
							Peq_SSE);


						if (best_mapp_occ1 == 0)
						{
							i++;
							continue;
						}

						///用1的结果过滤2
						///所以1用的是best_mapp_occ1, 而2用的是candidates_votes_length2
						filter_pairs_single_side(
							best_mapp_occ1,
							candidates_votes_length2,
							&candidates_votes1,
							&candidates_votes2,
							&candidates_votes_length2,
							inner_maxDistance_pair,
							inner_minDistance_pair
							);

						best_mapp_occ2 = verify_candidate_locations(
							error_threshold2,
							candidates_votes_length2,
							&read_batch2[i],
							bsSeq2,
							cigar,
							tmp_ref,
							candidates_votes2,
							Peq_SSE);
					}
					else
					{

						best_mapp_occ2 = verify_candidate_locations(
							error_threshold2,
							candidates_votes_length2,
							&read_batch2[i],
							bsSeq2,
							cigar,
							tmp_ref,
							candidates_votes2,
							Peq_SSE);


						if (best_mapp_occ2 == 0)
						{
							i++;
							continue;
						}



						///用2的结果过滤1
						///所以2用的是best_mapp_occ2, 而1用的是candidates_votes_length1
						filter_pairs_single_side(
							best_mapp_occ2,
							candidates_votes_length1,
							&candidates_votes2,
							&candidates_votes1,
							&candidates_votes_length1,
							inner_maxDistance_pair,
							inner_minDistance_pair
							);


						best_mapp_occ1 = verify_candidate_locations(
							error_threshold1,
							candidates_votes_length1,
							&read_batch1[i],
							bsSeq1,
							cigar,
							tmp_ref,
							candidates_votes1,
							Peq_SSE);
					}

				}
				else if (best_mapp_occ1 != -1)
				{

					if (candidates_votes_length1 < best_mapp_occ1)
					{
						best_mapp_occ1 = candidates_votes_length1;
					}

					if (best_mapp_occ2 == -1)
					{
						best_mapp_occ2 = verify_candidate_locations(
							error_threshold2,
							candidates_votes_length2,
							&read_batch2[i],
							bsSeq2,
							cigar,
							tmp_ref,
							candidates_votes2,
							Peq_SSE);
					}
				}
				else if (best_mapp_occ2 != -1)
				{
					if (candidates_votes_length2 < best_mapp_occ2)
					{
						best_mapp_occ2 = candidates_votes_length2;
					}



					if (best_mapp_occ1 == -1)
					{
						best_mapp_occ1 = verify_candidate_locations(
							error_threshold1,
							candidates_votes_length1,
							&read_batch1[i],
							bsSeq1,
							cigar,
							tmp_ref,
							candidates_votes1,
							Peq_SSE);
					}
				}

			}

			mapping_pair = 0;

			mapping_pair = new_faster_verify_pairs(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
					candidates_votes1, candidates_votes2,
					&best_pair_1_index, &best_pair_2_index,
					inner_maxDistance_pair, inner_minDistance_pair);



		end_i:
			if (mapping_pair == 1)
			{
				///unique_matched_read++;

				//****************第一个read的后处理
				result1.err = candidates_votes1[best_pair_1_index].err;
				result1.origin_site = candidates_votes1[best_pair_1_index].site;
				result1.end_site = candidates_votes1[best_pair_1_index].end_site;

				if (result1.err != 0)
				{
					calculate_best_map_cigar_end_to_end_return
						(&min_candidates_votes_length,
						read_batch1[i].length, error_threshold1, tmp_ref,
						read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual, cigar, path, matrix, matrix_bit,
						read_batch1[i].name,
						result1.cigar,
						&result1, &matched_length1);



				}
				else
				{

					output_sam_end_to_end_return(
						result1.origin_site,
						result1.end_site,
						result1.end_site + 1 - read_batch1[i].length,
						read_batch1[i].length,
						&result1
						);

					sprintf(result1.cigar, "%dM", read_batch1[i].length);
					matched_length1 = read_batch1[i].length;
				}
				//****************第一个read的后处理




				//****************第二个read的后处理
				result2.err = candidates_votes2[best_pair_2_index].err;
				result2.origin_site = candidates_votes2[best_pair_2_index].site;
				result2.end_site = candidates_votes2[best_pair_2_index].end_site;

				if (result2.err != 0)
				{
					calculate_best_map_cigar_end_to_end_return
						(&min_candidates_votes_length,
						read_batch2[i].length, error_threshold2, tmp_ref,
						read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual, cigar, path, matrix, matrix_bit,
						read_batch2[i].name,
						result2.cigar,
						&result2, &matched_length2);
				}
				else
				{

					output_sam_end_to_end_return(
						result2.origin_site,
						result2.end_site,
						result2.end_site + 1 - read_batch2[i].length,
						read_batch2[i].length,
						&result2
						);

					sprintf(result2.cigar, "%dM", read_batch2[i].length);
					matched_length2 = read_batch2[i].length;
				}
				//****************第二个read的后处理



				if (result2.site >= result1.site)
				{
					paired_end_distance = result2.site - result1.site + matched_length2;

				}
				else
				{
					paired_end_distance = result1.site - result2.site + matched_length1;
				}

				if (paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
				{
					unique_matched_read++;
					/**
					directly_output_read1_return_result(read_batch1[i].name, read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual,
						&result1, &result2, read_batch1[i].length,
						matched_length1, matched_length2, result_string);

					directly_output_read2_return_result(read_batch2[i].name, read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual,
						&result2, &result1, read_batch2[i].length,
						matched_length2, matched_length1, result_string);
					**/

					directly_output_read1_return_buffer(read_batch1[i].name, read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual,
						&result1, &result2, read_batch1[i].length,
						matched_length1, matched_length2, &current_sub_buffer);

					directly_output_read2_return_buffer(read_batch2[i].name, read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual,
						&result2, &result1, read_batch2[i].length,
						matched_length2, matched_length1, &current_sub_buffer);

				}

				

			}
			else if (mapping_pair > 1)
			{
				ambious_matched_read++;

				
			}


			i++;
		}

		current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

		push_results_to_buffer(&current_sub_buffer);

		current_sub_buffer.length = 0;

	}

	finish_output_buffer();


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;

	fprintf(stdout, "thread %d completed!\n", thread_id);


}



///输出成对匹配
///和一端unique一端不匹配的结果
void* Map_Pair_Seq_split(void* arg)
{

	int thread_id = *((int *)arg);

	FILE* output_file = get_Ouput_Dec();

	bwt_locate_queue get_queue;

	init_locate_queue_muti_thread(&get_queue);



	long long enq_i = 0;
	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;


	char** ref = (char**)malloc(sizeof(char*)* 8);

	for (i = 0; i<8; i++)
	{
		ref[i] = (char*)malloc(sizeof(char)*(SEQ_LENGTH + 2 * thread_e));
	}


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;

	///Read current_read1, current_read2;

	//bitmapper_bs_iter number_of_hits = 0;
	bitmapper_bs_iter number_of_hits1, number_of_hits2;

	unsigned int bsSeq_max_length = 1000;
	///char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	///bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter top, bot, match_step;
	bitmapper_bs_iter total_match_length1, total_match_length2;
	bitmapper_bs_iter match_length1, match_length2;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;
	///bitmapper_bs_iter available_seed_length = 30;

	///min_seed_length这个值不能改，改了就会报错
	///bitmapper_bs_iter min_seed_length = 17;
	///bitmapper_bs_iter min_seed_max_seed_matches = 100;

	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	bitmapper_bs_iter short_seed_max_seed_matches = 20;

	available_seed_length = over_all_seed_length;

	/**
	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	**/

	//bitmapper_bs_iter current_seed_length;
	bitmapper_bs_iter current_seed_length1, current_seed_length2;

	///bitmapper_bs_iter* locates;
	bitmapper_bs_iter *locates1, *locates2;
	///bitmapper_bs_iter* candidates;
	bitmapper_bs_iter *candidates1, *candidates2;
	///bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter candidate_length1, candidate_length2;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	///seed_votes* candidates_votes;
	seed_votes *candidates_votes1, *candidates_votes2;
	///bitmapper_bs_iter candidates_votes_length;
	bitmapper_bs_iter candidates_votes_length1, candidates_votes_length2;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	//bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter max_seed_number1, max_seed_number2;
	///bitmapper_bs_iter seed_id;
	bitmapper_bs_iter seed_id1, seed_id2;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter ambious_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	///bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;
	///int is_duplicate;
	///bitmapper_bs_iter min_index;

	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;
	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
		error_threshold = thread_e;
	}
	**/


	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	///max_seed_number = 9;
	max_seed_number1 = 25;
	max_seed_number2 = 25;
	int max_max_seed_number1 = max_seed_number1;
	int max_max_seed_number2 = max_seed_number2;



	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(max_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	double re_seed_threshold = 0.80;
	int reseed_length;
	bitmapper_bs_iter re_seed_error_threshold;

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	long long total_number_of_hit = 0;


	int is_mutiple_map1, is_mutiple_map2;

	long long distance;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	long long direct_cut = 0;



	int add_empty_read[4] = { 1, 0, 0, 0 };
	int add_matched_read[4] = { 0, 1, 1, 0 };
	int add_unique_matched_read[4] = { 0, 1, 0, 0 };
	int add_unmatched_read[4] = { 0, 0, 0, 1 };
	int best_mapp_occ1 = 0;
	int best_mapp_occ2 = 0;
	int inner_i, inner_j;
	long long mapping_pair;
	long long best_pair_1_index;
	long long best_pair_2_index;
	long long best_sum_err;
	long long current_sum_err;

	long long unique_best_read1;
	long long unique_best_read2;
	map_result result1, result2;


	int max_candidates_occ = max_seed_number1 * max_seed_matches *  2.5;


	candidates1 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates2 = (bitmapper_bs_iter *)malloc(max_candidates_occ * sizeof(bitmapper_bs_iter));
	candidates_votes1 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));
	candidates_votes2 = (seed_votes *)malloc(max_candidates_occ * sizeof(seed_votes));

	///这是个特殊情况，仅在很偶尔的地方才有用
	///就是一边read有超过max_candidates_occ个误配的时候，就不要了
	if (max_candidates_occ > 10000)
	{
		max_candidates_occ = 10000;
	}


	int matched_length1, matched_length2;



	char* bsSeq1 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	char* bsSeq2 = (char*)malloc(sizeof(char)* bsSeq_max_length);
	int C_site1;
	int C_site2;
	int direct_jump1 = 0;
	int direct_jump2 = 0;
	long long current_best_err;

	int first_seed_read1_occ, first_seed_read2_occ;
	int first_seed_read1_length, first_seed_read2_length;






	///fprintf(stdout, "max_candidates_occ: %lld\n", max_candidates_occ);



	Read* read_batch1;
	Read* read_batch2;

	Read_buffer_pe_sub_block curr_sub_block;

	init_single_sub_block_pe(&curr_sub_block);

	Output_buffer_sub_block current_sub_buffer;
	init_buffer_sub_block(&current_sub_buffer);

	


	int* read1_seed_start = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* read1_seed_length = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* read2_seed_start = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* read2_seed_length = (int*)malloc(sizeof(int)* max_max_seed_number1);

	int* result_seed_start = (int*)malloc(sizeof(int)* max_max_seed_number1);
	int* result_seed_length = (int*)malloc(sizeof(int)* max_max_seed_number1);




	double startTime = Get_T();

	int file_flag;
	int return_flag1, return_flag2;
	int read_1_jump, read_2_jump;

	int inner_maxDistance_pair;
	int inner_minDistance_pair;

	int obtain_reads_num;
	file_flag = 1;
	bitmapper_bs_iter full_seed_id1, full_seed_id2;


	int get_error1, get_error2;
	long long first_seed_match_length1, first_seed_match_length2;
	long long second_seed_length1, second_seed_length2;
	int extra_seed_flag1, extra_seed_flag2;
	int paired_end_distance;

	//正向模式
	i = 0;
	while (file_flag != 0)
	{


		file_flag = get_pe_reads_mul_thread(&curr_sub_block);
		obtain_reads_num = curr_sub_block.sub_block_read_number;
		read_batch1 = curr_sub_block.read1;
		read_batch2 = curr_sub_block.read2;




		enq_i = enq_i + obtain_reads_num;


		i = 0;


		while (i < obtain_reads_num)
		{
			error_threshold1 = thread_e_f * read_batch1[i].length;
			if (error_threshold1 >= max_error_threshold)
			{
				error_threshold1 = max_error_threshold;
			}


			error_threshold2 = thread_e_f * read_batch2[i].length;
			if (error_threshold2 >= max_error_threshold)
			{
				error_threshold2 = max_error_threshold;
			}

			if (error_threshold1 > error_threshold2)
			{
				large_error_threshold = error_threshold1;
			}
			else
			{
				large_error_threshold = error_threshold2;
			}

			best_sum_err = 2 * large_error_threshold + 1;

			



			inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
			inner_minDistance_pair = minDistance_pair - large_error_threshold * 2;




			/*************************************处理第一个种子********************************************/

			read_1_jump = 0;
			read_2_jump = 0;


			//*************************************第一个read的第一个种子********************************************





			//*******控制种子数量************************
			max_seed_number1 = read_batch1[i].length / 10 - 1;
			if (max_seed_number1 > max_max_seed_number1)
			{
				max_seed_number1 = max_max_seed_number1;
			}
			//*******控制种子数量************************



			





			unique_best_read1 = -1;

			direct_jump1 = 0;

			best_mapp_occ1 = 0;

			full_seed_id1 = 0;


			C_to_T_forward(read_batch1[i].seq, bsSeq1, read_batch1[i].length, &C_site1);

			total_match_length1 = 0;

			seed_id1 = 0;


			locates1 = candidates1;

			candidate_length1 = 0;

			is_mutiple_map1 = 0;


			///这里要改
			get_error1 = -1;
			///这里要改
			extra_seed_flag1 = 1;


			////5'端第一个seed单独处理
			if (seed_id1 < max_seed_number1 && total_match_length1 < read_batch1[i].length)
			{
				current_seed_length1 = read_batch1[i].length - total_match_length1;


				number_of_hits1 =
					count_backward_as_much_1_terminate
					(bsSeq1, current_seed_length1, &top, &bot, &pre_top, &pre_bot, &match_length1);

				///这里要改
				first_seed_match_length1 = match_length1;


				/**
				有两种情况能直接跳转
				1. 找到一个精确的唯一匹配
				2. 找到一堆精确的唯一匹配
				这两种情况best_mapp_occ1都不是0
				所以reseed不用担心这两种情况
				**/
				////新加的
				read1_seed_start[full_seed_id1] = total_match_length1;
				read1_seed_length[full_seed_id1] = match_length1;
				full_seed_id1++;
				////新加的




				///这里要改
				if (number_of_hits1 == 1 && try_process_unique_mismatch_end_to_end_return_site(
					read_batch1[i].length, read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual, cigar,
					read_batch1[i].name,
					locates1,
					top, tmp_ref, &match_length1, read_batch1[i].length - C_site1 - 1,
					&candidates_votes1[0].site, &candidates_votes1[0].err, &get_error1))
				{
					candidates_votes1[0].end_site = read_batch1[i].length - 1;
					best_mapp_occ1 = 1;
					///唯一匹配
					return_flag1 = 1;
					direct_jump1 = 1;
					read_1_jump = 1;
					goto read1_end;
				}






				///只有第一个seed才有可能出现这种情况
				if (match_length1 == read_batch1[i].length
					&&
					number_of_hits1 > 1
					&&
					number_of_hits1 <= max_candidates_occ)
				{
					is_mutiple_map1 = 1;

					///等于-1说明read里面没有C
					if (C_site1 == -1)
					{
						locate_muti_thread(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length, &get_queue);
						reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

						///排序是为了最后好合并双端的结果
						std::sort(locates1, locates1 + tmp_SA_length);

						best_mapp_occ1 = tmp_SA_length;

						for (inner_i = 0; inner_i < tmp_SA_length; inner_i++)
						{
							candidates_votes1[inner_i].site = locates1[inner_i];
							candidates_votes1[inner_i].err = 0;
							candidates_votes1[inner_i].end_site = read_batch1[i].length - 1;
						}




						///多匹配
						return_flag1 = 2;
						direct_jump1 = 1;
						read_1_jump = 1;
						goto read1_end;
					}
				}




				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
				if (number_of_hits1 == 1)
				{

					tmp_SA_length = 1;

					candidate_length1 = candidate_length1 + tmp_SA_length;

					locates1 = candidates1 + candidate_length1;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length1 >= available_seed_length && number_of_hits1 <= max_seed_matches)
				{

					if (number_of_hits1 != 0)
					{

						tmp_SA_length = 0;

						locate_muti_thread(bsSeq1 + current_seed_length1 - match_length1, top, bot, pre_top, pre_bot, locates1, match_length1, &tmp_SA_length, &get_queue);
						reverse_and_adjust_site(locates1, &tmp_SA_length, match_length1, total_match_length1);

						candidate_length1 = candidate_length1 + tmp_SA_length;
						locates1 = candidates1 + candidate_length1;
					}

				}
				else
				{
					////新加的
					full_seed_id1--;
					////新加的
				}











				if (match_length1 == 0)
				{
					total_match_length1 = total_match_length1 + match_step;
				}
				else
				{
					total_match_length1 = total_match_length1 + match_length1 / 2;
				}


				seed_id1++;

			}
		read1_end:



			/*************************************第二个read的第一个种子********************************************/

			/*******控制种子数量************************/
			max_seed_number2 = read_batch2[i].length / 10 - 1;
			if (max_seed_number2 > max_max_seed_number2)
			{
				max_seed_number2 = max_max_seed_number2;
			}
			/*******控制种子数量************************/

			unique_best_read2 = -1;

			direct_jump2 = 0;

			best_mapp_occ2 = 0;


			C_to_T_forward(read_batch2[i].seq, bsSeq2, read_batch2[i].length, &C_site2);

			total_match_length2 = 0;

			seed_id2 = 0;

			full_seed_id2 = 0;


			locates2 = candidates2;

			candidate_length2 = 0;

			is_mutiple_map2 = 0;

			///这里要改
			get_error2 = -1;
			///这里要改
			extra_seed_flag2 = 1;


			////5'端第一个seed单独处理
			if (seed_id2 < max_seed_number2 && total_match_length2 < read_batch2[i].length)
			{
				current_seed_length2 = read_batch2[i].length - total_match_length2;


				number_of_hits2 =
					count_backward_as_much_1_terminate
					(bsSeq2, current_seed_length2, &top, &bot, &pre_top, &pre_bot, &match_length2);

				///这里要改
				first_seed_match_length2 = match_length2;


				/**
				有两种情况能直接跳转
				1. 找到一个精确的唯一匹配
				2. 找到一堆精确的唯一匹配
				这两种情况best_mapp_occ2都不是0
				所以reseed不用担心这两种情况
				**/
				////新加的
				read2_seed_start[full_seed_id2] = total_match_length2;
				read2_seed_length[full_seed_id2] = match_length2;
				full_seed_id2++;
				////新加的




				///这里要改
				if (number_of_hits2 == 1 && try_process_unique_mismatch_end_to_end_return_site(
					read_batch2[i].length, read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual, cigar,
					read_batch2[i].name,
					locates2,
					top, tmp_ref, &match_length2, read_batch2[i].length - C_site2 - 1,
					&candidates_votes2[0].site, &candidates_votes2[0].err, &get_error2))
				{
					candidates_votes2[0].end_site = read_batch2[i].length - 1;
					best_mapp_occ2 = 1;
					///唯一匹配
					return_flag2 = 1;
					direct_jump2 = 1;
					read_2_jump = 1;
					goto read2_end;
				}


				///只有第一个seed才有可能出现这种情况
				if (match_length2 == read_batch2[i].length
					&&
					number_of_hits2 > 1
					&&
					number_of_hits2 <= max_candidates_occ)
				{
					is_mutiple_map2 = 1;

					///等于-1说明read里面没有C
					if (C_site2 == -1)
					{
						locate_muti_thread(bsSeq2 + current_seed_length2 - match_length2, top, bot, pre_top, pre_bot, locates2, match_length2, &tmp_SA_length, &get_queue);
						reverse_and_adjust_site(locates2, &tmp_SA_length, match_length2, total_match_length2);

						///排序是为了最后好合并双端的结果
						std::sort(locates2, locates2 + tmp_SA_length);

						best_mapp_occ2 = tmp_SA_length;

						for (inner_i = 0; inner_i < tmp_SA_length; inner_i++)
						{
							candidates_votes2[inner_i].site = locates2[inner_i];
							candidates_votes2[inner_i].err = 0;
							candidates_votes2[inner_i].end_site = read_batch2[i].length - 1;
						}




						///多匹配
						return_flag2 = 2;
						direct_jump2 = 1;
						read_2_jump = 1;
						goto read2_end;
					}
				}




				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
				if (number_of_hits2 == 1)
				{

					tmp_SA_length = 1;

					candidate_length2 = candidate_length2 + tmp_SA_length;

					locates2 = candidates2 + candidate_length2;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length2 >= available_seed_length && number_of_hits2 <= max_seed_matches)
				{

					if (number_of_hits2 != 0)
					{

						tmp_SA_length = 0;

						locate_muti_thread(bsSeq2 + current_seed_length2 - match_length2, top, bot, pre_top, pre_bot, locates2, match_length2, &tmp_SA_length, &get_queue);
						reverse_and_adjust_site(locates2, &tmp_SA_length, match_length2, total_match_length2);

						candidate_length2 = candidate_length2 + tmp_SA_length;
						locates2 = candidates2 + candidate_length2;
					}

				}
				else
				{
					////新加的
					full_seed_id2--;
					////新加的
				}







				if (match_length2 == 0)
				{
					total_match_length2 = total_match_length2 + match_step;
				}
				else
				{
					total_match_length2 = total_match_length2 + match_length2 / 2;
				}


				seed_id2++;

			}

		read2_end:

			/*************************************处理第一个种子********************************************/



			///if (number_of_hits1 <= number_of_hits2 && candidate_length1 <= candidate_length2)
			if (candidate_length1 <= candidate_length2)
			{
				/*************************************第一个read********************************************/
				if (read_1_jump == 0)
				{

					best_mapp_occ1 = process_rest_seed_muti_thread(
						seed_id1,
						max_seed_number1,
						total_match_length1,
						&read_batch1[i],
						current_seed_length1,
						number_of_hits1,
						bsSeq1,
						match_length1,
						locates1,
						candidates1,
						candidate_length1,
						available_seed_length,
						max_seed_matches,
						match_step,
						candidates_votes1,
						candidates_votes_length1,
						error_threshold1,
						tmp_ref,
						Peq_SSE,
						min_err,
						&get_queue,
						get_error1,
						first_seed_match_length1,
						read1_seed_start,
						read1_seed_length,
						&full_seed_id1);

				}


				if (best_mapp_occ1 == 0)
				{
					i++;
					continue;
				}




				//*************************************第一个read********************************************






				/*************************************第二个read********************************************/

				if (read_2_jump == 0)
				{

					best_mapp_occ2 = process_rest_seed_filter_muti_thread(
						seed_id2,
						max_seed_number2,
						total_match_length2,
						&read_batch2[i],
						current_seed_length2,
						number_of_hits2,
						bsSeq2,
						match_length2,
						locates2,
						candidates2,
						candidate_length2,
						available_seed_length,
						max_seed_matches,
						match_step,
						candidates_votes2,
						candidates_votes_length2,
						error_threshold2,
						tmp_ref,
						Peq_SSE,
						min_err,
						candidates_votes1,
						best_mapp_occ1,
						inner_maxDistance_pair,
						inner_minDistance_pair,
						&get_queue,
						get_error2,
						first_seed_match_length2,
						read2_seed_start,
						read2_seed_length,
						&full_seed_id2);



				}





				/*************************************第二个read********************************************/
			}
			else
			{


				/*************************************第二个read********************************************/

				if (read_2_jump == 0)
				{

					best_mapp_occ2 = process_rest_seed_muti_thread(
						seed_id2,
						max_seed_number2,
						total_match_length2,
						&read_batch2[i],
						current_seed_length2,
						number_of_hits2,
						bsSeq2,
						match_length2,
						locates2,
						candidates2,
						candidate_length2,
						available_seed_length,
						max_seed_matches,
						match_step,
						candidates_votes2,
						candidates_votes_length2,
						error_threshold2,
						tmp_ref,
						Peq_SSE,
						min_err,
						&get_queue,
						get_error2,
						first_seed_match_length2,
						read2_seed_start,
						read2_seed_length,
						&full_seed_id2);

				}


				if (best_mapp_occ2 == 0)
				{
					i++;
					continue;
				}





				/*************************************第二个read********************************************/












				/*************************************第一个read********************************************/
				if (read_1_jump == 0)
				{

					best_mapp_occ1 = process_rest_seed_filter_muti_thread(
						seed_id1,
						max_seed_number1,
						total_match_length1,
						&read_batch1[i],
						current_seed_length1,
						number_of_hits1,
						bsSeq1,
						match_length1,
						locates1,
						candidates1,
						candidate_length1,
						available_seed_length,
						max_seed_matches,
						match_step,
						candidates_votes1,
						candidates_votes_length1,
						error_threshold1,
						tmp_ref,
						Peq_SSE,
						min_err,
						candidates_votes2,
						best_mapp_occ2,
						inner_maxDistance_pair,
						inner_minDistance_pair,
						&get_queue,
						get_error1,
						first_seed_match_length1,
						read1_seed_start,
						read1_seed_length,
						&full_seed_id1);

				}




				//*************************************第一个read********************************************






			}






			



			///if (best_mapp_occ1 != 0 && best_mapp_occ2 == 0)
			if (best_mapp_occ2 == 0)
			{


				best_mapp_occ2 = reseed_filter_muti_thread(
					seed_id2,
					max_seed_number2,
					total_match_length2,
					&read_batch2[i],
					current_seed_length2,
					number_of_hits2,
					bsSeq2,
					match_length2,
					locates2,
					candidates2,
					candidate_length2,
					available_seed_length,
					max_seed_matches,
					match_step,
					candidates_votes2,
					candidates_votes_length2,
					error_threshold2,
					tmp_ref,
					Peq_SSE,
					min_err,
					candidates_votes1,
					best_mapp_occ1,
					inner_maxDistance_pair,
					inner_minDistance_pair,
					&get_queue,
					get_error2,
					first_seed_match_length2,
					read2_seed_start,
					read2_seed_length,
					full_seed_id2,
					result_seed_start,
					result_seed_length);
			}
			else if (best_mapp_occ1 == 0)  ///if (best_mapp_occ1 == 0 && best_mapp_occ2 != 0)
			{
				best_mapp_occ1 = reseed_filter_muti_thread(
					seed_id1,
					max_seed_number1,
					total_match_length1,
					&read_batch1[i],
					current_seed_length1,
					number_of_hits1,
					bsSeq1,
					match_length1,
					locates1,
					candidates1,
					candidate_length1,
					available_seed_length,
					max_seed_matches,
					match_step,
					candidates_votes1,
					candidates_votes_length1,
					error_threshold1,
					tmp_ref,
					Peq_SSE,
					min_err,
					candidates_votes2,
					best_mapp_occ2,
					inner_maxDistance_pair,
					inner_minDistance_pair,
					&get_queue,
					get_error1,
					first_seed_match_length1,
					read1_seed_start,
					read1_seed_length,
					full_seed_id1,
					result_seed_start,
					result_seed_length);
			}


			
			








			mapping_pair = 0;

			if (unique_best_read1 != -1 && unique_best_read2 != -1)
			{

				if (candidates_votes1[unique_best_read1].site > candidates_votes2[unique_best_read2].site)
				{
					distance = candidates_votes1[unique_best_read1].site - candidates_votes2[unique_best_read2].site;
				}
				else
				{
					distance = candidates_votes2[unique_best_read2].site - candidates_votes1[unique_best_read1].site;
				}

				if (distance <= inner_maxDistance_pair
					&&
					distance >= inner_minDistance_pair)
				{
					best_pair_1_index = unique_best_read1;
					best_pair_2_index = unique_best_read2;
					mapping_pair = 1;
				}
			}




			if (mapping_pair == 0)
			{

				mapping_pair = new_faster_verify_pairs(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
					candidates_votes1, candidates_votes2,
					&best_pair_1_index, &best_pair_2_index,
					inner_maxDistance_pair, inner_minDistance_pair);
			}


		end_i:
			if (mapping_pair == 1)
			{
				///unique_matched_read++;

				//****************第一个read的后处理
				result1.err = candidates_votes1[best_pair_1_index].err;
				result1.origin_site = candidates_votes1[best_pair_1_index].site;
				result1.end_site = candidates_votes1[best_pair_1_index].end_site;

				if (result1.err != 0)
				{
					calculate_best_map_cigar_end_to_end_return
						(&min_candidates_votes_length,
						read_batch1[i].length, error_threshold1, tmp_ref,
						read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual, cigar, path, matrix, matrix_bit,
						read_batch1[i].name,
						result1.cigar,
						&result1, &matched_length1);



				}
				else
				{

					output_sam_end_to_end_return(
						result1.origin_site,
						result1.end_site,
						result1.end_site + 1 - read_batch1[i].length,
						read_batch1[i].length,
						&result1
						);

					sprintf(result1.cigar, "%dM", read_batch1[i].length);
					matched_length1 = read_batch1[i].length;
				}
				//****************第一个read的后处理




				//****************第二个read的后处理
				result2.err = candidates_votes2[best_pair_2_index].err;
				result2.origin_site = candidates_votes2[best_pair_2_index].site;
				result2.end_site = candidates_votes2[best_pair_2_index].end_site;

				if (result2.err != 0)
				{
					calculate_best_map_cigar_end_to_end_return
						(&min_candidates_votes_length,
						read_batch2[i].length, error_threshold2, tmp_ref,
						read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual, cigar, path, matrix, matrix_bit,
						read_batch2[i].name,
						result2.cigar,
						&result2, &matched_length2);
				}
				else
				{

					output_sam_end_to_end_return(
						result2.origin_site,
						result2.end_site,
						result2.end_site + 1 - read_batch2[i].length,
						read_batch2[i].length,
						&result2
						);

					sprintf(result2.cigar, "%dM", read_batch2[i].length);
					matched_length2 = read_batch2[i].length;
				}
				//****************第二个read的后处理



				if (result2.site >= result1.site)
				{
					paired_end_distance = result2.site - result1.site + matched_length2;

				}
				else
				{
					paired_end_distance = result1.site - result2.site + matched_length1;
				}

				if (paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
				{
					unique_matched_read++;
					/**
					directly_output_read1_return_result(read_batch1[i].name, read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual,
						&result1, &result2, read_batch1[i].length,
						matched_length1, matched_length2, result_string);

					directly_output_read2_return_result(read_batch2[i].name, read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual,
						&result2, &result1, read_batch2[i].length,
						matched_length2, matched_length1, result_string);
					**/

					directly_output_read1_return_buffer(read_batch1[i].name, read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual,
						&result1, &result2, read_batch1[i].length,
						matched_length1, matched_length2, &current_sub_buffer);

					directly_output_read2_return_buffer(read_batch2[i].name, read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual,
						&result2, &result1, read_batch2[i].length,
						matched_length2, matched_length1, &current_sub_buffer);

				}



			}
			else if (mapping_pair > 1)
			{
				ambious_matched_read++;
			}


			i++;
		}



		current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

		push_results_to_buffer(&current_sub_buffer);

		current_sub_buffer.length = 0;

	}


	finish_output_buffer();

	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;

	fprintf(stdout, "thread %d completed!\n", thread_id);


}





int Map_Single_Seq_end_to_end(int thread_id)
{
	long long enq_i = 0;

	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;

	Read current_read;

	init_single_read(&current_read);


	bitmapper_bs_iter number_of_hits = 0;

	unsigned int bsSeq_max_length = 1000;
	char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;
	///bitmapper_bs_iter available_seed_length = 30;
	
	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;

	available_seed_length = over_all_seed_length;
	
	/**
	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	**/

	bitmapper_bs_iter current_seed_length;

	bitmapper_bs_iter* locates;
	bitmapper_bs_iter* candidates;
	bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	seed_votes* candidates_votes;
	bitmapper_bs_iter candidates_votes_length;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter seed_id;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;
	///int is_duplicate;
	///bitmapper_bs_iter min_index;

	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
		error_threshold = thread_e;
	}
	**/
	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;



	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	///max_seed_number = 9;
	max_seed_number = 25;
	int max_max_seed_number = max_seed_number;


	candidates = (bitmapper_bs_iter *)malloc(max_seed_number * max_seed_matches * sizeof(bitmapper_bs_iter));
	candidates_votes = (seed_votes *)malloc(max_seed_number * max_seed_matches * sizeof(seed_votes));

	bitmapper_bs_iter* re_seed_candidates
		= (bitmapper_bs_iter *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(bitmapper_bs_iter));

	seed_votes* re_seed_candidates_votes
		= (seed_votes *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(seed_votes));


	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(max_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int C_site_array[SEQ_MAX_LENGTH];
	int C_site_number;

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	fprintf(stdout, "Welcome to BitMapperBS!\n");

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}


	int is_mutiple_map = 0;

	int C_site;


	int return_flag;

	double startTime = Get_T();

	int file_flag;

	int get_error;

	//long long debug_0_mismatch = 0;
	long long debug_1_mismatch = 0;
	//long long debug_seed_match_length = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;
	///int processed = 0;

	//正向模式
	i = 0;
	while (1)
	{
		file_flag = inputReads_single_directly(&current_read);
		///如果等于0, 估计说明read读完了吧
		if (file_flag == 0)
		{
			break;
		}


		enq_i++;


		///等于3应该是有过多的N吧
		if (file_flag == 3)
		{
			continue;
		}


	    C_to_T_forward(current_read.seq, bsSeq, current_read.length, &C_site);


		error_threshold1 = thread_e_f * current_read.length;
		if (error_threshold1 >= max_error_threshold)
		{
			error_threshold1 = max_error_threshold;
		}



		total_match_length = 0;

		seed_id = 0;


		locates = candidates;

		candidate_length = 0;

		is_mutiple_map = 0;

		get_error = -1;

		extra_seed_flag = 1;

		/*******控制种子数量************************/
		max_seed_number = current_read.length / 10 - 1;
		///多用一个不重叠的种子
		///max_seed_number = current_read.length / 10;
		if (max_seed_number > max_max_seed_number)
		{
			max_seed_number = max_max_seed_number;
		}
		/*******控制种子数量************************/


		////5'端第一个seed单独处理
		if (seed_id < max_seed_number && total_match_length < current_read.length)
		{
			current_seed_length = current_read.length - total_match_length;

			number_of_hits =
				count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);
			

			first_seed_match_length = match_length;

			
			/**
			if (number_of_hits == 1 && try_process_unique_exact_match_end_to_end(
				current_read.length, current_read.seq, current_read.rseq, current_read.qual, cigar, current_read.name,
				locates,
				top, tmp_ref, &match_length, current_read.length - C_site - 1))**/
			if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end(
				current_read.length, current_read.seq, current_read.rseq, current_read.qual, cigar, current_read.name,
				locates,
				top, tmp_ref, &match_length, current_read.length - C_site - 1, &get_error))
			{

				///debug_0_mismatch++;

				unique_matched_read++;
				matched_read++;
				continue;
			}


			///只有第一个seed才有可能出现这种情况
			if (match_length == current_read.length
				&&
				number_of_hits > 1)
			{
				is_mutiple_map = 1;

				///等于-1说明read里面没有C
				if (C_site == -1)
				{
					matched_read++;
					i++;
					continue;
				}


			}



			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
			if (number_of_hits == 1)
			{
				tmp_SA_length = 1;

				locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

				candidate_length = candidate_length + tmp_SA_length;

				locates = candidates + candidate_length;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
			{
				if (number_of_hits != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length);
					reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;
					locates = candidates + candidate_length;

				}





			}



			if (match_length == 0)
			{
				total_match_length = total_match_length + match_step;
			}
			else
			{
				total_match_length = total_match_length + match_length / 2;
			}



			seed_id++;

		}




		
		///说明是目前的最优解是1-mismatch
		if (get_error == 1)
		{

			second_seed_length = current_read.length - first_seed_match_length;





			if (second_seed_length >= 17)
			{

				number_of_hits = count_hash_table(bsSeq, second_seed_length, &top, &bot, &pre_top, &pre_bot);
		
				if (number_of_hits == 1)
				{
					tmp_SA_length = 0;

					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, 
						second_seed_length, first_seed_match_length);

					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;

					extra_seed_flag = 0;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (number_of_hits <= max_seed_matches)
				{

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate(bsSeq, top, bot, pre_top, pre_bot, locates, second_seed_length, &tmp_SA_length);

						reverse_and_adjust_site(locates, &tmp_SA_length, second_seed_length, first_seed_match_length);


						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}

					extra_seed_flag = 0;

				}
				else if (number_of_hits > max_seed_matches)
				{
					extra_seed_flag = 1;
				}


			}
			else
			{
				extra_seed_flag = 1;
			}


		}




		

		if (extra_seed_flag == 1)
		{

		

			///这是剩下的seed
			while (seed_id < max_seed_number && total_match_length < current_read.length)
			{


				current_seed_length = current_read.length - total_match_length;

				number_of_hits =
					count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);





				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				if (number_of_hits == 1)
				{
					tmp_SA_length = 0;
					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, match_length, total_match_length);

					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
				{

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length);
						reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);
						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}





				}
				else if (current_seed_length == match_length)
				{
					break;
				}
			
			
		
				if (match_length == 0)
				{
					total_match_length = total_match_length + match_step;
				}
				else
				{
					total_match_length = total_match_length + match_length / 2;
				}

				seed_id++;

			}




		}

		///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
		if (extra_seed_flag == 0 
			&& 
			(candidate_length == 1 || 
			(candidate_length == 2 && candidates[0] == candidates[1])
			))
		{
			sprintf(cigar, "%dM", current_read.length);

			output_sam_end_to_end
				(
				current_read.name, 
				current_read.seq, 
				current_read.rseq, 
				current_read.qual,
				candidates[0],
				current_read.length - 1,
				0,
				1,
				cigar, 
				current_read.length);

			unique_matched_read++;
			matched_read++;

			debug_1_mismatch++;
		}
		else if (candidate_length != 0)
		{

			std::sort(candidates, candidates + candidate_length);
			////generate_candidate_votes_allow_error(candidates, candidate_length, candidates_votes, &candidates_votes_length);

			///generate_candidate_votes(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold);
			generate_candidate_votes_shift(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold1);

			///这个排序其实有助于剪枝
			///其实好像也没什么用..
			std::sort(candidates_votes, candidates_votes + candidates_votes_length, compare_seed_votes);


			///感觉我这样做的话，不可能出现有单独一个精确匹配的位置的情况了
			///如果一个read有唯一的精确匹配位置，那么他肯定在前面就被报告出来了
			if (is_mutiple_map == 0)
			{

				
				if (error_threshold1 <= 15)
				{
					map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);

					
				}
				else
				{
					map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);
				}
					
				
			}
			else
			{

				if (error_threshold1 <= 15)
				{
					map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);
					
				}
				else
				{
					map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);
				}
				
				
			}


			if (min_err_index >= 0)
			{

				matched_read++;
				calculate_best_map_cigar_end_to_end
					(candidates_votes, &min_candidates_votes_length,
					current_read.length, error_threshold1, tmp_ref,
					current_read.seq, current_read.rseq, current_read.qual, cigar, path, matrix, matrix_bit, current_read.name,
					best_cigar);

				unique_matched_read++;
			}
			else if (min_err_index != -1)
			{
				matched_read++;
			}

		}
	}


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;

	/**
	fprintf(stderr, "debug_0_mismatch: %lld\n", debug_0_mismatch);
	fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);
	fprintf(stderr, "debug_seed_match_length: %lld\n", debug_seed_match_length);
	**/
	
	fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);

	return 1;
}



int Map_Single_Seq_end_to_end_pbat(int thread_id)
{
	long long enq_i = 0;
	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;



	cand = (Candidate*)malloc(sizeof(Candidate)*(SEQ_LENGTH - WINDOW_SIZE + 1));
	vertify = (Candidate*)malloc(sizeof(Candidate)*((thread_e + 1 + 1)));
	all_vertify = (Candidate*)malloc(sizeof(Candidate)*(MERGE_STEP)*((thread_e + 1 + 1)));
	candidates_locs = (unsigned int**)malloc((thread_e + 1 + 1)*(MERGE_STEP)*sizeof(unsigned int*));
	listsites = (unsigned int*)malloc((thread_e + 1 + 1)*(MERGE_STEP)*sizeof(unsigned int));
	each_list_length = (unsigned int*)malloc((thread_e + 1 + 1)*(MERGE_STEP)*sizeof(unsigned int));


	char** ref = (char**)malloc(sizeof(char*)* 8);

	for (i = 0; i<8; i++)
	{
		ref[i] = (char*)malloc(sizeof(char)*(SEQ_LENGTH + 2 * thread_e));
	}


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;




	Read current_read;

	init_single_read(&current_read);





	bitmapper_bs_iter number_of_hits = 0;

	unsigned int bsSeq_max_length = 1000;
	char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;
	///bitmapper_bs_iter available_seed_length = 30;

	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;

	available_seed_length = over_all_seed_length;

	/**
	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;
	**/

	bitmapper_bs_iter current_seed_length;

	bitmapper_bs_iter* locates;
	bitmapper_bs_iter* candidates;
	bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	seed_votes* candidates_votes;
	bitmapper_bs_iter candidates_votes_length;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter seed_id;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;
	///int is_duplicate;
	///bitmapper_bs_iter min_index;

	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
		error_threshold = thread_e;
	}
	**/
	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;

	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	///max_seed_number = 9;
	max_seed_number = 25;
	int max_max_seed_number = max_seed_number;


	candidates = (bitmapper_bs_iter *)malloc(max_seed_number * max_seed_matches * sizeof(bitmapper_bs_iter));
	candidates_votes = (seed_votes *)malloc(max_seed_number * max_seed_matches * sizeof(seed_votes));

	bitmapper_bs_iter* re_seed_candidates
		= (bitmapper_bs_iter *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(bitmapper_bs_iter));

	seed_votes* re_seed_candidates_votes
		= (seed_votes *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(seed_votes));


	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(max_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	double re_seed_threshold = 0.80;
	int reseed_length;
	bitmapper_bs_iter re_seed_error_threshold;

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	long long total_number_of_hit = 0;




	int is_mutiple_map = 0;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	long long direct_cut = 0;
	///long long number_of_short_seeds = 0;
	///long long number_of_short_seeds_matches = 0;
	int C_site;

	///double total_c_time = 0;
	///double start_c_time;


	double startTime = Get_T();

	int file_flag;

	int get_error;
	long long debug_1_mismatch = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;



	//正向模式
	i = 0;
	///while (getReads_single(&current_read))
	while (1)
	{

		///int file_flag = inputReads_single(&current_read);

		file_flag = inputReads_single_directly_pbat(&current_read);

		

		///如果等于0, 估计说明read读完了吧
		if (file_flag == 0)
		{
			break;
		}

		enq_i++;


		///等于3应该是有过多的N吧
		if (file_flag == 3)
		{
			continue;
		}



		/**
		if (i % 100000 == 0)
		{
			fprintf(stdout, "i: %llu \n", i);

		}
		**/

		///fprintf(stderr, "*******\current_read.length: %llu\n", current_read.length);

		///start_c_time = Get_T();
		///因为参考组是反的, 那么查询的时候实际上read也要反着查
		///所以实际上对是把read的反向做bs转换
		///C_to_T(current_read.rseq, bsSeq, current_read.length, &C_site);
		C_to_T_forward(current_read.seq, bsSeq, current_read.length, &C_site);



		error_threshold1 = thread_e_f * current_read.length;
		if (error_threshold1 >= max_error_threshold)
		{
			error_threshold1 = max_error_threshold;
		}

		///total_c_time = total_c_time + Get_T() - start_c_time;


		total_match_length = 0;

		seed_id = 0;


		locates = candidates;

		candidate_length = 0;

		is_mutiple_map = 0;

		get_error = -1;

		extra_seed_flag = 1;



		/*******控制种子数量************************/
		max_seed_number = current_read.length / 10 - 1;
		if (max_seed_number > max_max_seed_number)
		{
			max_seed_number = max_max_seed_number;
		}
		/*******控制种子数量************************/




		////5'端第一个seed单独处理
		if (seed_id < max_seed_number && total_match_length < current_read.length)
		{
			current_seed_length = current_read.length - total_match_length;

			number_of_hits =
				count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);


			first_seed_match_length = match_length;



















			/**
			if (number_of_hits == 1 && try_process_unique_exact_match_end_to_end_pbat(
				current_read.length, current_read.seq, current_read.rseq, current_read.qual, cigar, current_read.name,
				locates,
				top, tmp_ref, &match_length, current_read.length - C_site - 1))**/
			if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end_pbat(
				current_read.length, current_read.seq, current_read.rseq, current_read.qual, cigar, current_read.name,
				locates,
				top, tmp_ref, &match_length, current_read.length - C_site - 1, &get_error))
			{

				unique_matched_read++;
				matched_read++;
				i++;
				continue;


			}



			











			///fprintf(stderr, "seed_id: %llu\n", seed_id);
			///fprintf(stderr, "match_length: %llu\n", match_length);

			///只有第一个seed才有可能出现这种情况
			if (match_length == current_read.length
				&&
				number_of_hits > 1)
			{
				is_mutiple_map = 1;

				///等于-1说明read里面没有C
				if (C_site == -1)
				{
					matched_read++;
					i++;
					direct_cut++;



					continue;
				}
			}




			///如果当前seed的只有一个匹配位置
			///则不管多长都是有效的
			///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
			if (number_of_hits == 1)
			{
				total_seed_matches = total_seed_matches + number_of_hits;
				total_seed_number++;

				tmp_SA_length = 1;


				locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

				candidate_length = candidate_length + tmp_SA_length;

				locates = candidates + candidate_length;


			}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
			else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
			{

				total_seed_matches = total_seed_matches + number_of_hits;
				total_seed_number++;

				if (number_of_hits != 0)
				{

					tmp_SA_length = 0;

					locate(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length);
					reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);
					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;
					locates = candidates + candidate_length;

				}





			}


			///total_match_length = total_match_length + match_step;

			if (match_length == 0)
			{
				total_match_length = total_match_length + match_step;
			}
			else
			{
				total_match_length = total_match_length + match_length / 2;
			}



			seed_id++;

		}

		








		///说明是目前的最优解是1-mismatch
		if (get_error == 1)
		{

			second_seed_length = current_read.length - first_seed_match_length;





			if (second_seed_length >= 17)
			{

				number_of_hits = count_hash_table(bsSeq, second_seed_length, &top, &bot, &pre_top, &pre_bot);

				if (number_of_hits == 1)
				{
					tmp_SA_length = 0;

					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length,
						second_seed_length, first_seed_match_length);

					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;

					extra_seed_flag = 0;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (number_of_hits <= max_seed_matches)
				{

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate(bsSeq, top, bot, pre_top, pre_bot, locates, second_seed_length, &tmp_SA_length);

						reverse_and_adjust_site(locates, &tmp_SA_length, second_seed_length, first_seed_match_length);


						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}

					extra_seed_flag = 0;

				}
				else if (number_of_hits > max_seed_matches)
				{
					extra_seed_flag = 1;
				}


			}
			else
			{
				extra_seed_flag = 1;
			}


		}
















		if (extra_seed_flag == 1)
		{



			///这是剩下的seed
			while (seed_id < max_seed_number && total_match_length < current_read.length)
			{



				current_seed_length = current_read.length - total_match_length;


				number_of_hits =
					count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);


				///fprintf(stderr, "seed_id: %llu\n", seed_id);
				///fprintf(stderr, "match_length: %llu\n", match_length);


				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				if (number_of_hits == 1)
				{
					total_seed_matches = total_seed_matches + number_of_hits;
					total_seed_number++;



					tmp_SA_length = 0;

					/**
					locate_one_position(locates, top, &tmp_SA_length);

					reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);
					**/

					locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, match_length, total_match_length);


					/**
					tmp_SA_length = 1;

					fprintf(stderr, "*******match_length: %llu\n", match_length);

					try_process_unique_exact_match_end_to_end_second_seed(
					current_seed_length, current_read.seq + total_match_length,
					locates, top, tmp_ref, &match_length);
					fprintf(stderr, "match_length: %llu\n", match_length);
					**/

					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
				{

					total_seed_matches = total_seed_matches + number_of_hits;
					total_seed_number++;

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length);
						reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);
						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}





				}
				else if (current_seed_length == match_length)
				{
					break;
				}






				///total_match_length = total_match_length + match_step;

				if (match_length == 0)
				{
					total_match_length = total_match_length + match_step;
				}
				else
				{
					total_match_length = total_match_length + match_length / 2;
				}

				seed_id++;


			}
		}

        ///fprintf(stderr, "candidate_length: %lld\n", candidate_length);









		///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
		if (extra_seed_flag == 0
			&&
			(candidate_length == 1 ||
			(candidate_length == 2 && candidates[0] == candidates[1])
			))
		{
			sprintf(cigar, "%dM", current_read.length);

			output_sam_end_to_end_pbat
				(
				current_read.name,
				current_read.seq,
				current_read.rseq,
				current_read.qual,
				candidates[0],
				current_read.length - 1,
				0,
				1,
				cigar,
				current_read.length);

			unique_matched_read++;
			matched_read++;

			debug_1_mismatch++;
		}
		else if(candidate_length != 0)
		{

			std::sort(candidates, candidates + candidate_length);
			////generate_candidate_votes_allow_error(candidates, candidate_length, candidates_votes, &candidates_votes_length);

			///generate_candidate_votes(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold);
			generate_candidate_votes_shift(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold1);

			total_candidate_length = total_candidate_length + candidate_length;
			total_candidates_votes_length = total_candidates_votes_length + candidates_votes_length;

			///fprintf(stderr, "candidates_votes_length: %lld\n", candidates_votes_length);
			///fprintf(stderr, "is_mutiple_map: %lld\n", is_mutiple_map);

			///这个排序其实有助于剪枝
			///其实好像也没什么用..
			std::sort(candidates_votes, candidates_votes + candidates_votes_length, compare_seed_votes);
			/**
			if (i == 12427839)
			{
			for (size_t debug_i = 0; debug_i < candidates_votes_length; debug_i++)
			{
			fprintf(stderr, "candidates_votes[%u].site: %llu\n", debug_i, candidates_votes[debug_i].site);
			fflush(stderr);

			}
			}
			**/

			///感觉我这样做的话，不可能出现有单独一个精确匹配的位置的情况了
			///如果一个read有唯一的精确匹配位置，那么他肯定在前面就被报告出来了
			if (is_mutiple_map == 0)
			{


				if (error_threshold1 <= 15)
				{
					map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);


				}
				else
				{
					map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);
				}


			}
			else
			{

				if (error_threshold1 <= 15)
				{
					map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);

				}
				else
				{
					map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index);
				}


			}





			min_candidates_votes_length = 0;

			if (min_err_index >= 0)
			{

				total_best_mapping_site++;
				matched_read++;
				min_candidates_votes_length = 0;

				calculate_best_map_cigar_end_to_end_pbat
					(candidates_votes, &min_candidates_votes_length,
					current_read.length, error_threshold1, tmp_ref,
					current_read.seq, current_read.rseq, current_read.qual, cigar, path, matrix, matrix_bit, current_read.name,
					best_cigar);

				unique_matched_read++;


			}
			else if (min_err_index == -1)
			{

				///fprintf(stderr, "name: %s\n", current_read.name);
				///fprintf(stderr, "length: %u\n", current_read.length);
				///fprintf(stderr, "max_seed_number: %u\n", max_seed_number);
				///fprintf(stderr, "seq: %s\n", current_read.seq);
				///fprintf(stderr, "bsSeq: %s\n", bsSeq);
				///fprintf(stderr, "rseq: %s\n\n\n", current_read.rseq);
				
				

				unmatched_read++;
			}
			else
			{
				matched_read++;
			}


		}
		else
		{
			empty_read++;
		}


		i++;
	}


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;


	fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);


	return 1;
}


int Map_Pair_Seq(int thread_id)
{
	if (is_local == 0)
	{
		Map_Pair_Seq_end_to_end(thread_id);
	}
	else
	{
		Map_Pair_Seq_end_to_end_fast(thread_id);
	}
	
}



int Map_Single_Seq(int thread_id)
{
	Map_Single_Seq_end_to_end(thread_id);
}


int Map_Pair_Seq_muti_thread(int thread_id)
{
	///build_back_up_set();

	init_Pair_Seq_input_buffer(THREAD_COUNT);
	init_output_buffer(THREAD_COUNT);

	pthread_t inputReadsHandle;
	pthread_create(&inputReadsHandle, NULL, input_pe_reads_muti_threads, NULL);



	pthread_t outputResultSinkHandle;
	pthread_create(&outputResultSinkHandle, NULL, pop_buffer, NULL);




	int ijkijk = 0;


	pthread_t *_r_threads;


	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*THREAD_COUNT);


	int i = 0;

	for (i = 0; i < THREAD_COUNT; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;

		if (is_local == 0)
		{
			pthread_create(_r_threads + i, NULL, Map_Pair_Seq_split, (void*)arg);
		}
		else
		{
			pthread_create(_r_threads + i, NULL, Map_Pair_Seq_split_fast, (void*)arg);
		}
		///pthread_create(_r_threads + i, NULL, Map_Single_Seq_split_one_by_one, (void*)arg);

	}

	pthread_join(inputReadsHandle, NULL);

	for (i = 0; i<THREAD_COUNT; i++)
		pthread_join(_r_threads[i], NULL);

	pthread_join(outputResultSinkHandle, NULL);

	

	free(_r_threads);

	
}


int Map_Single_Seq_muti_thread(int thread_id)
{

	init_Single_Seq_input_buffer(THREAD_COUNT);
	init_output_buffer(THREAD_COUNT);

	pthread_t inputReadsHandle;
	pthread_create(&inputReadsHandle, NULL, input_single_reads_muti_threads, NULL);



	pthread_t outputResultSinkHandle;
	pthread_create(&outputResultSinkHandle, NULL, pop_buffer, NULL);




	int ijkijk = 0;


	pthread_t *_r_threads;


	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*THREAD_COUNT);


	int i = 0;

	for (i = 0; i < THREAD_COUNT; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;


		pthread_create(_r_threads + i, NULL, Map_Single_Seq_split, (void*)arg);

		///pthread_create(_r_threads + i, NULL, Map_Single_Seq_split_one_by_one, (void*)arg);

	}

	pthread_join(inputReadsHandle, NULL);

	for (i = 0; i<THREAD_COUNT; i++)
		pthread_join(_r_threads[i], NULL);

	pthread_join(outputResultSinkHandle, NULL);



	free(_r_threads);


}


int Map_Single_Seq_pbat_muti_thread(int thread_id)
{

	init_Single_Seq_input_buffer(THREAD_COUNT);
	init_output_buffer(THREAD_COUNT);

	pthread_t inputReadsHandle;
	pthread_create(&inputReadsHandle, NULL, input_single_reads_muti_threads_pbat, NULL);



	pthread_t outputResultSinkHandle;
	pthread_create(&outputResultSinkHandle, NULL, pop_buffer, NULL);




	int ijkijk = 0;


	pthread_t *_r_threads;


	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*THREAD_COUNT);


	int i = 0;

	for (i = 0; i < THREAD_COUNT; i++)
	{
		int *arg = (int*)malloc(sizeof(*arg));
		*arg = i;


		pthread_create(_r_threads + i, NULL, Map_Single_Seq_split_pbat, (void*)arg);

		///pthread_create(_r_threads + i, NULL, Map_Single_Seq_split_one_by_one, (void*)arg);

	}

	pthread_join(inputReadsHandle, NULL);

	for (i = 0; i<THREAD_COUNT; i++)
		pthread_join(_r_threads[i], NULL);

	pthread_join(outputResultSinkHandle, NULL);



	free(_r_threads);


}


int Map_Single_Seq_pbat(int thread_id)
{

	Map_Single_Seq_end_to_end_pbat(thread_id);
	
}


inline void load_batch_read(Read* read_batch, int batch_read_size, 
	int* return_file_flag, int* obtain_reads_num)
{
	int inner_i = 0;
	int file_flag = 1;

	while (inner_i<batch_read_size)
	{
		file_flag = inputReads_single_directly(&read_batch[inner_i]);

		if (file_flag == 1)
		{
			inner_i++;
		}
		else if (file_flag == 0)
		{
			break;
		}
	}

	*return_file_flag = file_flag;
	*obtain_reads_num = inner_i;

}


inline void load_batch_read_pbat(Read* read_batch, int batch_read_size,
	int* return_file_flag, int* obtain_reads_num)
{
	int inner_i = 0;
	int file_flag = 1;

	while (inner_i<batch_read_size)
	{
		///file_flag = inputReads_single_directly(&read_batch[inner_i]);
		file_flag = inputReads_single_directly_pbat(&read_batch[inner_i]);


		if (file_flag == 1)
		{
			inner_i++;
		}
		else if (file_flag == 0)
		{
			break;
		}
	}

	*return_file_flag = file_flag;
	*obtain_reads_num = inner_i;

}


void* Map_Single_Seq_split(void* arg)
{
	int thread_id = *((int *)arg);

	FILE* output_file = get_Ouput_Dec();

	

	long long enq_i = 0;
	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;

	bwt_locate_queue get_queue;

	init_locate_queue_muti_thread(&get_queue);


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;

	///Read* current_read;

	bitmapper_bs_iter number_of_hits = 0;

	unsigned int bsSeq_max_length = 1000;
	char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;

	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;


	available_seed_length = over_all_seed_length;


	bitmapper_bs_iter current_seed_length;

	bitmapper_bs_iter* locates;
	bitmapper_bs_iter* candidates;
	bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	seed_votes* candidates_votes;
	bitmapper_bs_iter candidates_votes_length;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter seed_id;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;

	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
		error_threshold = thread_e;
	}
	**/


	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;



	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	max_seed_number = 25;
	int max_max_seed_number = max_seed_number;


	candidates = (bitmapper_bs_iter *)malloc(max_seed_number * max_seed_matches * sizeof(bitmapper_bs_iter));
	candidates_votes = (seed_votes *)malloc(max_seed_number * max_seed_matches * sizeof(seed_votes));

	bitmapper_bs_iter* re_seed_candidates
		= (bitmapper_bs_iter *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(bitmapper_bs_iter));

	seed_votes* re_seed_candidates_votes
		= (seed_votes *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(seed_votes));


	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(max_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int C_site_array[SEQ_MAX_LENGTH];
	int C_site_number;

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	double re_seed_threshold = 0.80;
	int reseed_length;
	bitmapper_bs_iter re_seed_error_threshold;

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	long long total_number_of_hit = 0;

	/**
	fprintf(stdout, "error_threshold: %llu\n", error_threshold);
	fprintf(stdout, "score_threshold: %f\n", score_threshold);
	fprintf(stdout, "edit_distance_threshold: %f\n", edit_distance_threshold);
	fprintf(stdout, "available_seed_length: %llu\n", available_seed_length);
	**/

	int is_mutiple_map = 0;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	long long direct_cut = 0;
	long long direct_match = 0;
	int C_site;
	///std::stringstream result_string;
	///std::string tmp_result_string;
	///char* tmp_result_string_char_array = NULL;
	///long long tmp_result_string_char_array_length;

	///batch_read_size

	///Read* read_batch = (Read*)malloc(sizeof(Read)*batch_read_size);
	Read* read_batch;

	Read_buffer_single_sub_block curr_sub_block;

	init_single_sub_block_single(&curr_sub_block);

	Output_buffer_sub_block current_sub_buffer;
	init_buffer_sub_block(&current_sub_buffer);



	///int inner_i;

	double startTime = Get_T();

	int file_flag = 1;
	int obtain_reads_num;


	int get_error;

	long long debug_1_mismatch = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;



	file_flag = 1;
	while (file_flag != 0)
	{


		///file_flag = inputReads_single_directly(&current_read);

		///batch_read_size



		///file_flag = get_reads_mul_thread(&current_read);
		/**
		pthread_mutex_lock(&readinputMutex);


		load_batch_read(read_batch, batch_read_size, &file_flag, &obtain_reads_num);

		pthread_mutex_unlock(&readinputMutex);
		**/


		file_flag = get_single_reads_mul_thread(&curr_sub_block);
		obtain_reads_num = curr_sub_block.sub_block_read_number;
		read_batch = curr_sub_block.read;



		
		
		enq_i = enq_i + obtain_reads_num;

		i = 0;
		while (i < obtain_reads_num)
		{
			
			///因为参考组是反的, 那么查询的时候实际上read也要反着查
			///所以实际上对是把read的反向做bs转换
			///C_to_T(current_read.rseq, bsSeq, current_read.length, &C_site);
			C_to_T_forward(read_batch[i].seq, bsSeq, read_batch[i].length, &C_site);


			////fprintf(stderr, "i: %llu, thread_id: %d, obtain_reads_num: %d\n", thread_id, obtain_reads_num);



			error_threshold1 = thread_e_f * read_batch[i].length;
			if (error_threshold1 >= max_error_threshold)
			{
				error_threshold1 = max_error_threshold;
			}



			total_match_length = 0;

			seed_id = 0;


			locates = candidates;

			candidate_length = 0;

			is_mutiple_map = 0;

			///这里要改
			get_error = -1;
			///这里要改
			extra_seed_flag = 1;

			/*******控制种子数量************************/
			max_seed_number = read_batch[i].length / 10 - 1;
			if (max_seed_number > max_max_seed_number)
			{
				max_seed_number = max_max_seed_number;
			}
			/*******控制种子数量************************/


			////5'端第一个seed单独处理
			if (seed_id < max_seed_number && total_match_length < read_batch[i].length)
			{
				current_seed_length = read_batch[i].length - total_match_length;

				number_of_hits =
					count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);


				///这里要改
				first_seed_match_length = match_length;








				


					/**
				if (number_of_hits == 1 && try_process_unique_exact_match_end_to_end_get_result_string(
					read_batch[i].length, read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual, cigar, read_batch[i].name,
					locates,
					top, tmp_ref, &match_length, read_batch[i].length - C_site - 1, result_string))**/
				/**
				if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end_get_result_string(
					read_batch[i].length, read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual, cigar, read_batch[i].name,
					locates,
					top, tmp_ref, &match_length, read_batch[i].length - C_site - 1, result_string, &get_error))**/
				if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end_output_buffer(
					read_batch[i].length, read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual, cigar, read_batch[i].name,
					locates,
					top, tmp_ref, &match_length, read_batch[i].length - C_site - 1, &current_sub_buffer, &get_error))
				{
					direct_match++;

					unique_matched_read++;
					matched_read++;
					i++;
					continue;


				}





				///只有第一个seed才有可能出现这种情况
				if (match_length == read_batch[i].length
					&&
					number_of_hits > 1)
				{
					is_mutiple_map = 1;

					///等于-1说明read里面没有C
					if (C_site == -1)
					{
						matched_read++;
						i++;
						direct_cut++;
						continue;
					}


				}




				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
				if (number_of_hits == 1)
				{
					total_seed_matches = total_seed_matches + number_of_hits;
					total_seed_number++;

					tmp_SA_length = 1;


					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
				{

					total_seed_matches = total_seed_matches + number_of_hits;
					total_seed_number++;

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate_muti_thread(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length, &get_queue);
						reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}





				}



				if (match_length == 0)
				{
					total_match_length = total_match_length + match_step;
				}
				else
				{
					total_match_length = total_match_length + match_length / 2;
				}






				seed_id++;

			}






			///这里要改
			///说明是目前的最优解是1-mismatch
			if (get_error == 1)
			{

				second_seed_length = read_batch[i].length - first_seed_match_length;





				if (second_seed_length >= 17)
				{

					number_of_hits = count_hash_table(bsSeq, second_seed_length, &top, &bot, &pre_top, &pre_bot);

					if (number_of_hits == 1)
					{
						tmp_SA_length = 0;

						locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length,
							second_seed_length, first_seed_match_length);

						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;

						locates = candidates + candidate_length;

						extra_seed_flag = 0;


					}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
					else if (number_of_hits <= max_seed_matches)
					{

						if (number_of_hits != 0)
						{

							tmp_SA_length = 0;

							locate_muti_thread(bsSeq, top, bot, pre_top, pre_bot, locates, second_seed_length, &tmp_SA_length, &get_queue);

							reverse_and_adjust_site(locates, &tmp_SA_length, second_seed_length, first_seed_match_length);


							locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

							candidate_length = candidate_length + tmp_SA_length;
							locates = candidates + candidate_length;

						}

						extra_seed_flag = 0;

					}
					else if (number_of_hits > max_seed_matches)
					{
						extra_seed_flag = 1;
					}


				}
				else
				{
					extra_seed_flag = 1;
				}


			}































			///这里要改
			if (extra_seed_flag == 1)
			{



				///这是剩下的seed
				while (seed_id < max_seed_number && total_match_length < read_batch[i].length)
				{



					current_seed_length = read_batch[i].length - total_match_length;

					number_of_hits =
						count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);


					///如果当前seed的只有一个匹配位置
					///则不管多长都是有效的
					if (number_of_hits == 1)
					{
						total_seed_matches = total_seed_matches + number_of_hits;
						total_seed_number++;

						tmp_SA_length = 0;
						locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, match_length, total_match_length);

						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;

						locates = candidates + candidate_length;


					}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
					else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
					{

						total_seed_matches = total_seed_matches + number_of_hits;
						total_seed_number++;

						if (number_of_hits != 0)
						{

							tmp_SA_length = 0;

							locate_muti_thread(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length, &get_queue);
							reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);
							locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

							candidate_length = candidate_length + tmp_SA_length;
							locates = candidates + candidate_length;

						}





					}
					else if (current_seed_length == match_length)
					{
						break;
					}



					if (match_length == 0)
					{
						total_match_length = total_match_length + match_step;
					}
					else
					{
						total_match_length = total_match_length + match_length / 2;
					}

					seed_id++;

				}

			}









			///这里要改
			///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
			if (extra_seed_flag == 0
				&&
				(candidate_length == 1 ||
				(candidate_length == 2 && candidates[0] == candidates[1])
				))
			{
				sprintf(cigar, "%dM", read_batch[i].length);
				/**
				output_sam_end_to_end_result_string
					(
					read_batch[i].name,
					read_batch[i].seq,
					read_batch[i].rseq,
					read_batch[i].qual,
					candidates[0],
					read_batch[i].length - 1,
					0,
					1,
					cigar,
					read_batch[i].length, result_string);
				**/

				output_sam_end_to_end_output_buffer(
					read_batch[i].name,
					read_batch[i].seq,
					read_batch[i].rseq,
					read_batch[i].qual,
					candidates[0],
					read_batch[i].length - 1,
					0,
					1,
					cigar,
					read_batch[i].length, &current_sub_buffer);



				unique_matched_read++;
				matched_read++;

				debug_1_mismatch++;
			}
			else if (candidate_length != 0)
			{

				std::sort(candidates, candidates + candidate_length);
				////generate_candidate_votes_allow_error(candidates, candidate_length, candidates_votes, &candidates_votes_length);

				///generate_candidate_votes(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold);
				generate_candidate_votes_shift(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold1);

				total_candidate_length = total_candidate_length + candidate_length;
				total_candidates_votes_length = total_candidates_votes_length + candidates_votes_length;



				///这个排序其实有助于剪枝
				///其实好像也没什么用..
				std::sort(candidates_votes, candidates_votes + candidates_votes_length, compare_seed_votes);


				///感觉我这样做的话，不可能出现有单独一个精确匹配的位置的情况了
				///如果一个read有唯一的精确匹配位置，那么他肯定在前面就被报告出来了
				if (is_mutiple_map == 0)
				{


					if (error_threshold1 <= 15)
					{
						map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);


					}
					else
					{
						map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);
					}


				}
				else
				{
					total_is_mutiple_map++;


					if (error_threshold1 <= 15)
					{
						map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);

					}
					else
					{
						map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);
					}


				}


				min_candidates_votes_length = 0;

				if (min_err_index >= 0)
				{

					total_best_mapping_site++;
					matched_read++;
					min_candidates_votes_length = 0;

					/**
					calculate_best_map_cigar_end_to_end_result_string
						(candidates_votes, &min_candidates_votes_length,
						read_batch[i].length, error_threshold1, tmp_ref,
						read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual,
						cigar, path, matrix, matrix_bit, read_batch[i].name,
						best_cigar, result_string);
					**/

					calculate_best_map_cigar_end_to_end_output_buffer
						(candidates_votes, &min_candidates_votes_length,
						read_batch[i].length, error_threshold1, tmp_ref,
						read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual,
						cigar, path, matrix, matrix_bit, read_batch[i].name,
						best_cigar, &current_sub_buffer);

					unique_matched_read++;


				}
				else if (min_err_index == -1)
				{
					unmatched_read++;
				}
				else
				{
					matched_read++;
				}







			}
			else
			{
				empty_read++;
			}

			i++;

		}

	/**
		pthread_mutex_lock(&queueMutex);

		fprintf(output_file, "%s", result_string.str().c_str());
		///out << result_string.str();

		pthread_mutex_unlock(&queueMutex);

		///stringstream和string都要清空
		result_string.clear();
		result_string.str("");

		
	}
	**/

		current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

		push_results_to_buffer(&current_sub_buffer);

		current_sub_buffer.length = 0;

	}

	finish_output_buffer();



	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;

	fprintf(stdout, "thread %d completed!\n", thread_id);

	///这里要改
	fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);

}


void* Map_Single_Seq_split_pbat(void* arg)
{
	int thread_id = *((int *)arg);

	FILE* output_file = get_Ouput_Dec();



	long long enq_i = 0;
	mappingCnt[thread_id] = 0;
	mappedSeqCnt[thread_id] = 0;
	completedSeqCnt[thread_id] = 0;

	bwt_locate_queue get_queue;

	init_locate_queue_muti_thread(&get_queue);


	unsigned int i;
	unsigned int j;
	int j_m = 0;
	unsigned int tmp_j;
	unsigned int *locs = NULL;
	unsigned int **candidates_locs = NULL;
	unsigned int *listsites = NULL;
	unsigned int *M_lists = NULL;
	unsigned int *swap_list = NULL;
	unsigned int *each_list_length = NULL;
	unsigned int sum_lists_locs;

	Candidate *cand = NULL;
	Candidate *vertify = NULL;
	Candidate *all_vertify = NULL;
	unsigned int length_prefix;
	unsigned int hash_key = 0;
	unsigned int N_site;
	int anti_1 = 1;
	int candidates_length;


	char* tmp_ref = (char*)malloc(sizeof(char)*(SEQ_MAX_LENGTH)* 16);

	memset(tmp_ref, 0, (SEQ_MAX_LENGTH)* 16);


	unsigned int pre_length = 0;
	int ijk = 0;

	///Read* current_read;

	bitmapper_bs_iter number_of_hits = 0;

	unsigned int bsSeq_max_length = 1000;
	char* bsSeq = (char*)malloc(sizeof(char)* bsSeq_max_length);

	bitmapper_bs_iter top, bot, match_length, total_match_length, match_step;
	bitmapper_bs_iter pre_top, pre_bot;


	bitmapper_bs_iter total_seed_length = 0;
	bitmapper_bs_iter total_seed_matches = 0;
	bitmapper_bs_iter total_seed_number = 0;
	bitmapper_bs_iter total_best_mapping_site = 0;

	bitmapper_bs_iter available_seed_length = 30;
	bitmapper_bs_iter max_seed_matches = 1000;



	available_seed_length = over_all_seed_length;


	bitmapper_bs_iter current_seed_length;

	bitmapper_bs_iter* locates;
	bitmapper_bs_iter* candidates;
	bitmapper_bs_iter candidate_length;
	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter locate_number_of_locations = 0;
	seed_votes* candidates_votes;
	bitmapper_bs_iter candidates_votes_length;

	bitmapper_bs_iter total_candidate_length = 0;
	bitmapper_bs_iter total_candidates_votes_length = 0;


	bitmapper_bs_iter max_seed_number;
	bitmapper_bs_iter seed_id;
	bitmapper_bs_iter empty_read = 0;
	bitmapper_bs_iter matched_read = 0;
	bitmapper_bs_iter unique_matched_read = 0;
	bitmapper_bs_iter unmatched_read = 0;
	bitmapper_bs_iter inner_i;

	bitmapper_bs_iter total_all_candidates_length = 0;
	unsigned int min_err;
	int min_err_index;

	/**
	bitmapper_bs_iter error_threshold = 31;
	if (thread_e != 255)
	{
		error_threshold = thread_e;
	}
	**/

	bitmapper_bs_iter max_error_threshold = 31;
	bitmapper_bs_iter error_threshold1;
	bitmapper_bs_iter error_threshold2;
	bitmapper_bs_iter large_error_threshold;



	bitmapper_bs_iter re_seed_length = 25;
	bitmapper_bs_iter re_seed_max_seed_matches = 2000;
	bitmapper_bs_iter re_seed_max_seed_number = 10;

	bitmapper_bs_iter min_candidates_votes_length;




	match_step = 8;
	max_seed_number = 25;
	int max_max_seed_number = max_seed_number;


	candidates = (bitmapper_bs_iter *)malloc(max_seed_number * max_seed_matches * sizeof(bitmapper_bs_iter));
	candidates_votes = (seed_votes *)malloc(max_seed_number * max_seed_matches * sizeof(seed_votes));

	bitmapper_bs_iter* re_seed_candidates
		= (bitmapper_bs_iter *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(bitmapper_bs_iter));

	seed_votes* re_seed_candidates_votes
		= (seed_votes *)malloc(re_seed_max_seed_number * re_seed_max_seed_matches * sizeof(seed_votes));


	uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(max_error_threshold * 2 + 1)*(SEQ_MAX_LENGTH));

	char* cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	char* best_cigar = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);

	char* path = (char*)malloc(sizeof(char)* SEQ_MAX_LENGTH);
	int* local_score = (int*)malloc(sizeof(int)* SEQ_MAX_LENGTH);

	Word* matrix_bit = (Word*)malloc(sizeof(Word)* SEQ_MAX_LENGTH * 8);

	int C_site_array[SEQ_MAX_LENGTH];
	int C_site_number;

	int start_site;

	double check_error_threshold = 0.1;
	bitmapper_bs_iter check_error;
	double second_best_check = 1.3;
	double score_threshold = 0.3;
	double edit_distance_threshold = 0.1;


	if (bs_score_threshold != -1)
	{
		score_threshold = bs_score_threshold;
	}

	if (bs_edit_distance_threshold != -1)
	{
		edit_distance_threshold = edit_distance_threshold;
	}

	if (bs_available_seed_length != -1)
	{
		available_seed_length = bs_available_seed_length;
	}

	double re_seed_threshold = 0.80;
	int reseed_length;
	bitmapper_bs_iter re_seed_error_threshold;

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	long long total_number_of_hit = 0;


	int is_mutiple_map = 0;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	long long direct_cut = 0;
	long long direct_match = 0;
	int C_site;
	///std::stringstream result_string;
	///std::string tmp_result_string;
	///char* tmp_result_string_char_array = NULL;
	///long long tmp_result_string_char_array_length;

	///batch_read_size

	///Read* read_batch = (Read*)malloc(sizeof(Read)*batch_read_size);




	Read* read_batch;

	Read_buffer_single_sub_block curr_sub_block;

	init_single_sub_block_single(&curr_sub_block);

	Output_buffer_sub_block current_sub_buffer;
	init_buffer_sub_block(&current_sub_buffer);





	///int inner_i;

	double startTime = Get_T();

	int file_flag = 1;
	int obtain_reads_num;


	int get_error;
	long long debug_1_mismatch = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;




	file_flag = 1;
	while (file_flag != 0)
	{



		/**
		pthread_mutex_lock(&readinputMutex);


		load_batch_read_pbat(read_batch, batch_read_size, &file_flag, &obtain_reads_num);

		pthread_mutex_unlock(&readinputMutex);
		**/

		file_flag = get_single_reads_mul_thread_pbat(&curr_sub_block);
		obtain_reads_num = curr_sub_block.sub_block_read_number;
		read_batch = curr_sub_block.read;



		enq_i = enq_i + obtain_reads_num;

		i = 0;
		while (i < obtain_reads_num)
		{

			///因为参考组是反的, 那么查询的时候实际上read也要反着查
			///所以实际上对是把read的反向做bs转换
			///C_to_T(current_read.rseq, bsSeq, current_read.length, &C_site);
			C_to_T_forward(read_batch[i].seq, bsSeq, read_batch[i].length, &C_site);



			error_threshold1 = thread_e_f * read_batch[i].length;
			if (error_threshold1 >= max_error_threshold)
			{
				error_threshold1 = max_error_threshold;
			}


			////fprintf(stderr, "i: %llu, thread_id: %d, obtain_reads_num: %d\n", thread_id, obtain_reads_num);


			total_match_length = 0;

			seed_id = 0;


			locates = candidates;

			candidate_length = 0;

			is_mutiple_map = 0;

			///这里要改
			get_error = -1;
			///这里要改
			extra_seed_flag = 1;



			/*******控制种子数量************************/
			max_seed_number = read_batch[i].length / 10 - 1;
			if (max_seed_number > max_max_seed_number)
			{
				max_seed_number = max_max_seed_number;
			}
			/*******控制种子数量************************/


			////5'端第一个seed单独处理
			if (seed_id < max_seed_number && total_match_length < read_batch[i].length)
			{
				current_seed_length = read_batch[i].length - total_match_length;

				number_of_hits =
					count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);

				///这里要改
				first_seed_match_length = match_length;






				/**
				if (number_of_hits == 1 && try_process_unique_exact_match_end_to_end_pbat_get_result_string(
					read_batch[i].length, read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual, cigar, read_batch[i].name,
					locates,
					top, tmp_ref, &match_length, read_batch[i].length - C_site - 1, result_string))**/
				if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end_pbat_get_output_buffer(
					read_batch[i].length, read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual, cigar, read_batch[i].name,
					locates,
					top, tmp_ref, &match_length, read_batch[i].length - C_site - 1, &current_sub_buffer, &get_error))
				{
					direct_match++;
					unique_matched_read++;
					matched_read++;
					i++;
					continue;


				}





				///只有第一个seed才有可能出现这种情况
				if (match_length == read_batch[i].length
					&&
					number_of_hits > 1)
				{
					is_mutiple_map = 1;

					///等于-1说明read里面没有C
					if (C_site == -1)
					{
						matched_read++;
						i++;
						direct_cut++;
						continue;
					}


				}




				///如果当前seed的只有一个匹配位置
				///则不管多长都是有效的
				///注意，如果number_of_hits == 1，则这个位置已经被拿到过了
				if (number_of_hits == 1)
				{
					total_seed_matches = total_seed_matches + number_of_hits;
					total_seed_number++;

					tmp_SA_length = 1;


					locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

					candidate_length = candidate_length + tmp_SA_length;

					locates = candidates + candidate_length;


				}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
				else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
				{

					total_seed_matches = total_seed_matches + number_of_hits;
					total_seed_number++;

					if (number_of_hits != 0)
					{

						tmp_SA_length = 0;

						locate_muti_thread(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length, &get_queue);
						reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);

						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;
						locates = candidates + candidate_length;

					}





				}



				if (match_length == 0)
				{
					total_match_length = total_match_length + match_step;
				}
				else
				{
					total_match_length = total_match_length + match_length / 2;
				}






				seed_id++;

			}











			///这里要改
			///说明是目前的最优解是1-mismatch
			if (get_error == 1)
			{

				second_seed_length = read_batch[i].length - first_seed_match_length;





				if (second_seed_length >= 17)
				{

					number_of_hits = count_hash_table(bsSeq, second_seed_length, &top, &bot, &pre_top, &pre_bot);

					if (number_of_hits == 1)
					{
						tmp_SA_length = 0;

						locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length,
							second_seed_length, first_seed_match_length);

						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;

						locates = candidates + candidate_length;

						extra_seed_flag = 0;


					}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
					else if (number_of_hits <= max_seed_matches)
					{

						if (number_of_hits != 0)
						{

							tmp_SA_length = 0;

							locate_muti_thread(bsSeq, top, bot, pre_top, pre_bot, locates, second_seed_length, &tmp_SA_length, &get_queue);

							reverse_and_adjust_site(locates, &tmp_SA_length, second_seed_length, first_seed_match_length);


							locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

							candidate_length = candidate_length + tmp_SA_length;
							locates = candidates + candidate_length;

						}

						extra_seed_flag = 0;

					}
					else if (number_of_hits > max_seed_matches)
					{
						extra_seed_flag = 1;
					}


				}
				else
				{
					extra_seed_flag = 1;
				}


			}

























			///这里要改
			if (extra_seed_flag == 1)
			{


				///这是剩下的seed
				while (seed_id < max_seed_number && total_match_length < read_batch[i].length)
				{



					current_seed_length = read_batch[i].length - total_match_length;

					number_of_hits =
						count_backward_as_much_1_terminate(bsSeq, current_seed_length, &top, &bot, &pre_top, &pre_bot, &match_length);


					///如果当前seed的只有一个匹配位置
					///则不管多长都是有效的
					if (number_of_hits == 1)
					{
						total_seed_matches = total_seed_matches + number_of_hits;
						total_seed_number++;

						tmp_SA_length = 0;
						locate_one_position_direct(locates, top, &tmp_SA_length, total_SA_length, match_length, total_match_length);

						locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

						candidate_length = candidate_length + tmp_SA_length;

						locates = candidates + candidate_length;


					}///如果当前有多个匹配位置，那么我们只把长度超过30的当做有效seed
					else if (match_length >= available_seed_length && number_of_hits <= max_seed_matches)
					{

						total_seed_matches = total_seed_matches + number_of_hits;
						total_seed_number++;

						if (number_of_hits != 0)
						{

							tmp_SA_length = 0;

							locate_muti_thread(bsSeq + current_seed_length - match_length, top, bot, pre_top, pre_bot, locates, match_length, &tmp_SA_length, &get_queue);
							reverse_and_adjust_site(locates, &tmp_SA_length, match_length, total_match_length);
							locate_number_of_locations = locate_number_of_locations + tmp_SA_length;

							candidate_length = candidate_length + tmp_SA_length;
							locates = candidates + candidate_length;

						}





					}
					else if (current_seed_length == match_length)
					{
						break;
					}



					if (match_length == 0)
					{
						total_match_length = total_match_length + match_step;
					}
					else
					{
						total_match_length = total_match_length + match_length / 2;
					}

					seed_id++;

				}
			}








			///这里要改
			///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
			if (extra_seed_flag == 0
				&&
				(candidate_length == 1 ||
				(candidate_length == 2 && candidates[0] == candidates[1])
				))
			{
				sprintf(cigar, "%dM", read_batch[i].length);

				output_sam_end_to_end_pbat_output_buffer
					(
					read_batch[i].name,
					read_batch[i].seq,
					read_batch[i].rseq,
					read_batch[i].qual,
					candidates[0],
					read_batch[i].length - 1,
					0,
					1,
					cigar,
					read_batch[i].length,
					&current_sub_buffer);

				unique_matched_read++;
				matched_read++;

				debug_1_mismatch++;
			}
			else if (candidate_length != 0)
			{

				std::sort(candidates, candidates + candidate_length);
				////generate_candidate_votes_allow_error(candidates, candidate_length, candidates_votes, &candidates_votes_length);

				///generate_candidate_votes(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold);
				generate_candidate_votes_shift(candidates, candidate_length, candidates_votes, &candidates_votes_length, error_threshold1);

				total_candidate_length = total_candidate_length + candidate_length;
				total_candidates_votes_length = total_candidates_votes_length + candidates_votes_length;



				///这个排序其实有助于剪枝
				///其实好像也没什么用..
				std::sort(candidates_votes, candidates_votes + candidates_votes_length, compare_seed_votes);


				///感觉我这样做的话，不可能出现有单独一个精确匹配的位置的情况了
				///如果一个read有唯一的精确匹配位置，那么他肯定在前面就被报告出来了
				if (is_mutiple_map == 0)
				{


					if (error_threshold1 <= 15)
					{
						map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);


					}
					else
					{
						map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);
					}


				}
				else
				{
					total_is_mutiple_map++;


					if (error_threshold1 <= 15)
					{
						map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);

					}
					else
					{
						map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index);
					}


				}


				min_candidates_votes_length = 0;

				if (min_err_index >= 0)
				{

					total_best_mapping_site++;
					matched_read++;
					min_candidates_votes_length = 0;

					calculate_best_map_cigar_end_to_end_pbat_output_buffer
						(candidates_votes, &min_candidates_votes_length,
						read_batch[i].length, error_threshold1, tmp_ref,
						read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual,
						cigar, path, matrix, matrix_bit, read_batch[i].name,
						best_cigar, &current_sub_buffer);

					unique_matched_read++;


				}
				else if (min_err_index == -1)
				{
					unmatched_read++;
				}
				else
				{
					matched_read++;
				}







			}
			else
			{
				empty_read++;
			}

			i++;

		}

	/**
		pthread_mutex_lock(&queueMutex);

		fprintf(output_file, "%s", result_string.str().c_str());
		///out << result_string.str();

		pthread_mutex_unlock(&queueMutex);

		///stringstream和string都要清空
		result_string.clear();
		result_string.str("");


	}
	**/
		current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

		push_results_to_buffer(&current_sub_buffer);

		current_sub_buffer.length = 0;

	}

	finish_output_buffer();



	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;

	fprintf(stdout, "thread %d completed!\n", thread_id);

	///这里要改
	fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);
}






