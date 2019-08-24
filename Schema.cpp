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
#include "ksw.h"

#include "Ref_Genome.h"
#include "bam_prase.h"

Inpute_PE_methy_alignment_buffer methy_input_buffer;
PE_methy_parameters parameters;
char methy_hash[5] = { 'A', 'C', 'G', 'T', 'N' };

FILE* input_methy_alignment_file;

char char_to_3bit[128] = {
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4
};


int8_t mat[25] = {
	0, -1, -1, -1, -1,
	-1, 0, -1, -1, -1,
	-1, -1, 0, -1, -1,
	-1, 0, -1, 0, -1,
	-1, -1, -1, -1, -1
};

int8_t mat_diff[25] = {
	0, -1, -1, -1, -1,
	-1, 0, -1, -1, -1,
	-1, -1, 0, -1, -1,
	-1, 0, -1, 0, -1,
	-1, -1, -1, -1, -1
};





_rg_name_l  *_ih_refGenName;
int refChromeCont;

char *versionN = "1.0.2.1";
long long mappingCnt[MAX_Thread];
unsigned int done;
long long mappedSeqCnt[MAX_Thread];
long long completedSeqCnt[MAX_Thread];
///pair_distance_count PE_distance[MAX_Thread];
pair_distance_count* PE_distance;
///char *Mapped_File;
char *_msf_refGen = NULL;
unsigned char *_BS_refGen = NULL;
bitmapper_bs_iter _msf_refGenLength = 0;
bitmapper_bs_iter _msf_refGenLength_2_bit = 0;
Commpress_HashTable  chhy_HashTable;
bitmapper_bs_iter total_SA_length;

long long unique_mapped_read[MAX_Thread];
long long ambiguous_mapped_read[MAX_Thread];

long long mapped_bases[MAX_Thread];
long long error_mapped_bases[MAX_Thread];
///long long unmapped_read[MAX_Thread];

int batch_read_size = 50000;


OPT_FIELDS **_chhy_optionalFields;

FILE *schema_out_fp;
FILE *schema_input_fp;

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

int read_format = 1;


Methylation_Output methy_out;



///bitmapper_bs_iter debug_pos;


inline int MAP_Calculation(unsigned int second_best_diff, unsigned int error_threshold,
int best_score,int GapOpenPenalty, int GapExtensionPenalty, int MistMatchPenaltyMax)
{

	int scoreMax = GapOpenPenalty + GapExtensionPenalty;
	if(scoreMax < MistMatchPenaltyMax)
	{
		scoreMax = MistMatchPenaltyMax;
	}
	///scoreMax < 0, is a negative number
	scoreMax = -scoreMax * error_threshold;

	///the max score and min score are 0 and scoreMax, respectively
	///and scoreMax is a negative number
	///thus scoreMaxRange = -scoreMax
	int scoreMaxRange = -scoreMax;
	///both best_score and scoreMax are negative numbers
	int score_diff = best_score - scoreMax;
	if(score_diff < 0)
	{
		fprintf(stderr, "error best_score: %d, scoreMax: %d\n", best_score, scoreMax);
		score_diff = 0;
	}
	int error_diff = second_best_diff;
	///only one alignment
	///in this case, second_best_diff is very large
	if(second_best_diff > error_threshold)
	{
		error_diff = error_threshold + 1;
	}

	double rank;
	double rank_error;

	///only one alignment
	if(error_diff > error_threshold)
	{
		rank = (double)((double)score_diff/(double)scoreMaxRange);
		if(rank >= 0.8)
		{
			return 42;
		}
		else if(rank >= 0.7)
		{
			return 40;
		}
		else if(rank >= 0.6)
		{
			return 24;
		}
		else if(rank >= 0.5)
		{
			return 23;
		}
		else if(rank >= 0.4)
		{
			return 8;
		}
		else if(rank >= 0.3)
		{
			return 3;
		}
		else
		{
			return 0;
		}
		
	}
	else
	{
		rank_error = (double)((double)error_diff / (double)error_threshold);
		rank = (double)((double)score_diff/(double)scoreMaxRange);

		if(rank_error >= 0.9)
		{
			if(best_score == 0)
			{
				return 39;
			}
			else
			{
				return 33;
			}
		}
		else if(rank_error >= 0.8)
		{
			if(best_score == 0)
			{
				return 38;
			}
			else
			{
				return 27;
			}
		}
		else if(rank_error >= 0.7)
		{
			if(best_score == 0)
			{
				return 37;
			}
			else
			{
				return 26;
			}
		}
		else if(rank_error >= 0.6)
		{
			if(best_score == 0)
			{
				return 36;
			}
			else
			{
				return 22;
			}
		}
		else if(rank_error >= 0.5)
		{
			if(best_score == 0)
			{
				return 35;
			}
			else if(rank >= 0.84)
			{
				return 25;
			}
			else if(rank >= 0.68)
			{
				return 16;
			}
			else
			{
				return 5;
			}
		}
		else if(rank_error >= 0.4)
		{
			if(best_score == 0)
			{
				return 34;
			}
			else if(rank >= 0.84)
			{
				return 21;
			}
			else if(rank >= 0.68)
			{
				return 14;
			}
			else
			{
				return 4;
			}
		}
		else if(rank_error >= 0.3)
		{
			if(best_score == 0)
			{
				return 32;
			}
			else if(rank >= 0.88)
			{
				return 18;
			}
			else if(rank >= 0.67)
			{
				return 15;
			}
			else
			{
				return 3;
			}
		}
		else if(rank_error >= 0.2)
		{
			if(best_score == 0)
			{
				return 31;
			}
			else if(rank >= 0.88)
			{
				return 17;
			}
			else if(rank >= 0.67)
			{
				return 11;
			}
			else
			{
				return 0;
			}
		}
		else if(rank_error >= 0.1)
		{
			if(best_score == 0)
			{
				return 30;
			}
			else if(rank >= 0.88)
			{
				return 12;
			}
			else if(rank >= 0.67)
			{
				return 7;
			}
			else
			{
				return 0;
			}
		}
		else if(error_diff == 0)
		{
			if(rank >= 0.67)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
		else
		{
			if(rank >= 0.67)
			{
				return 6;
			}
			else
			{
				return 2;
			}
		}
	}


}




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
long long* number_of_ambiguous_mapped_read, long long* number_of_unmapped_read, 
long long* number_of_mapped_bases,
long long* number_of_mapped_errors)
{
	*number_of_read = 0;
	*number_of_unique_mapped_read = 0;
	*number_of_ambiguous_mapped_read = 0;
	*number_of_mapped_bases = 0;
	*number_of_mapped_errors = 0;

	int i = 0;
	for (i = 0; i<THREAD_COUNT; i++)
	{
		*number_of_read += completedSeqCnt[i];
		*number_of_unique_mapped_read += unique_mapped_read[i];
		*number_of_ambiguous_mapped_read += ambiguous_mapped_read[i];
		///*number_of_unmapped_read += unmapped_read[i];
		*number_of_mapped_bases += mapped_bases[i];
		*number_of_mapped_errors += error_mapped_bases[i];

	}

	*number_of_unmapped_read = *number_of_read - *number_of_unique_mapped_read - *number_of_ambiguous_mapped_read;
}


bool compare_pe_distance(const pe_distance_info& s1, const pe_distance_info& s2)
{
	return s1.count > s2.count; //从大到小排序
}

void out_paired_distance_statistic()
{
	int i = 0;
	///这里一定一定要+1
	int length = maxDistance_pair + 1;
	
	
	fprintf(stderr, "@%llu\t%llu\n", maxDistance_pair, minDistance_pair);


	pair_distance_result sort_results;
	sort_results.count = (pe_distance_info*)malloc(sizeof(pe_distance_info)*length);

	
	int j;
	
	for (j = 0; j < genome_cuts; j++)
	{
		int k = 0;

		for (k = 0; k < length; k++)
		{
			long long tmp = 0;
			for (i = 0; i < THREAD_COUNT; i++)
			{
				tmp = tmp + PE_distance[i].count[j][k];
			}


			sort_results.count[k].count = tmp;
			sort_results.count[k].index = k;
			
		
		}


		std::sort(sort_results.count, sort_results.count + length, compare_pe_distance);


		double total_frequence = 0;

		for (k = 0; k < length; k++)
		{
			total_frequence += sort_results.count[k].count;
		}

		double current_frequence = 0;


		for (k = 0; k < length; k++)
		{
			current_frequence += sort_results.count[k].count;

			double prec = current_frequence / total_frequence;
			prec = prec * 100;



			if ((long long)prec % 5 == 0)
			{
				fprintf(stderr, "%.2f, %lld, %llu\t", prec, (long long)current_frequence, k);
			}

		}

		fprintf(stderr, "\n\n");
	}
	

	
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

void open_methy_output(Methylation_Output* methy_out, char* file_name, bitmapper_bs_iter cut_length)
{
	(*methy_out).genome_cut = genome_cuts;
	(*methy_out).files = (FILE**)malloc(sizeof(FILE*)*(*methy_out).genome_cut);
	(*methy_out).cut_index = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*((*methy_out).genome_cut + 1));
	int i = 0;
	char tmp_file_name[SEQ_MAX_LENGTH];

	for (i = 0; i < (*methy_out).genome_cut; i++)
	{
		sprintf(tmp_file_name, "%s_%d.bmm", file_name, i);
		(*methy_out).files[i] = fopen(tmp_file_name, "w");
		(*methy_out).cut_index[i] = 0;

		fprintf((*methy_out).files[i], "%llu\t%d\t%d\t", (*methy_out).genome_cut, maxDistance_pair, is_pairedEnd);

		fprintf((*methy_out).files[i], "%llu\t", cut_length*i);

		if (cut_length*(i + 1) - 1 > _msf_refGenLength - 1)
		{
			fprintf((*methy_out).files[i], "%llu\n", _msf_refGenLength - 1);

			if (i != (*methy_out).genome_cut - 1)
			{
				fprintf(stderr, "ERROR!\n");
				exit(0);

			}
		}
		else
		{
			fprintf((*methy_out).files[i], "%llu\n", cut_length*(i + 1) - 1);
		}
	}
	(*methy_out).cut_index[(*methy_out).genome_cut] = 0;

}

void Prepare_alignment(char* outputFileName, char *genFileName, _rg_name_l *chhy_ih_refGenName, int chhy_refChromeCont, int i_read_format, int is_pairedEnd)
{
	read_format = i_read_format;

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

	if (output_methy == 0)
	{
		schema_out_fp = get_Ouput_Dec();
	}

	map_cigar[0] = 'M';
	map_cigar[1] = 'M';
	map_cigar[2] = 'D';
	map_cigar[3] = 'I';

	map_cigar[5] = 'M';
	map_cigar[6] = 'M';
	map_cigar[7] = 'D';
	map_cigar[8] = 'I';

	if (output_methy == 1)
	{
		if (is_pairedEnd)
		{

			long long total_bit = maxDistance_pair*genome_cuts;

			genome_cuts = total_bit / 500;

			if (total_bit % 500 != 0)
			{
				genome_cuts++;
			}

			
		}

		///genome_cuts = 64;
		fprintf(stderr, "genome_cuts: %d\n", genome_cuts);

		///cut_length = (_msf_refGenLength * 2) / genome_cuts;
		cut_length = _msf_refGenLength / genome_cuts;
		////if ((_msf_refGenLength * 2) % genome_cuts != 0)
		if (_msf_refGenLength % genome_cuts != 0)
		{
			cut_length = cut_length + 1;
		}

		///outputFileName一定是bmm_folder/output
		open_methy_output(&methy_out, outputFileName, cut_length);


		if (is_pairedEnd)
		{
			PE_distance = (pair_distance_count*)malloc(sizeof(pair_distance_count)*THREAD_COUNT);
		}

		
	}

	///exit(0);

	/**
	int8_t mat[25] = {
	0, -1, -1, -1, -1,
	-1, 0, -1, -1, -1,
	-1, -1, 0, -1, -1,
	-1, 0, -1, 0, -1,
	-1, -1, -1, -1, -1
	};
	**/
	
	int j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
		{
			mat_diff[k] = i == j? 0 : (MistMatchPenaltyMax - MistMatchPenaltyMin);
			mat[k++] = i == j? 0 : -MistMatchPenaltyMin;
		}
		
		mat_diff[k] = 0;
		mat[k++] = -N_Penalty; // ambiguous base
	}

	
	for (j = 0; j < 5; ++j) 
	{
		mat_diff[k] = 0;
		mat[k++] = -N_Penalty;
	}

	mat_diff[16] = 0;
	mat[16] = 0;
	
	/**
	fprintf(stderr, "mat: ");
	for (i = 0; i < 25; i++)
	{
		fprintf(stderr, "%d, ", mat[i]);
	}


	fprintf(stderr, "\nmat_diff: ");
	for (i = 0; i < 25; i++)
	{
		fprintf(stderr, "%d, ", mat_diff[i]);
	}
	**/
	
	

}

void get_genome_cuts(char* file_name)
{
	char tmp_file_name[SEQ_MAX_LENGTH];

	sprintf(tmp_file_name, "%s_%d.bmm", file_name, 0);
	FILE* read_file = fopen(tmp_file_name, "r");

	bitmapper_bs_iter start_pos, end_pos;

	fscanf(read_file, "%llu\t%d\t%d\t%llu\t%llu\n", &genome_cuts, &maxDistance_pair, &is_pairedEnd, &start_pos, &end_pos);

	fclose(read_file);
}

void Prepare_methy(char *genFileName, _rg_name_l *chhy_ih_refGenName, int chhy_refChromeCont)
{
	_ih_refGenName = chhy_ih_refGenName;
	refChromeCont = chhy_refChromeCont;
	_msf_refGenLength = getRefGenomeLength();

	schema_input_fp = get_index_file();

}

int get_single_result(int* flag, bitmapper_bs_iter pos, char* buffer)
{
	

}






int init_genome_cut_PE(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	char* ref_genome, bitmapper_bs_iter length,
	bitmapper_bs_iter index, bitmapper_bs_iter extra_length,
	char* last_two_bases,
	deduplicate_PE* PE_de)
{
	if (index == 0)
	{
		memset(methy, 0, sizeof(uint16_t)*(length + extra_length));
		memset(total, 0, sizeof(uint16_t)*(length + extra_length));
		memset(correct_methy, 0, sizeof(uint16_t)*(length + extra_length));
		memset(correct_total, 0, sizeof(uint16_t)*(length + extra_length));


		///如果第一个block也是最后一个
		///那么读参考基因组的时候不应该+extra_length
		if (index != genome_cuts - 1)
		{
			length = length + extra_length;
		}


	}
	else if (index == genome_cuts - 1)
	{
		memset(methy, 0, sizeof(uint16_t)*length);
		memset(total, 0, sizeof(uint16_t)*length);
		memset(correct_methy, 0, sizeof(uint16_t)*length);
		memset(correct_total, 0, sizeof(uint16_t)*length);


		memcpy(methy, methy + length, extra_length*sizeof(uint16_t));
		memcpy(total, total + length, extra_length*sizeof(uint16_t));
		memcpy(correct_methy, correct_methy + length, extra_length*sizeof(uint16_t));
		memcpy(correct_total, correct_total + length, extra_length*sizeof(uint16_t));

		memset(methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(total + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_total + length, 0, sizeof(uint16_t)*extra_length);

		///注意对最后一个元素,这东西不能要
		///length = length + extra_length;

		fseek(schema_input_fp, extra_length*-1, SEEK_CUR);
	}
	else
	{

		memset(methy, 0, sizeof(uint16_t)*length);
		memset(total, 0, sizeof(uint16_t)*length);
		memset(correct_methy, 0, sizeof(uint16_t)*length);
		memset(correct_total, 0, sizeof(uint16_t)*length);


		memcpy(methy, methy + length, extra_length*sizeof(uint16_t));
		memcpy(total, total + length, extra_length*sizeof(uint16_t));
		memcpy(correct_methy, correct_methy + length, extra_length*sizeof(uint16_t));
		memcpy(correct_total, correct_total + length, extra_length*sizeof(uint16_t));

		memset(methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(total + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_total + length, 0, sizeof(uint16_t)*extra_length);

		length = length + extra_length;

		fseek(schema_input_fp, extra_length*-1, SEEK_CUR);
	}

	re_init_deduplicate_PE(PE_de);


	if (fread(ref_genome, 1, length, schema_input_fp) != length)
	{
		fprintf(stderr, "index: %llu\n", index);

		return 0;
	}
	else
	{
		return 1;
	}

}















int init_genome_cut(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	char* ref_genome, bitmapper_bs_iter length,
	bitmapper_bs_iter index, bitmapper_bs_iter extra_length,
	char* last_two_bases)
{
	if (index == 0)
	{
		memset(methy, 0, sizeof(uint16_t)*(length + extra_length));
		memset(total, 0, sizeof(uint16_t)*(length + extra_length));
		memset(correct_methy, 0, sizeof(uint16_t)*(length + extra_length));
		memset(correct_total, 0, sizeof(uint16_t)*(length + extra_length));


		///如果第一个block也是最后一个
		///那么读参考基因组的时候不应该+extra_length
		if (index != genome_cuts - 1)
		{
			length = length + extra_length;
		}

		
	}
	else if (index == genome_cuts - 1)
	{
		memset(methy, 0, sizeof(uint16_t)*length);
		memset(total, 0, sizeof(uint16_t)*length);
		memset(correct_methy, 0, sizeof(uint16_t)*length);
		memset(correct_total, 0, sizeof(uint16_t)*length);


		memcpy(methy, methy + length, extra_length*sizeof(uint16_t));
		memcpy(total, total + length, extra_length*sizeof(uint16_t));
		memcpy(correct_methy, correct_methy + length, extra_length*sizeof(uint16_t));
		memcpy(correct_total, correct_total + length, extra_length*sizeof(uint16_t));

		memset(methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(total + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_total + length, 0, sizeof(uint16_t)*extra_length);

		///注意对最后一个元素,这东西不能要
		///length = length + extra_length;

		fseek(schema_input_fp, extra_length*-1, SEEK_CUR);
	}
	else
	{
		
		memset(methy, 0, sizeof(uint16_t)*length);
		memset(total, 0, sizeof(uint16_t)*length);
		memset(correct_methy, 0, sizeof(uint16_t)*length);
		memset(correct_total, 0, sizeof(uint16_t)*length);
		

		memcpy(methy, methy + length, extra_length*sizeof(uint16_t));
		memcpy(total, total + length, extra_length*sizeof(uint16_t));
		memcpy(correct_methy, correct_methy + length, extra_length*sizeof(uint16_t));
		memcpy(correct_total, correct_total + length, extra_length*sizeof(uint16_t));

		memset(methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(total + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_methy + length, 0, sizeof(uint16_t)*extra_length);
		memset(correct_total + length, 0, sizeof(uint16_t)*extra_length);

		length = length + extra_length;

		fseek(schema_input_fp, extra_length*-1, SEEK_CUR);
	}

	
	if (fread(ref_genome, 1, length, schema_input_fp) != length)
	{
		fprintf(stderr, "index: %llu\n",index);

		return 0;
	}
	else
	{
		return 1;
	}
		
}
//C_or_G， 1代表C, 2代表G
//direction, 1代表forward, 2代表reverse comp
inline void get_direction(uint16_t flag, int* C_or_G, int* direction)
{
	if (flag&IS_PAIRED) ///双端
	{
		if ((flag & 0x50) == 0x50)//Read1, reverse comp. == OB, flag是83
		{
			(*C_or_G) = 2;
			(*direction) = 2;
		}
		else if (flag & 0x40)//Read1, forward == OT, flag是99
		{
			(*C_or_G) = 1;
			(*direction) = 1;
		}
		else if ((flag & 0x90) == 0x90)//Read2, reverse comp. == OT, flag是147
		{
			(*C_or_G) = 1;
			(*direction) = 2;
		}
		else if (flag & 0x80)//Read2, forward == OB, flag是163
		{
			(*C_or_G) = 2;
			(*direction) = 1;
		}
		else
		{
			(*C_or_G) = 0;
			(*direction) = 0;
		}
	}
	else   ///单端
	{
		if (flag & 0x10)//Reverse comp. == OB
		{
			(*C_or_G) = 2;
			(*direction) = 2;
		}
		else  //OT
		{
			(*C_or_G) = 1;
			(*direction) = 1;
		}
		
	}

}
///0就是重复了, 1就是没重复
inline int remove_dup(char* ref_genome, int direction, bitmapper_bs_iter pos, int seq_length)
{

	///return 1;


	bitmapper_bs_iter dup_pos;

	


	///说明这是rc链
	if (direction == 2)
	{
		dup_pos = pos + seq_length - 1;

		/**
		if (pos == 976903)
		{
			fprintf(stderr, "pos: %llu, dup_pos: %llu, direction: %d\n", pos, dup_pos, direction);
			fprintf(stderr, "ref_genome[dup_pos]: %u\n", ref_genome[dup_pos]);

		}

		if (dup_pos == 977002)
		{
			fprintf(stderr, "***************pos: %llu, seq_length: %llu\n", pos, seq_length);
		}
		**/


		///这就说明这个是重复
		if (ref_genome[dup_pos] & RC_DUP_MASK)
		{
			return 0;
		}

		ref_genome[dup_pos] = ref_genome[dup_pos] | RC_DUP_MASK;

		return 1;

	}
	else   ///正向链
	{
		dup_pos = pos;
		///这就说明这个是重复
		if (ref_genome[dup_pos] & FORWARD_DUP_MASK)
		{
			return 0;
		}

		ref_genome[dup_pos] = ref_genome[dup_pos] | FORWARD_DUP_MASK;

		return 1;
	}

	return 1;
}


///char methy_hash[5] = { 'A', 'C', 'G', 'T', 'N' };
inline int convert_3bit_to_char(bitmapper_bs_iter* tmp_buffer, int read_length, char* read)
{
	int i = 0;
	int word_21_index = 0;
	int word_21_inner_index = 0;
	bitmapper_bs_iter mask = 7;
	bitmapper_bs_iter tmp;
	for (i = 0; i < read_length; i++)
	{
		tmp = tmp_buffer[word_21_index] >> (word_21_inner_index * 3);
		tmp = tmp & 7;
		read[i] = methy_hash[tmp];

		word_21_inner_index++;

		if (word_21_inner_index == 21)
		{
			word_21_inner_index = 0;
			word_21_index++;
		}
	}

	read[read_length] = '\0';
}


inline int input_single_alignment(FILE* read_file,
	uint16_t* flag, bitmapper_bs_iter* pos, int* seq_length, bitmapper_bs_iter* read)
{


	fread(flag, sizeof((*flag)), 1, read_file);

	if (feof(read_file))
	{
		return 0;

	}

	fread(pos, sizeof((*pos)), 1, read_file);
	fread(seq_length, sizeof((*seq_length)), 1, read_file);

	/**********这里要改**********/
	///fread(read, 1, (*seq_length), read_file);
	///read[(*seq_length)] = '\0';


	fread(read, sizeof(bitmapper_bs_iter),
		((*seq_length) / 21 + ((*seq_length) % 21 != 0)), read_file);
	
	/**********这里要改**********/

	return 1;
}


///1是拿到一个有效的alignment
///0是读到文件末尾了
///-1是这个alignment重复了, 拿到了一个空的alignment
inline int get_alignment(FILE* read_file, char* ref_genome, bitmapper_bs_iter start_pos, bitmapper_bs_iter* read,
	int* g_seq_length, int* g_C_or_G, int* g_direction, bitmapper_bs_iter* g_pos)
{


	uint16_t flag;
	bitmapper_bs_iter pos;
	
	int C_or_G;
	int direction;
	int seq_length;

	

	///if (fscanf(read_file, "%llu\t%llu\t%s\n", &flag, &pos, read) != EOF)
	if(input_single_alignment(read_file,
		&flag, &pos, &seq_length, read))
	{
		


		get_direction(flag, &C_or_G, &direction);

		/**
		if (pos == 976934)
		{
			fprintf(stderr, "####################pos: %llu, direction: %d, seq_length: %llu\n", pos, direction, seq_length);
		}
		**/
		/**
		if (pos >= 976914 - seq_length + 1 && pos <= 976914 + 1)
		{
			fprintf(stderr, "+pos: %llu, direction: %d, inner_pos: %llu\n", pos, direction, pos - start_pos);
			char hahahh_seq[1000];
			convert_3bit_to_char(read, seq_length, hahahh_seq);

			fprintf(stderr, "%s\n", hahahh_seq);
		}
		**/
		
		


		///seq_length = strlen(read);
		///如果重复了
		if (!remove_dup(ref_genome, direction, pos - start_pos, seq_length))
		{
			return -1;
		}

		/**
		if (pos >= 976914 - seq_length + 1 && pos <= 976914 + 1)
		{
			fprintf(stderr, "-pos: %llu, direction: %d\n", pos, direction);
			
			char hahahh_seq[1000];
			convert_3bit_to_char(read, seq_length, hahahh_seq);

			fprintf(stderr, "%s\n", hahahh_seq);
			
		}
		**/


		
		(*g_C_or_G) = C_or_G;
		(*g_direction) = direction;
		(*g_seq_length) = seq_length;
		(*g_pos) = pos;
		return 1;
	}
	else
	{
		return 0;
	}

}



///0就是重复了, 1就是没重复，2是异常值
inline int remove_dup_PE(char* ref_genome, 
	int direction1, bitmapper_bs_iter pos1, int seq_length1,
	int direction2, bitmapper_bs_iter pos2, int seq_length2,
	deduplicate_PE* PE_de, int PE_distance)
{

	///return 1;


	bitmapper_bs_iter dup_pos1;
	bitmapper_bs_iter dup_pos2;

	bitmapper_bs_iter distance;

	///说明read1比对到rc链, read2比对到正向链
	if (direction1 == 2)
	{
		if (direction2 != 1)
		{
			fprintf(stderr, "ERROR!\n");
		}

		dup_pos1 = pos1 + seq_length1 - 1;
		dup_pos2 = pos2;


		///按道理来说, 这种情况dup_pos2必须小于等于dup_pos1, 但总有特殊情况
		///这个后面要单独处理
		if (dup_pos2 > dup_pos1)
		{
			return 2;
		}

		///这个时候distance绝对是>=0
		distance = dup_pos1 - dup_pos2;

		if (distance>PE_distance)
		{
			fprintf(stderr, "ERROR!\n");
		}
		

		///判断重复
		if ((ref_genome[dup_pos2] & RC_DUP_MASK)
			&&
			get_deduplicate_PE_mask(PE_de, PE_de->rc_bits, dup_pos2, distance))
		{

			return 0;
		}
		else
		{
			///没重复则要置位为重复
			ref_genome[dup_pos2] = ref_genome[dup_pos2] | RC_DUP_MASK;
			set_deduplicate_PE_mask(PE_de, PE_de->rc_bits, dup_pos2, distance);
		}
		

		return 1;

	}
	else   ///说明read1比对到正向链, read2比对到rc链
	{

		if (direction2 != 2)
		{
			fprintf(stderr, "ERROR!\n");
		}


		dup_pos1 = pos1;
		dup_pos2 = pos2 + seq_length2 - 1;


		///按道理来说, 这种情况dup_pos1必须小于等于dup_pos2, 但总有特殊情况
		///这个后面要单独处理
		if (dup_pos1 > dup_pos2)
		{
			return 2;
		}


		///这个时候distance绝对是>=0
		distance = dup_pos2 - dup_pos1;

		if (distance>PE_distance)
		{
			fprintf(stderr, "ERROR!\n");
		}


		///判断重复
		if ((ref_genome[dup_pos1] & FORWARD_DUP_MASK)
			&&
			get_deduplicate_PE_mask(PE_de, PE_de->forward_bits, dup_pos1, distance))
		{

			return 0;
		}
		else
		{
			///没重复则要置位为重复
			ref_genome[dup_pos1] = ref_genome[dup_pos1] | FORWARD_DUP_MASK;
			set_deduplicate_PE_mask(PE_de, PE_de->forward_bits, dup_pos1, distance);
		}

		return 1;
	}

	return 1;
}




///0就是重复了, 1就是没重复，2是异常值
inline int remove_abnormal_dup_PE(char* ref_genome,
	int direction1, bitmapper_bs_iter pos1, int seq_length1,
	int direction2, bitmapper_bs_iter pos2, int seq_length2,
	deduplicate_PE* PE_de, int PE_distance)
{

	///return 1;


	bitmapper_bs_iter dup_pos1;
	bitmapper_bs_iter dup_pos2;

	bitmapper_bs_iter distance;

	///说明read1比对到rc链, read2比对到正向链
	if (direction1 == 2)
	{
		if (direction2 != 1)
		{
			fprintf(stderr, "ERROR!\n");
		}

		dup_pos1 = pos1 + seq_length1 - 1;
		dup_pos2 = pos2;


		///按道理来说, 这种情况dup_pos2必须小于等于dup_pos1
		///但这个函数本身就是在处理异常情况了
		///所以dup_pos2必须大于dup_pos1
		if (dup_pos2 <= dup_pos1)
		{
			fprintf(stderr, "ERROR!\n");
		}

		///这个时候distance绝对是>=0
		distance = dup_pos2 - dup_pos1;

		if (distance>PE_distance)
		{
			fprintf(stderr, "ERROR!\n");
		}


		/**
		///大的那个是0是1无所谓，小的这个变成0就好了
		unset_deduplicate_PE_mask(PE_de, PE_de->rc_bits, dup_pos1, distance);
		**/


		///判断重复
		if ((ref_genome[dup_pos1] & RC_DUP_MASK)
			&&
			get_deduplicate_PE_mask(PE_de, PE_de->rc_bits, dup_pos1, distance))
		{

			return 0;
		}
		else
		{
			///没重复则要置位为重复
			ref_genome[dup_pos1] = ref_genome[dup_pos1] | RC_DUP_MASK;
			set_deduplicate_PE_mask(PE_de, PE_de->rc_bits, dup_pos1, distance);
		}









		return 1;

	}
	else   ///说明read1比对到正向链, read2比对到rc链
	{

		if (direction2 != 2)
		{
			fprintf(stderr, "ERROR!\n");
		}


		dup_pos1 = pos1;
		dup_pos2 = pos2 + seq_length2 - 1;


		///按道理来说, 这种情况dup_pos1必须小于等于dup_pos2
		///但因为我们本来就是处理异常情况
		///所以dup_pos1必须大于dup_pos2
		if (dup_pos1 <= dup_pos2)
		{
			fprintf(stderr, "ERROR!\n");
		}


		///这个时候distance绝对是>=0
		distance = dup_pos1 - dup_pos2;

		if (distance>PE_distance)
		{
			fprintf(stderr, "ERROR!\n");
		}


		///大的那个是0是1无所谓，小的这个变成0就好了
		///unset_deduplicate_PE_mask(PE_de, PE_de->forward_bits, dup_pos2, distance);


		///判断重复
		if ((ref_genome[dup_pos2] & FORWARD_DUP_MASK)
			&&
			get_deduplicate_PE_mask(PE_de, PE_de->forward_bits, dup_pos2, distance))
		{

			return 0;
		}
		else
		{
			///没重复则要置位为重复
			ref_genome[dup_pos2] = ref_genome[dup_pos2] | FORWARD_DUP_MASK;
			set_deduplicate_PE_mask(PE_de, PE_de->forward_bits, dup_pos2, distance);
		}


		return 1;
	}

	return 1;
}








///0就是重复了, 1就是没重复，2是异常值
inline int unset_abnormal_dup_PE(char* ref_genome,
	int direction1, bitmapper_bs_iter pos1, int seq_length1,
	int direction2, bitmapper_bs_iter pos2, int seq_length2,
	deduplicate_PE* PE_de, int PE_distance)
{

	///return 1;


	bitmapper_bs_iter dup_pos1;
	bitmapper_bs_iter dup_pos2;

	bitmapper_bs_iter distance;

	///说明read1比对到rc链, read2比对到正向链
	if (direction1 == 2)
	{
		if (direction2 != 1)
		{
			fprintf(stderr, "unset ERROR!\n");
		}

		dup_pos1 = pos1 + seq_length1 - 1;
		dup_pos2 = pos2;


		///按道理来说, 这种情况dup_pos2必须小于等于dup_pos1
		///但这个函数本身就是在处理异常情况了
		///所以dup_pos2必须大于dup_pos1
		if (dup_pos2 <= dup_pos1)
		{
			fprintf(stderr, "unset ERROR!\n");
		}

		///这个时候distance绝对是>=0
		distance = dup_pos2 - dup_pos1;

		if (distance>PE_distance)
		{
			fprintf(stderr, "unset ERROR!\n");
		}

		///大的那个是0是1无所谓，小的这个变成0就好了
		unset_deduplicate_PE_mask(PE_de, PE_de->rc_bits, dup_pos1, distance);

		/**
		if (get_deduplicate_PE_mask(PE_de, PE_de->rc_bits, dup_pos1, distance) == 1)
		{
			fprintf(stderr, "SBSBSBSBS!\n");
		}
		**/

		return 1;

	}
	else   ///说明read1比对到正向链, read2比对到rc链
	{

		if (direction2 != 2)
		{
			fprintf(stderr, "unset ERROR!\n");
		}


		dup_pos1 = pos1;
		dup_pos2 = pos2 + seq_length2 - 1;


		///按道理来说, 这种情况dup_pos1必须小于等于dup_pos2
		///但因为我们本来就是处理异常情况
		///所以dup_pos1必须大于dup_pos2
		if (dup_pos1 <= dup_pos2)
		{
			fprintf(stderr, "unset ERROR!\n");
		}


		///这个时候distance绝对是>=0
		distance = dup_pos1 - dup_pos2;

		if (distance>PE_distance)
		{
			fprintf(stderr, "unset ERROR!\n");
		}


		///大的那个是0是1无所谓，小的这个变成0就好了
		unset_deduplicate_PE_mask(PE_de, PE_de->forward_bits, dup_pos2, distance);

		/**
		if (get_deduplicate_PE_mask(PE_de, PE_de->forward_bits, dup_pos2, distance) == 1)
		{
			fprintf(stderr, "SBSBSBSBS!\n");
		}
		**/

		return 1;
	}

	return 1;
}


inline int input_PE_alignment(FILE* read_file, 
	uint16_t* flag1, bitmapper_bs_iter* pos1, int* seq_length1, bitmapper_bs_iter* read1,
	uint16_t* flag2, bitmapper_bs_iter* pos2, int* seq_length2, bitmapper_bs_iter* read2)
{

		fread(flag1, sizeof((*flag1)), 1, read_file);


		if (feof(read_file))
		{
			return 0;

		}

		fread(pos1, sizeof((*pos1)), 1, read_file);
		fread(seq_length1, sizeof((*seq_length1)), 1, read_file);
		/**********这里要改**********/
		///fread(read1, 1, (*seq_length1), read_file);
		///read1[(*seq_length1)] = '\0';
		fread(read1, sizeof(bitmapper_bs_iter),
			((*seq_length1) / 21 + ((*seq_length1) % 21 != 0)), read_file);
		///convert_3bit_to_char(tmp_buffer, (*seq_length1), read1);
		/**********这里要改**********/

		fread(flag2, sizeof((*flag2)), 1, read_file);
		fread(pos2, sizeof((*pos2)), 1, read_file);
		fread(seq_length2, sizeof((*seq_length2)), 1, read_file);
		/**********这里要改**********/
		///fread(read2, 1, (*seq_length2), read_file);
		///read2[(*seq_length2)] = '\0';
		fread(read2, sizeof(bitmapper_bs_iter),
			((*seq_length2) / 21 + ((*seq_length2) % 21 != 0)), read_file);
		///convert_3bit_to_char(tmp_buffer, (*seq_length2), read2);
		/**********这里要改**********/

		return 1;

}




inline int input_PE_alignment_multiple_threads(FILE* read_file,
	uint16_t* flag1, bitmapper_bs_iter* pos1, int* seq_length1, buffer_3_bit* read1,
	uint16_t* flag2, bitmapper_bs_iter* pos2, int* seq_length2, buffer_3_bit* read2, int thread_id)
{

	buffer_3_bit tmp_swap;

	
	pthread_mutex_lock(&methy_input_buffer.Mutex[thread_id]);

	///如果缓冲区为0且read还没从文件里读完
	///则消费者要等待，并且通知生成者进行生产
	while (methy_input_buffer.intervals[thread_id].current_size == 0 && methy_input_buffer.end == 0)
	{
		/**
		fprintf(stderr, "R:current_size: %llu\n", methy_input_buffer.intervals[thread_id].current_size);
		fprintf(stderr, "R:total_size: %llu\n", methy_input_buffer.intervals[thread_id].total_size);
		fprintf(stderr, "R:end: %llu\n", methy_input_buffer.end);
		fflush(stderr);
		**/

		///按道理这个信号量似乎不用发
		///因为队列不可能一边满一边空
		pthread_cond_signal(&methy_input_buffer.flushCond[thread_id]);
		pthread_cond_wait(&methy_input_buffer.stallCond[thread_id], &methy_input_buffer.Mutex[thread_id]);
	}

	if (methy_input_buffer.intervals[thread_id].current_size > 0)
	{

		methy_input_buffer.intervals[thread_id].current_size--;


		/************************read1*****************************/
		(*flag1) = methy_input_buffer.intervals[thread_id].R[0].r_length[methy_input_buffer.intervals[thread_id].current_size];
		(*pos1) = methy_input_buffer.intervals[thread_id].R[0].sites[methy_input_buffer.intervals[thread_id].current_size];
		(*seq_length1) = methy_input_buffer.intervals[thread_id].R[0].r_real_length[methy_input_buffer.intervals[thread_id].current_size];


		tmp_swap.buffer = methy_input_buffer.intervals[thread_id].R[0].reads_3_bit[methy_input_buffer.intervals[thread_id].current_size];
		tmp_swap.size = methy_input_buffer.intervals[thread_id].R[0].r_size_3_bit[methy_input_buffer.intervals[thread_id].current_size];

		methy_input_buffer.intervals[thread_id].R[0].reads_3_bit[methy_input_buffer.intervals[thread_id].current_size] = read1->buffer;
		methy_input_buffer.intervals[thread_id].R[0].r_size_3_bit[methy_input_buffer.intervals[thread_id].current_size] = read1->size;

		*read1 = tmp_swap;
		/************************read1*****************************/

		/************************read2*****************************/
		(*flag2) = methy_input_buffer.intervals[thread_id].R[1].r_length[methy_input_buffer.intervals[thread_id].current_size];
		(*pos2) = methy_input_buffer.intervals[thread_id].R[1].sites[methy_input_buffer.intervals[thread_id].current_size];
		(*seq_length2) = methy_input_buffer.intervals[thread_id].R[1].r_real_length[methy_input_buffer.intervals[thread_id].current_size];

		tmp_swap.buffer = methy_input_buffer.intervals[thread_id].R[1].reads_3_bit[methy_input_buffer.intervals[thread_id].current_size];
		tmp_swap.size = methy_input_buffer.intervals[thread_id].R[1].r_size_3_bit[methy_input_buffer.intervals[thread_id].current_size];

		methy_input_buffer.intervals[thread_id].R[1].reads_3_bit[methy_input_buffer.intervals[thread_id].current_size] = read2->buffer;
		methy_input_buffer.intervals[thread_id].R[1].r_size_3_bit[methy_input_buffer.intervals[thread_id].current_size] = read2->size;


		*read2 = tmp_swap;
		/************************read2*****************************/


		/*******************代码在这里*********************/

		pthread_cond_signal(&methy_input_buffer.flushCond[thread_id]);
		pthread_mutex_unlock(&methy_input_buffer.Mutex[thread_id]);

		return 1;

	}
	else
	{

		pthread_cond_signal(&methy_input_buffer.stallCond[thread_id]);
		pthread_mutex_unlock(&methy_input_buffer.Mutex[thread_id]);

		return 0;
	}



	return 1;

}



inline int input_single_alignment_multiple_threads(FILE* read_file,
	uint16_t* flag, bitmapper_bs_iter* pos, int* seq_length, buffer_3_bit* read,int thread_id)
{

	buffer_3_bit tmp_swap;


	pthread_mutex_lock(&methy_input_buffer.Mutex[thread_id]);

	///如果缓冲区为0且read还没从文件里读完
	///则消费者要等待，并且通知生成者进行生产
	while (methy_input_buffer.intervals[thread_id].current_size == 0 && methy_input_buffer.end == 0)
	{
		///按道理这个信号量似乎不用发
		///因为队列不可能一边满一边空
		pthread_cond_signal(&methy_input_buffer.flushCond[thread_id]);
		pthread_cond_wait(&methy_input_buffer.stallCond[thread_id], &methy_input_buffer.Mutex[thread_id]);
	}

	if (methy_input_buffer.intervals[thread_id].current_size > 0)
	{

		methy_input_buffer.intervals[thread_id].current_size--;


		/************************read1*****************************/
		(*flag) = methy_input_buffer.intervals[thread_id].R[0].r_length[methy_input_buffer.intervals[thread_id].current_size];
		(*pos) = methy_input_buffer.intervals[thread_id].R[0].sites[methy_input_buffer.intervals[thread_id].current_size];
		(*seq_length) = methy_input_buffer.intervals[thread_id].R[0].r_real_length[methy_input_buffer.intervals[thread_id].current_size];


		tmp_swap.buffer = methy_input_buffer.intervals[thread_id].R[0].reads_3_bit[methy_input_buffer.intervals[thread_id].current_size];
		tmp_swap.size = methy_input_buffer.intervals[thread_id].R[0].r_size_3_bit[methy_input_buffer.intervals[thread_id].current_size];

		methy_input_buffer.intervals[thread_id].R[0].reads_3_bit[methy_input_buffer.intervals[thread_id].current_size] = read->buffer;
		methy_input_buffer.intervals[thread_id].R[0].r_size_3_bit[methy_input_buffer.intervals[thread_id].current_size] = read->size;

		*read = tmp_swap;
		/************************read1*****************************/

		


		/*******************代码在这里*********************/

		pthread_cond_signal(&methy_input_buffer.flushCond[thread_id]);
		pthread_mutex_unlock(&methy_input_buffer.Mutex[thread_id]);

		return 1;

	}
	else
	{

		pthread_cond_signal(&methy_input_buffer.stallCond[thread_id]);
		pthread_mutex_unlock(&methy_input_buffer.Mutex[thread_id]);

		return 0;
	}



	return 1;

}




///1是拿到一个有效的alignment
///0是读到文件末尾了
///-1是这个alignment重复了, 拿到了一个空的alignment
////2是异常情况, 就是读了个5'-end位置有问题的结果
inline int get_alignment_PE(FILE* read_file, char* ref_genome, bitmapper_bs_iter start_pos, 
	bitmapper_bs_iter* read1, int* g_seq_length1, int* g_C_or_G1, int* g_direction1, bitmapper_bs_iter* g_pos1,
	bitmapper_bs_iter* read2, int* g_seq_length2, int* g_C_or_G2, int* g_direction2, bitmapper_bs_iter* g_pos2,
	deduplicate_PE* PE_de, int PE_distance)
{


	uint16_t flag1;
	bitmapper_bs_iter pos1;
	int C_or_G1;
	int direction1;
	int seq_length1;


	uint16_t flag2;
	bitmapper_bs_iter pos2;
	int C_or_G2;
	int direction2;
	int seq_length2;

	int return_value;

	if (input_PE_alignment(read_file,
		&flag1, &pos1, &seq_length1, read1,
		&flag2, &pos2, &seq_length2, read2))
	{


		get_direction(flag1, &C_or_G1, &direction1);
		get_direction(flag2, &C_or_G2, &direction2);
	

		return_value = remove_dup_PE(ref_genome, direction1, pos1 - start_pos, seq_length1,
			direction2, pos2 - start_pos, seq_length2, PE_de, PE_distance);


		///如果重复了
		if (!return_value)
		{
			return -1;
		}
		
		(*g_C_or_G1) = C_or_G1;
		(*g_direction1) = direction1;
		(*g_seq_length1) = seq_length1;
		(*g_pos1) = pos1;


		(*g_C_or_G2) = C_or_G2;
		(*g_direction2) = direction2;
		(*g_seq_length2) = seq_length2;
		(*g_pos2) = pos2;

		
	
		///return 1;
		return return_value;
	}
	else
	{


		return 0;
	}

}





///1是拿到一个有效的alignment
///0是读到文件末尾了
///-1是这个alignment重复了, 拿到了一个空的alignment
////2是异常情况, 就是读了个5'-end位置有问题的结果
inline int get_alignment_PE_multiple_threads(FILE* read_file, char* ref_genome, bitmapper_bs_iter start_pos,
	buffer_3_bit* read1, int* g_seq_length1, int* g_C_or_G1, int* g_direction1, bitmapper_bs_iter* g_pos1,
	buffer_3_bit* read2, int* g_seq_length2, int* g_C_or_G2, int* g_direction2, bitmapper_bs_iter* g_pos2,
	deduplicate_PE* PE_de, int PE_distance, int thread_id)
{


	uint16_t flag1;
	bitmapper_bs_iter pos1;
	int C_or_G1;
	int direction1;
	int seq_length1;


	uint16_t flag2;
	bitmapper_bs_iter pos2;
	int C_or_G2;
	int direction2;
	int seq_length2;

	int return_value;

	if (input_PE_alignment_multiple_threads(read_file,
		&flag1, &pos1, &seq_length1, read1,
		&flag2, &pos2, &seq_length2, read2, thread_id))
	{


		get_direction(flag1, &C_or_G1, &direction1);
		get_direction(flag2, &C_or_G2, &direction2);


		return_value = remove_dup_PE(ref_genome, direction1, pos1 - start_pos, seq_length1,
			direction2, pos2 - start_pos, seq_length2, PE_de, PE_distance);


		///如果重复了
		if (!return_value)
		{
			return -1;
		}

		(*g_C_or_G1) = C_or_G1;
		(*g_direction1) = direction1;
		(*g_seq_length1) = seq_length1;
		(*g_pos1) = pos1;


		(*g_C_or_G2) = C_or_G2;
		(*g_direction2) = direction2;
		(*g_seq_length2) = seq_length2;
		(*g_pos2) = pos2;



		///return 1;
		return return_value;
	}
	else
	{


		return 0;
	}

}






///1是拿到一个有效的alignment
///0是读到文件末尾了
///-1是这个alignment重复了, 拿到了一个空的alignment
inline int get_alignment_single_multiple_threads(FILE* read_file, char* ref_genome, bitmapper_bs_iter start_pos,
	buffer_3_bit* read, int* g_seq_length, int* g_C_or_G, int* g_direction, bitmapper_bs_iter* g_pos, int thread_id)
{


	uint16_t flag;
	bitmapper_bs_iter pos;
	int C_or_G;
	int direction;
	int seq_length;



	int return_value;


	if (input_single_alignment_multiple_threads(read_file,
		&flag, &pos, &seq_length, read, thread_id))
	{


		get_direction(flag, &C_or_G, &direction);



		if (!remove_dup(ref_genome, direction, pos - start_pos, seq_length))
		{
			return -1;
		}


		(*g_C_or_G) = C_or_G;
		(*g_direction) = direction;
		(*g_seq_length) = seq_length;
		(*g_pos) = pos;




		return 1;
	}
	else
	{


		return 0;
	}

}











inline bitmapper_bs_iter calculate_end5_pos(int direction, bitmapper_bs_iter pos, int seq_length)
{


	if (direction == 2)
	{
		return pos + seq_length - 1;
	}
	else   ///正向链
	{
		return pos;
	}

}

/**
///0就是重复了, 1就是没重复
inline int remove_dup_PE(char* ref_genome, deduplicate_PE* de_PE, 
	int direction1, bitmapper_bs_iter pos1, int seq_length1,
	int direction2, bitmapper_bs_iter pos2, int seq_length2)
{

	///return 1;

	bitmapper_bs_iter end_5_pos_1;
	bitmapper_bs_iter end_5_pos_2;

	end_5_pos_1 = calculate_end5_pos(direction1, pos1, seq_length1);
	end_5_pos_2 = calculate_end5_pos(direction2, pos2, seq_length2);

	bitmapper_bs_iter pair_distance;

	if (end_5_pos_1 <= end_5_pos_2)   ///检查read1的mask
	{

		pair_distance = end_5_pos_2 - end_5_pos_1;

		///说明这是rc链
		if (direction1 == 2)
		{
			///这就说明这个是重复
			if ((ref_genome[end_5_pos_1] & RC_DUP_MASK)
				&&
				(get_deduplicate_PE_mask(de_PE, de_PE->read1_bits, end_5_pos_1)))
			{
				return 0;
			}

			ref_genome[end_5_pos_1] = ref_genome[end_5_pos_1] | RC_DUP_MASK;
			set_deduplicate_PE_mask(de_PE, de_PE->read1_bits, end_5_pos_1);

			return 1;

		}
		else   ///正向链
		{
			///这就说明这个是重复
			if ((ref_genome[end_5_pos_1] & FORWARD_DUP_MASK)
				&&
				(get_deduplicate_PE_mask(de_PE, de_PE->read1_bits, end_5_pos_1)))
			{
				return 0;
			}

			ref_genome[end_5_pos_1] = ref_genome[end_5_pos_1] | FORWARD_DUP_MASK;
			set_deduplicate_PE_mask(de_PE, de_PE->read1_bits, end_5_pos_1);

			return 1;
		}



	}
	else   ///检查read2的mask
	{

	}



	return 1;
}



///1是拿到一个有效的alignment
///0是读到文件末尾了
///-1是这个alignment重复了, 拿到了一个空的alignment
inline int get_alignment_PE(FILE* read_file, char* ref_genome, bitmapper_bs_iter start_pos, 
	char* read1,int* g_seq_length1, int* g_C_or_G1, int* g_direction1, bitmapper_bs_iter* g_pos1,
	char* read2, int* g_seq_length2, int* g_C_or_G2, int* g_direction2, bitmapper_bs_iter* g_pos2,
	deduplicate_PE* de_PE)
{


	bitmapper_bs_iter flag1;
	bitmapper_bs_iter pos1;

	bitmapper_bs_iter flag2;
	bitmapper_bs_iter pos2;

	int C_or_G1;
	int direction1;
	int seq_length1;

	int C_or_G2;
	int direction2;
	int seq_length2;


	if (fscanf(read_file, "%llu\t%llu\t%s\t%llu\t%llu\t%s\n", &flag1, &pos1, read1, &flag2, &pos2, read2) != EOF)
	{
		get_direction(flag1, &C_or_G1, &direction1);
		seq_length1 = strlen(read1);

		get_direction(flag2, &C_or_G2, &direction2);
		seq_length2 = strlen(read2);

		///如果重复了
		if (!remove_dup(ref_genome, direction, pos - start_pos, seq_length))
		{
			return -1;
		}

		(*g_C_or_G) = C_or_G;
		(*g_direction) = direction;
		(*g_seq_length) = seq_length;
		(*g_pos) = pos;
		return 1;
	}
	else
	{
		return 0;
	}

}


**/



inline void update_methylation(int C_or_G, char* ref_genome, bitmapper_bs_iter* read, int read_length,
	uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total)
{


	char match, convert, rc_match, rc_convert;
	Get_BS_conversion(C_or_G);
	///bitmapper_bs_iter i;
	int i = 0;
	int word_21_index = 0;
	int word_21_inner_index = 0;
	bitmapper_bs_iter mask = 7;
	bitmapper_bs_iter read_tmp;
	bitmapper_bs_iter ref_tmp;


	for (i = 0; i < read_length; i++)
	{
		read_tmp = read[word_21_index] >> (word_21_inner_index * 3);
		read_tmp = read_tmp & 7;
		ref_tmp = ref_genome[i] & POS_MASK;

		/**
		if (debug_pos + i == 976914)
		{
			fprintf(stderr, "debug_pos: %llu\n", debug_pos);
			fprintf(stderr, "i: %llu\n", i);
		}
		**/

		if (ref_tmp == match)
		{
			if (read_tmp == convert && total[i] < 65535)
			{
				total[i]++;
			}
			else if (read_tmp == match && total[i] < 65535)
			{
				total[i]++;
				methy[i]++;
			}
		}
		else if (ref_tmp == rc_match)
		{
			if (read_tmp == rc_convert && total[i] < 65535)
			{
				correct_total[i]++;
			}
			else if (read_tmp == rc_match && total[i] < 65535)
			{
				correct_total[i]++;
				correct_methy[i]++;
			}
		}


		word_21_inner_index++;

		if (word_21_inner_index == 21)
		{
			word_21_inner_index = 0;
			word_21_index++;
		}
	}

}







inline void get_chrome_position(bitmapper_bs_iter* position, char** chrome_name, int current_chrome_id, int* chrome_id)
{
	bitmapper_bs_iter map_location = *position;

	int now_ref_name = 0;


	for (now_ref_name = current_chrome_id; now_ref_name < refChromeCont; ++now_ref_name)
	{
		if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
			break;
	}

	if (now_ref_name >= refChromeCont)
	{
		for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
		{
			if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
				break;
		}
	}


	///map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;
	///bed_Graph是0-based
	map_location = map_location - _ih_refGenName[now_ref_name].start_location;


	(*position) = map_location;
	(*chrome_name) = _ih_refGenName[now_ref_name]._rg_chrome_name;
	(*chrome_id) = now_ref_name;
	///fprintf(stderr, "%s\t%llu\n", chrome_name, (*position));

}



inline int check_CpG(long long ref_pos, long long part_pos,
	int chrome_id, char* ref_genome, char direction, char* last_two_bases)
{

	if (direction == 1) ///C
	{
		///第一个判断条件是不能跨染色体
		///第二个是判断是不是CpG
		///不用担心左边会不会超过当前block的长度, 因为后面会多取一点出来(注意反方向需要考虑这个问题)
		///即使是最后一个block, 也受第一个条件限制
		if ((ref_pos + 1 <= _ih_refGenName[chrome_id].end_location)
			&&
			((ref_genome[part_pos + 1] & POS_MASK) == 2))  ///G
		{
			return 1;
		}
	}
	else if (direction == 2) ///G
	{
		///这个是防止溢出
		if (part_pos > 0)
		{

			///第一个判断条件是不能跨染色体
			///第二个是判断是不是CpG
			if ((ref_pos >= _ih_refGenName[chrome_id].start_location + 1)
				&&
				((ref_genome[part_pos - 1] & POS_MASK) == 1)) //C
			{
				return 1;
			}
		}
		else   ///这个是为了处理part_pos = 0的情况
		{
			///第一个判断条件是不能跨染色体
			///第二个是判断是不是CpG
			if ((ref_pos >= _ih_refGenName[chrome_id].start_location + 1)
				&&
				((last_two_bases[1] & POS_MASK) == 1)) //C
			{
				return 1;
			}
		}

		
	}

	return 0;
	
}




////这个必须和check_CpG连在一起用才正确
inline int check_CHG(long long ref_pos, long long part_pos,
	int chrome_id, char* ref_genome, char direction, char* last_two_bases)
{

	if (direction == 1) ///C
	{
		///第一个判断条件是不能跨染色体
		///第二个是判断是不是CpG
		///不用担心左边会不会超过当前block的长度, 因为后面会多取一点出来(注意反方向需要考虑这个问题)
		///即使是最后一个block, 也受第一个条件限制
		if ((ref_pos + 2 <= _ih_refGenName[chrome_id].end_location)
			&&
			((ref_genome[part_pos + 2] & POS_MASK) == 2))  ///G
		{
			return 1;
		}
	}
	else if (direction == 2) ///G
	{


		///这个是防止溢出, >1就至少是2
		///if (part_pos > 0)
		if (part_pos > 1)
		{
			///第一个判断条件是不能跨染色体
			///第二个是判断是不是CpG
			if ((ref_pos >= _ih_refGenName[chrome_id].start_location + 2)
				&&
				((ref_genome[part_pos - 2] & POS_MASK) == 1)) //C
			{
				return 1;
			}
		}
		else   ///这个是为了处理part_pos = 0和1的情况
		{


			///第一个判断条件是不能跨染色体
			///第二个是判断是不是CpG
			if ((ref_pos >= _ih_refGenName[chrome_id].start_location + 2)
				&&
				((last_two_bases[part_pos] & POS_MASK) == 1)) //C   ///当part_pos = 1时，应该检查last_two_bases[1]; 当part_pos = 0时，应该检查last_two_bases[0]
			{
				return 1;
			}
		}


	}

	return 0;

}


inline int get_context(long long ref_pos, long long part_pos,
	int chrome_id, char* ref_genome, char direction, char* last_two_bases)
{
	if (check_CpG(ref_pos, part_pos, chrome_id, ref_genome, direction, last_two_bases))
	{
		return 0;
	}
	else if (check_CHG(ref_pos, part_pos, chrome_id, ref_genome, direction, last_two_bases))
	{
		return 1;
	}
	else
	{
		return 2;
	}
}



inline void check_print_methylation(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	bitmapper_bs_iter start_pos, bitmapper_bs_iter end_pos, char* ref_genome, 
	char* last_two_bases, int* need_context)
{
	bitmapper_bs_iter length = end_pos - start_pos + 1;

	bitmapper_bs_iter i;
	bitmapper_bs_iter tmp_pos;
	char* chrome_name;
	int chrome_id = 0;

	unsigned char tmp;

	int current_context;

	for (i = 0; i < length; i++, start_pos++)
	{


		tmp = ref_genome[i] & POS_MASK;

		if (tmp == 1 || tmp == 2) //C或者G
		{
			/**
			///if (correct_methy[i] != correct_total[i])
			if (correct_total[i] > minVariantDepth && ((double)(correct_methy[i]) / (double)(correct_total[i])) < maxVariantFrac)
			{
				continue;
			}

			if (total[i]==0)
			{
				continue;
			}
			**/




			tmp_pos = start_pos;
			get_chrome_position(&tmp_pos, &chrome_name, chrome_id, &chrome_id);

			
			current_context = 
				get_context(start_pos, i, chrome_id, ref_genome, tmp, last_two_bases);
			
			tmp = ref_genome[i] & CONTEXT_MASK;
			tmp = tmp >> 4;

			///tmp = (ref_genome[i] & CONTEXT_MASK)>>4;

			if (tmp != current_context)
			{
				fprintf(stderr, "tmp: %llu, current_context: %llu, ref_genome[i]: %llu, start_pos: %llu, i: %llu\n", 
					tmp, current_context, ref_genome[i], start_pos, i);

				fprintf(stderr, "chrome_id: %llu, chrome_name: %s\n",
					chrome_id, chrome_name);

				fprintf(stderr, "start_location: %llu, end_location: %llu\n",
					_ih_refGenName[chrome_id].start_location, _ih_refGenName[chrome_id].end_location);
			}

			///current_context = tmp;


			/**
			if (need_context[current_context])
			{

				if (current_context == 0)
				{
					output_single_methy_CpG(tmp_pos, chrome_name, methy[i], total[i]);
				}
				else if (current_context == 1)
				{
					output_single_methy_CHG(tmp_pos, chrome_name, methy[i], total[i]);
				}
				else if (current_context == 2)
				{
					output_single_methy_CHH(tmp_pos, chrome_name, methy[i], total[i]);
				}

			}
			**/

			if (need_context[current_context])
			{

				if (current_context == 0)
				{
					output_single_methy_CpG(tmp_pos, chrome_name, methy[i], total[i]);
				}

				if (current_context == 1)
				{
					output_single_methy_CHG(tmp_pos, chrome_name, methy[i], total[i]);
				}

				if (current_context == 2)
				{
					output_single_methy_CHH(tmp_pos, chrome_name, methy[i], total[i]);
				}

			}


		}
		


		

	}
}



inline void print_methylation_back(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	bitmapper_bs_iter start_pos, bitmapper_bs_iter end_pos, char* ref_genome,
	char* last_two_bases, int* need_context)
{
	bitmapper_bs_iter length = end_pos - start_pos + 1;

	bitmapper_bs_iter i;
	bitmapper_bs_iter tmp_pos;
	char* chrome_name;
	int chrome_id = 0;

	unsigned char tmp;

	int current_context;

	for (i = 0; i < length; i++, start_pos++)
	{


		tmp = ref_genome[i] & POS_MASK;

		if (tmp == 1 || tmp == 2) //C或者G
		{
			
			///if (correct_methy[i] != correct_total[i])
			if (correct_total[i] > minVariantDepth && ((double)(correct_methy[i]) / (double)(correct_total[i])) < maxVariantFrac)
			{
				continue;
			}

			if (total[i]==0)
			{
				continue;
			}
				




			tmp_pos = start_pos;
			get_chrome_position(&tmp_pos, &chrome_name, chrome_id, &chrome_id);


			current_context =
				get_context(start_pos, i, chrome_id, ref_genome, tmp, last_two_bases);



			if (need_context[current_context])
			{

				if (current_context == 0)
				{
					output_single_methy_CpG(tmp_pos, chrome_name, methy[i], total[i]);
				}

				if (current_context == 1)
				{
					output_single_methy_CHG(tmp_pos, chrome_name, methy[i], total[i]);
				}

				if (current_context == 2)
				{
					output_single_methy_CHH(tmp_pos, chrome_name, methy[i], total[i]);
				}

			}


		}





	}
}



inline void print_sub_methylation_multiple_thread(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	bitmapper_bs_iter start_pos, bitmapper_bs_iter end_pos, char* ref_genome, int* need_context, int thread_id)
{
	bitmapper_bs_iter length = end_pos - start_pos + 1;

	bitmapper_bs_iter i;
	bitmapper_bs_iter tmp_pos;
	char* chrome_name;
	int chrome_id = 0;

	unsigned char tmp;


	for (i = 0; i < length; i++, start_pos++)
	{

		tmp = (ref_genome[i] & CONTEXT_MASK) >> 4;

		if (tmp <= 2)
		{

			if (total[i] == 0)
			{
				continue;
			}

			if (correct_total[i] > minVariantDepth && ((double)(correct_methy[i]) / (double)(correct_total[i])) < maxVariantFrac)
			{
				continue;
			}


			if (need_context[tmp])
			{

				tmp_pos = start_pos;
				get_chrome_position(&tmp_pos, &chrome_name, chrome_id, &chrome_id);

				switch (tmp)
				{
				case 0: ///CpG
					output_single_methy_multiple_thread(tmp_pos, chrome_name, methy[i], total[i], 
						&(methy_input_buffer.CpG_buffer[thread_id]));
					break;
				case 1: ///CHG
					output_single_methy_multiple_thread(tmp_pos, chrome_name, methy[i], total[i],
						&(methy_input_buffer.CHG_buffer[thread_id]));
					break;
				case 2:  ///CHH
					output_single_methy_multiple_thread(tmp_pos, chrome_name, methy[i], total[i],
						&(methy_input_buffer.CHH_buffer[thread_id]));
					break;

				}

			}
		}
	}
}



void* extract_methylation_multiple_thread(void* arg)
{
	int thread_id = *((int *)arg);
	uint16_t* methy;
	uint16_t* total;
	uint16_t* correct_methy;
	uint16_t* correct_total;
	bitmapper_bs_iter start_pos;
	bitmapper_bs_iter end_pos;
	char* ref_genome;
	bitmapper_bs_iter length;


	bitmapper_bs_iter i;
	bitmapper_bs_iter tmp_pos;
	char* chrome_name;
	int chrome_id = 0;

	unsigned char tmp;

	/**
	fprintf(stderr, "extract1\n");
	fflush(stderr);
	**/

	/**
	///先把自己阻塞住，因为此时数据没准备好
	pthread_mutex_lock(&methy_input_buffer.process_Mutex_methy[thread_id]);
	pthread_cond_wait(&methy_input_buffer.process_stallCond_methy[thread_id],
		&methy_input_buffer.process_Mutex_methy[thread_id]);
	pthread_mutex_unlock(&methy_input_buffer.process_Mutex_methy[thread_id]);
	**/
	/**
	fprintf(stderr, "extract2\n");
	fflush(stderr);
	**/

	while (methy_input_buffer.all_end == 0)
	{
		methy = methy_input_buffer.M_methy[thread_id];
		total = methy_input_buffer.M_total[thread_id];
		correct_methy = methy_input_buffer.M_correct_methy[thread_id];
		correct_total = methy_input_buffer.M_correct_total[thread_id];
		start_pos = methy_input_buffer.M_start_pos[thread_id];
		end_pos = methy_input_buffer.M_end_pos[thread_id];
		ref_genome = methy_input_buffer.M_ref_genome[thread_id];

		////这是为了防止溢出啊
		length = end_pos + 1;
		length = length - start_pos;
		////这是为了防止溢出啊

		for (i = 0; i < length; i++, start_pos++)
		{

			tmp = (ref_genome[i] & CONTEXT_MASK) >> 4;

			if (tmp <= 2)
			{

				if (total[i] == 0)
				{
					continue;
				}

				if (correct_total[i] > minVariantDepth && ((double)(correct_methy[i]) / (double)(correct_total[i])) < maxVariantFrac)
				{
					continue;
				}


				if (methy_input_buffer.need_context[tmp])
				{

					tmp_pos = start_pos;
					get_chrome_position(&tmp_pos, &chrome_name, chrome_id, &chrome_id);

					switch (tmp)
					{
					case 0: ///CpG
						output_single_methy_multiple_thread(tmp_pos, chrome_name, methy[i], total[i],
							&(methy_input_buffer.CpG_buffer[thread_id]));
						break;
					case 1: ///CHG
						output_single_methy_multiple_thread(tmp_pos, chrome_name, methy[i], total[i],
							&(methy_input_buffer.CHG_buffer[thread_id]));
						break;
					case 2:  ///CHH
						output_single_methy_multiple_thread(tmp_pos, chrome_name, methy[i], total[i],
							&(methy_input_buffer.CHH_buffer[thread_id]));
						break;

					}

				}
			}
		}

		/**
		fprintf(stderr, "extract3\n");
		fflush(stderr);
		**/

		pthread_mutex_lock(&methy_input_buffer.all_completed_Mutex_methy);
		///注意这里和methy_input_buffer.all_completed不一样, 那个变量还有个input线程控制
		///这里要先+1，再判断
		methy_input_buffer.all_completed_methy++;
		if (methy_input_buffer.all_completed_methy == THREAD_COUNT)
		{
			pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex_methy);
			pthread_cond_signal(&methy_input_buffer.main_thread_flushCond_methy);
			pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex_methy);
		}
		pthread_mutex_unlock(&methy_input_buffer.all_completed_Mutex_methy);

		/**
		fprintf(stderr, "extract4\n");
		fflush(stderr);
		**/

		pthread_mutex_lock(&methy_input_buffer.process_Mutex_methy[thread_id]);
		pthread_cond_wait(&methy_input_buffer.process_stallCond_methy[thread_id], 
			&methy_input_buffer.process_Mutex_methy[thread_id]);
		pthread_mutex_unlock(&methy_input_buffer.process_Mutex_methy[thread_id]);

		/**
		fprintf(stderr, "extract5\n");
		fflush(stderr);
		**/

	}

	/**
	fprintf(stderr, "thread_id: %d\n", thread_id);
	fflush(stderr);
	**/
}




inline void print_methylation_multiple_thread(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	long long start_pos, long long end_pos, char* ref_genome,
	char* last_two_bases, int* need_context, pthread_t *_r_threads_output, int is_update)
{
	long long length = end_pos - start_pos + 1;

	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter tmp_pos;
	char* chrome_name;
	int chrome_id = 0;

	unsigned char tmp;

	int j;
	bitmapper_bs_iter sub_step_length = 100000;
	bitmapper_bs_iter step_length = sub_step_length*THREAD_COUNT;

	long long inner_offset;


	/**
	fprintf(stderr, "******\n");

	fprintf(stderr, "i: %d\n", i);
	fprintf(stderr, "step_length: %d\n", step_length);
	**/

	while (i + step_length <= length)
	{
		/**
		fprintf(stderr, "0:inner_i; %d\n", i);
		fflush(stderr);
		**/

		///这个一定要有
		methy_input_buffer.all_completed_methy = 0;

		
		for (j = 0; j < THREAD_COUNT; j++)
		{
			inner_offset = i + j * sub_step_length;
			///准备数据
			methy_input_buffer.M_methy[j] = methy + inner_offset;
			methy_input_buffer.M_total[j] = total + inner_offset;
			methy_input_buffer.M_correct_methy[j] = correct_methy + inner_offset;
			methy_input_buffer.M_correct_total[j] = correct_total + inner_offset;
			methy_input_buffer.M_ref_genome[j] = ref_genome + inner_offset;
			methy_input_buffer.M_start_pos[j] = start_pos + inner_offset;
			methy_input_buffer.M_end_pos[j] = methy_input_buffer.M_start_pos[j] + sub_step_length - 1;

			if (is_update == 0)
			{
				int *arg = (int*)malloc(sizeof(*arg));
				*arg = j;

				pthread_create(_r_threads_output + j, NULL, extract_methylation_multiple_thread, (void*)arg);
			}
			else
			{
				///启动线程
				pthread_mutex_lock(&methy_input_buffer.process_Mutex_methy[j]);
				pthread_cond_signal(&methy_input_buffer.process_stallCond_methy[j]);
				pthread_mutex_unlock(&methy_input_buffer.process_Mutex_methy[j]);
			}

		
		}

		/**
		fprintf(stderr, "1:inner_i; %d\n", i);
		fflush(stderr);
		**/

		///阻塞自己
		pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex_methy);
		pthread_cond_wait(&methy_input_buffer.main_thread_flushCond_methy, &methy_input_buffer.main_thread_Mutex_methy);
		pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex_methy);

		/**
		fprintf(stderr, "2:inner_i; %d\n", i);
		fflush(stderr);
		**/

		///输出
		for (j = 0; j < THREAD_COUNT; j++)
		{
			output_methy_directly(&(methy_input_buffer.CpG_buffer[j]), &(methy_input_buffer.CHG_buffer[j]),
				&(methy_input_buffer.CHH_buffer[j]));
		}
		

		i = i + step_length;
	}

	/**
	fprintf(stderr, "i: %d\n", i);
	fprintf(stderr, "length: %d\n", length);
	**/

	/**
	if (i < length)
	{
		start_pos = start_pos + i;

		print_sub_methylation_multiple_thread(methy + i, total + i, correct_methy + i, correct_total + i,
			start_pos, start_pos + length -i - 1, ref_genome + i, need_context, 0);

		output_methy_directly(&(methy_input_buffer.CpG_buffer[0]), &(methy_input_buffer.CHG_buffer[0]),
			&(methy_input_buffer.CHH_buffer[0]));
	}
	**/
	///不管还剩多少数据，都得启动THREAD_COUNT这么多线程
	if (i < length)
	{
		///这个一定要有
		methy_input_buffer.all_completed_methy = 0;



		for (j = 0; j < THREAD_COUNT; j++)
		{
			inner_offset = i + j * sub_step_length;

			if (inner_offset < length)
			{
				///准备数据
				methy_input_buffer.M_methy[j] = methy + inner_offset;
				methy_input_buffer.M_total[j] = total + inner_offset;
				methy_input_buffer.M_correct_methy[j] = correct_methy + inner_offset;
				methy_input_buffer.M_correct_total[j] = correct_total + inner_offset;
				methy_input_buffer.M_ref_genome[j] = ref_genome + inner_offset;
				methy_input_buffer.M_start_pos[j] = start_pos + inner_offset;

				if (length - inner_offset >= sub_step_length)
				{
					methy_input_buffer.M_end_pos[j] = methy_input_buffer.M_start_pos[j] + sub_step_length - 1;
				}
				else
				{
					methy_input_buffer.M_end_pos[j] = methy_input_buffer.M_start_pos[j] + length - inner_offset - 1;
				}
				

				if (is_update == 0)
				{
					int *arg = (int*)malloc(sizeof(*arg));
					*arg = j;

					pthread_create(_r_threads_output + j, NULL, extract_methylation_multiple_thread, (void*)arg);
				}
				else
				{
					///启动线程
					pthread_mutex_lock(&methy_input_buffer.process_Mutex_methy[j]);
					pthread_cond_signal(&methy_input_buffer.process_stallCond_methy[j]);
					pthread_mutex_unlock(&methy_input_buffer.process_Mutex_methy[j]);
				}
			}
			else
			{
				///虽然这个线程分不到数据，但是也得启动，否则会死锁
				///置0和置1就是说明长度为0
				methy_input_buffer.M_start_pos[j] = 1;
				methy_input_buffer.M_end_pos[j] = 0;

				if (is_update == 0)
				{
					int *arg = (int*)malloc(sizeof(*arg));
					*arg = j;

					pthread_create(_r_threads_output + j, NULL, extract_methylation_multiple_thread, (void*)arg);
				}
				else
				{
					///启动线程
					pthread_mutex_lock(&methy_input_buffer.process_Mutex_methy[j]);
					pthread_cond_signal(&methy_input_buffer.process_stallCond_methy[j]);
					pthread_mutex_unlock(&methy_input_buffer.process_Mutex_methy[j]);
				}


			}

			


		}

		/**
		fprintf(stderr, "1:inner_i; %d\n", i);
		fflush(stderr);
		**/

		///阻塞自己
		pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex_methy);
		pthread_cond_wait(&methy_input_buffer.main_thread_flushCond_methy, &methy_input_buffer.main_thread_Mutex_methy);
		pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex_methy);

		/**
		fprintf(stderr, "2:inner_i; %d\n", i);
		fflush(stderr);
		**/

		///输出
		for (j = 0; j < THREAD_COUNT; j++)
		{
			output_methy_directly(&(methy_input_buffer.CpG_buffer[j]), &(methy_input_buffer.CHG_buffer[j]),
				&(methy_input_buffer.CHH_buffer[j]));
		}
	}


}




inline void print_methylation(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	bitmapper_bs_iter start_pos, bitmapper_bs_iter end_pos, char* ref_genome,
	char* last_two_bases, int* need_context)
{
	bitmapper_bs_iter length = end_pos - start_pos + 1;

	bitmapper_bs_iter i;
	bitmapper_bs_iter tmp_pos;
	char* chrome_name;
	int chrome_id = 0;

	unsigned char tmp;


	for (i = 0; i < length; i++, start_pos++)
	{

		tmp = (ref_genome[i] & CONTEXT_MASK) >> 4;

		if (tmp <= 2)
		{

			if (total[i] == 0)
			{
				continue;
			}

			if (correct_total[i] > minVariantDepth && ((double)(correct_methy[i]) / (double)(correct_total[i])) < maxVariantFrac)
			{
				continue;
			}


			if (need_context[tmp])
			{

				tmp_pos = start_pos;
				get_chrome_position(&tmp_pos, &chrome_name, chrome_id, &chrome_id);

				switch (tmp)
				{
				case 0: ///CpG
					output_single_methy_CpG(tmp_pos, chrome_name, methy[i], total[i]);
					break;
				case 1: ///CHG
					output_single_methy_CHG(tmp_pos, chrome_name, methy[i], total[i]);
					break;
				case 2:  ///CHH
					output_single_methy_CHH(tmp_pos, chrome_name, methy[i], total[i]);
					break;

				}

			}
		}
	}
}




inline void check_context(uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	bitmapper_bs_iter start_pos, bitmapper_bs_iter end_pos, char* ref_genome,
	char* last_two_bases, int* need_context)
{
	bitmapper_bs_iter length = end_pos - start_pos + 1;

	bitmapper_bs_iter i;
	bitmapper_bs_iter tmp_pos;
	char* chrome_name;
	int chrome_id = 0;

	char tmp;

	int current_context;

	for (i = 0; i < length; i++, start_pos++)
	{


		tmp = ref_genome[i] & POS_MASK;

		if (tmp == 1 || tmp == 2) //C或者G
		{
			




			tmp_pos = start_pos;
			get_chrome_position(&tmp_pos, &chrome_name, chrome_id, &chrome_id);


			current_context =
				get_context(start_pos, i, chrome_id, ref_genome, tmp, last_two_bases);



			tmp = (ref_genome[i] & CONTEXT_MASK) >> 4;

			if (tmp != current_context)
			{
				fprintf(stderr, "check***** tmp: %llu, current_context: %llu, ref_genome[i]: %llu\n", tmp, current_context, ref_genome[i]);
			}

			
		}

	}
}


inline int get_abnormal_PE_alignment(bitmapper_bs_iter* read1, bitmapper_bs_iter* read2, int* C_or_G1, int* C_or_G2,
	int* direction1, int* direction2, bitmapper_bs_iter* pos1, bitmapper_bs_iter* pos2, int* read1_length, int* read2_length,
	FILE* abnormal_file)
{
	fread(read1_length, sizeof(int), 1, abnormal_file);
	if (feof(abnormal_file))
	{
		return 0;
	}
	fread(C_or_G1, sizeof(int), 1, abnormal_file);
	fread(direction1, sizeof(int), 1, abnormal_file);
	fread(pos1, sizeof(bitmapper_bs_iter), 1, abnormal_file);
	fread(read1, sizeof(bitmapper_bs_iter),
		((*read1_length) / 21 + ((*read1_length) % 21 != 0)), abnormal_file);

	fread(read2_length, sizeof(int), 1, abnormal_file);
	fread(C_or_G2, sizeof(int), 1, abnormal_file);
	fread(direction2, sizeof(int), 1, abnormal_file);
	fread(pos2, sizeof(bitmapper_bs_iter), 1, abnormal_file);
	fread(read2, sizeof(bitmapper_bs_iter),
		((*read2_length) / 21 + ((*read2_length) % 21 != 0)), abnormal_file);

	return 1;
}



void phrase_abnormal_alignment(char* ref_genome, bitmapper_bs_iter start_pos,
	bitmapper_bs_iter* read1, bitmapper_bs_iter* read2, FILE* abnormal_file, deduplicate_PE* PE_de, int PE_distance,
	uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	long long* unique_read, long long* duplicate_read)
{
	int seq_length1, seq_length2;
	int C_or_G1, C_or_G2;
	int direction1, direction2;
	bitmapper_bs_iter pos1, pos2;
	int read1_length, read2_length;


	


	


	while (get_abnormal_PE_alignment(read1, read2, &C_or_G1, &C_or_G2,
		&direction1, &direction2, &pos1, &pos2, &read1_length, &read2_length,
		abnormal_file))
	{
		
		unset_abnormal_dup_PE(ref_genome, direction1, pos1 - start_pos, seq_length1,
			direction2, pos2 - start_pos, seq_length2, PE_de, PE_distance);
		
	}



	///回到文件开头
	fseek(abnormal_file, 0, SEEK_SET);


	while (get_abnormal_PE_alignment(read1, read2, &C_or_G1, &C_or_G2,
		&direction1, &direction2, &pos1, &pos2, &read1_length, &read2_length,
		abnormal_file))
	{

		if (remove_abnormal_dup_PE(ref_genome, direction1, pos1 - start_pos, seq_length1,
			direction2, pos2 - start_pos, seq_length2, PE_de, PE_distance))
		{
			update_methylation(C_or_G1, ref_genome + pos1 - start_pos, read1, read1_length,
				methy + pos1 - start_pos, total + pos1 - start_pos,
				correct_methy + pos1 - start_pos, correct_total + pos1 - start_pos);

			update_methylation(C_or_G2, ref_genome + pos2 - start_pos, read2, read2_length,
				methy + pos2 - start_pos, total + pos2 - start_pos,
				correct_methy + pos2 - start_pos, correct_total + pos2 - start_pos);

			(*unique_read)++;
		}
		else
		{
			(*duplicate_read)++;
		}
	}

}



inline int get_overlap_alignment(bitmapper_bs_iter* read, int* C_or_G, 
	int* direction, bitmapper_bs_iter* pos, int* read_length, FILE* overlap_file)
{
	fread(read_length, sizeof(int), 1, overlap_file);
	if (feof(overlap_file))
	{
		return 0;
	}
	fread(C_or_G, sizeof(int), 1, overlap_file);
	fread(direction, sizeof(int), 1, overlap_file);
	fread(pos, sizeof(bitmapper_bs_iter), 1, overlap_file);
	fread(read, sizeof(bitmapper_bs_iter),
		((*read_length) / 21 + ((*read_length) % 21 != 0)), overlap_file);

	return 1;
}


void phrase_overlap_alignment(char* ref_genome, bitmapper_bs_iter start_pos,
	bitmapper_bs_iter* read, FILE* overlap_file, 
	uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total, long long* overlap_read)
{
	int C_or_G;
	int direction;
	bitmapper_bs_iter pos;
	int read_length;








	while (get_overlap_alignment(read, &C_or_G, &direction, &pos, &read_length, overlap_file))
	{

		update_methylation(C_or_G, ref_genome + pos - start_pos, read, read_length,
			methy + pos - start_pos, total + pos - start_pos,
			correct_methy + pos - start_pos, correct_total + pos - start_pos);

		(*overlap_read)++;

	}

}




inline void write_abnormal_PE_alignment(
	bitmapper_bs_iter* read1, int read1_length, int C_or_G1, int direction1, bitmapper_bs_iter pos1,
	bitmapper_bs_iter* read2, int read2_length, int C_or_G2, int direction2, bitmapper_bs_iter pos2,
	FILE* abnormal_file)
{
	fwrite(&read1_length, sizeof(int), 1, abnormal_file);
	fwrite(&C_or_G1, sizeof(int), 1, abnormal_file);
	fwrite(&direction1, sizeof(int), 1, abnormal_file);
	fwrite(&pos1, sizeof(bitmapper_bs_iter), 1, abnormal_file);
	fwrite(read1, sizeof(bitmapper_bs_iter),
		(read1_length / 21 + (read1_length % 21 != 0)), abnormal_file);




	fwrite(&read2_length, sizeof(int), 1, abnormal_file);
	fwrite(&C_or_G2, sizeof(int), 1, abnormal_file);
	fwrite(&direction2, sizeof(int), 1, abnormal_file);
	fwrite(&pos2, sizeof(bitmapper_bs_iter), 1, abnormal_file);
	fwrite(read2, sizeof(bitmapper_bs_iter),
		(read2_length / 21 + (read2_length % 21 != 0)), abnormal_file);
}


inline void write_overlap_alignment(
	bitmapper_bs_iter* read, int read_length, int C_or_G, int direction, bitmapper_bs_iter pos,
	FILE* overlap_file)
{
	fwrite(&read_length, sizeof(int), 1, overlap_file);
	fwrite(&C_or_G, sizeof(int), 1, overlap_file);
	fwrite(&direction, sizeof(int), 1, overlap_file);
	fwrite(&pos, sizeof(bitmapper_bs_iter), 1, overlap_file);
	fwrite(read, sizeof(bitmapper_bs_iter),
		(read_length / 21 + (read_length % 21 != 0)), overlap_file);
}


/**
void process_PE_methy_alignment(FILE* read_file, FILE* abnormal_file, char* ref_genome, 
	uint16_t* methy, uint16_t* total, uint16_t* correct_methy, uint16_t* correct_total,
	bitmapper_bs_iter start_pos, int thread_id, deduplicate_PE* PE_de, int PE_distance, 
	long long* total_read, long long* duplicate_read, long long* unique_read, long long* abnormal_read)**/
void process_PE_methy_alignment(int thread_id)
{

	int return_value;

	buffer_3_bit read1;
	buffer_3_bit read2;

	Init_buffer_3_bit(read1);
	Init_buffer_3_bit(read2);

	int read1_length;
	int read2_length;

	int C_or_G1;
	int C_or_G2;
	int direction1;
	int direction2;

	bitmapper_bs_iter pos1;
	bitmapper_bs_iter pos2;



	///1是拿到一个有效的alignment
	///0是读到文件末尾了
	///-1是这个alignment重复了, 拿到了一个空的alignment
	////2是异常情况, 就是读了个5'-end位置有问题的结果
	while (return_value = get_alignment_PE_multiple_threads(
		parameters.read_file, parameters.ref_genome, parameters.start_pos,
		&read1, &read1_length, &C_or_G1, &direction1, &pos1,
		&read2, &read2_length, &C_or_G2, &direction2, &pos2,
		&(parameters.PE_de), parameters.PE_distance, thread_id))
	{

		parameters.total_read[thread_id]++;

		if (return_value == -1)
		{
			parameters.duplicate_read[thread_id]++;
			continue;
		}
		else if (return_value == 1)
		{
			parameters.unique_read[thread_id]++;
		}
		else   ///这是返回值为2的情况，也是5'端位置异常的情况  ///暂时还没处理
		{

			parameters.abnormal_read[thread_id]++;


			pthread_mutex_lock(&methy_input_buffer.abnormal_file_Mutex);

			///将异常结果写入临时文件
			write_abnormal_PE_alignment(
				read1.buffer, read1_length, C_or_G1, direction1, pos1,
				read2.buffer, read2_length, C_or_G2, direction2, pos2,
				parameters.abnormal_file);

			pthread_mutex_unlock(&methy_input_buffer.abnormal_file_Mutex);

		}


		update_methylation(C_or_G1, parameters.ref_genome + pos1 - parameters.start_pos, read1.buffer, read1_length,
			parameters.methy + pos1 - parameters.start_pos, parameters.total + pos1 - parameters.start_pos,
			parameters.correct_methy + pos1 - parameters.start_pos, parameters.correct_total + pos1 - parameters.start_pos);

		update_methylation(C_or_G2, parameters.ref_genome + pos2 - parameters.start_pos, read2.buffer, read2_length,
			parameters.methy + pos2 - parameters.start_pos, parameters.total + pos2 - parameters.start_pos,
			parameters.correct_methy + pos2 - parameters.start_pos, parameters.correct_total + pos2 - parameters.start_pos);



	}
}


void* process_PE_methy_alignment_multiple_threads(void* arg)
{
	int thread_id = *((int *)arg);
	int return_value;

	buffer_3_bit read1;
	buffer_3_bit read2;

	Init_buffer_3_bit(read1);
	Init_buffer_3_bit(read2);

	int read1_length;
	int read2_length;

	int C_or_G1;
	int C_or_G2;
	int direction1;
	int direction2;

	bitmapper_bs_iter pos1;
	bitmapper_bs_iter pos2;

	bitmapper_bs_iter key_pos;

	


	while (methy_input_buffer.all_end == 0)
	{

		///1是拿到一个有效的alignment
		///0是读到文件末尾了
		///-1是这个alignment重复了, 拿到了一个空的alignment
		////2是异常情况, 就是读了个5'-end位置有问题的结果
		while (return_value = get_alignment_PE_multiple_threads(
			parameters.read_file, parameters.ref_genome, parameters.start_pos,
			&read1, &read1_length, &C_or_G1, &direction1, &pos1,
			&read2, &read2_length, &C_or_G2, &direction2, &pos2,
			&(parameters.PE_de), parameters.PE_distance, thread_id))
		{

			parameters.total_read[thread_id]++;

			if (return_value == -1)
			{
				parameters.duplicate_read[thread_id]++;
				continue;
			}
			else if (return_value == 1)
			{
				parameters.unique_read[thread_id]++;
			}
			else   ///这是返回值为2的情况，也是5'端位置异常的情况  ///暂时还没处理
			{

				parameters.abnormal_read[thread_id]++;


				pthread_mutex_lock(&methy_input_buffer.abnormal_file_Mutex);

				///将异常结果写入临时文件
				write_abnormal_PE_alignment(
					read1.buffer, read1_length, C_or_G1, direction1, pos1,
					read2.buffer, read2_length, C_or_G2, direction2, pos2,
					parameters.abnormal_file);

				pthread_mutex_unlock(&methy_input_buffer.abnormal_file_Mutex);

				continue;

			}



			key_pos = pos1 <= pos2 ? pos1 : pos2;
			if (key_pos > methy_input_buffer.each_buffer_interval_end[thread_id] || key_pos < methy_input_buffer.each_buffer_interval_start[thread_id])
			{
				fprintf(stderr, "ERROR: process_PE_methy_alignment_multiple_threads\n");
				exit(0);
			}


			if (pos1 + read1_length - 1 <= methy_input_buffer.each_buffer_interval_end[thread_id])
			{
				update_methylation(C_or_G1, parameters.ref_genome + pos1 - parameters.start_pos, read1.buffer, read1_length,
					parameters.methy + pos1 - parameters.start_pos, parameters.total + pos1 - parameters.start_pos,
					parameters.correct_methy + pos1 - parameters.start_pos, parameters.correct_total + pos1 - parameters.start_pos);
			}
			else
			{
				pthread_mutex_lock(&methy_input_buffer.overlap_file_Mutex);

				write_overlap_alignment(read1.buffer, read1_length, C_or_G1, direction1, pos1,
					parameters.overlap_file);

				pthread_mutex_unlock(&methy_input_buffer.overlap_file_Mutex);
			}
			
			if (pos2 + read2_length - 1 <= methy_input_buffer.each_buffer_interval_end[thread_id])
			{
				update_methylation(C_or_G2, parameters.ref_genome + pos2 - parameters.start_pos, read2.buffer, read2_length,
					parameters.methy + pos2 - parameters.start_pos, parameters.total + pos2 - parameters.start_pos,
					parameters.correct_methy + pos2 - parameters.start_pos, parameters.correct_total + pos2 - parameters.start_pos);
			}
			else
			{
				pthread_mutex_lock(&methy_input_buffer.overlap_file_Mutex);

				write_overlap_alignment(read2.buffer, read2_length, C_or_G2, direction2, pos2,
					parameters.overlap_file);

				pthread_mutex_unlock(&methy_input_buffer.overlap_file_Mutex);
			}
			



		}

		pthread_mutex_lock(&methy_input_buffer.all_completed_Mutex);

		///这个是先判断，然后再+1，是因为除了这些处理线程之外，还有个input线程也要+1
		if (methy_input_buffer.all_completed == THREAD_COUNT)
		{

			pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex);
			pthread_cond_signal(&methy_input_buffer.main_thread_flushCond);
			pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex);
		}
		methy_input_buffer.all_completed++;
		pthread_mutex_unlock(&methy_input_buffer.all_completed_Mutex);

		
		

		pthread_mutex_lock(&methy_input_buffer.process_Mutex[thread_id]);
		pthread_cond_wait(&methy_input_buffer.process_stallCond[thread_id], &methy_input_buffer.process_Mutex[thread_id]);
		pthread_mutex_unlock(&methy_input_buffer.process_Mutex[thread_id]);
	}

}



void* process_single_methy_alignment_multiple_threads(void* arg)
{
	int thread_id = *((int *)arg);
	int return_value;
	buffer_3_bit read;
	Init_buffer_3_bit(read);
	int read_length;
	int C_or_G;
	int direction;
	bitmapper_bs_iter pos;




	while (methy_input_buffer.all_end == 0)
	{

		///1是拿到一个有效的alignment
		///0是读到文件末尾了
		///-1是这个alignment重复了, 拿到了一个空的alignment
		while (return_value = get_alignment_single_multiple_threads(
			parameters.read_file, parameters.ref_genome, parameters.start_pos,
			&read, &read_length, &C_or_G, &direction, &pos, thread_id))
		{

			parameters.total_read[thread_id]++;

			if (return_value == -1)
			{
				parameters.duplicate_read[thread_id]++;
				continue;
			}
			else if (return_value == 1)
			{
				parameters.unique_read[thread_id]++;
			}



			if (pos > methy_input_buffer.each_buffer_interval_end[thread_id] || pos < methy_input_buffer.each_buffer_interval_start[thread_id])
			{
				fprintf(stderr, "ERROR: process_PE_methy_alignment_multiple_threads\n");
				fprintf(stderr, "pos: %llu\n", pos);
				fprintf(stderr, "start: %llu\n", methy_input_buffer.each_buffer_interval_start[thread_id]);
				fprintf(stderr, "end: %llu\n", methy_input_buffer.each_buffer_interval_end[thread_id]);
				exit(0);
			}


			if (pos + read_length - 1 <= methy_input_buffer.each_buffer_interval_end[thread_id])
			{
				update_methylation(C_or_G, parameters.ref_genome + pos - parameters.start_pos, read.buffer, read_length,
					parameters.methy + pos - parameters.start_pos, parameters.total + pos - parameters.start_pos,
					parameters.correct_methy + pos - parameters.start_pos, parameters.correct_total + pos - parameters.start_pos);
			}
			else
			{
				pthread_mutex_lock(&methy_input_buffer.overlap_file_Mutex);

				write_overlap_alignment(read.buffer, read_length, C_or_G, direction, pos,
					parameters.overlap_file);

				pthread_mutex_unlock(&methy_input_buffer.overlap_file_Mutex);
			}

		}

		pthread_mutex_lock(&methy_input_buffer.all_completed_Mutex);

		///这个是先判断，然后再+1，是因为除了这些处理线程之外，还有个input线程也要+1
		if (methy_input_buffer.all_completed == THREAD_COUNT)
		{

			pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex);
			pthread_cond_signal(&methy_input_buffer.main_thread_flushCond);
			pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex);
		}
		methy_input_buffer.all_completed++;
		pthread_mutex_unlock(&methy_input_buffer.all_completed_Mutex);




		pthread_mutex_lock(&methy_input_buffer.process_Mutex[thread_id]);
		pthread_cond_wait(&methy_input_buffer.process_stallCond[thread_id], &methy_input_buffer.process_Mutex[thread_id]);
		pthread_mutex_unlock(&methy_input_buffer.process_Mutex[thread_id]);
	}


}




void methy_extract_PE_mutiple_thread(int thread_id, char* file_name, int PE_distance, int* need_context)
{

	cut_length = _msf_refGenLength / genome_cuts;
	////if ((_msf_refGenLength * 2) % genome_cuts != 0)
	if (_msf_refGenLength % genome_cuts != 0)
	{
		cut_length = cut_length + 1;
	}


	int i = 0;
	char tmp_file_name[SEQ_MAX_LENGTH];
	char abnormal_file_name[SEQ_MAX_LENGTH];
	char overlap_file_name[SEQ_MAX_LENGTH];
	
	bitmapper_bs_iter L_read1[SEQ_MAX_LENGTH];
	bitmapper_bs_iter L_read2[SEQ_MAX_LENGTH];
	

	bitmapper_bs_iter flag;
	bitmapper_bs_iter pos1;
	bitmapper_bs_iter pos2;
	int C_or_G1;
	int C_or_G2;
	int direction1;
	int direction2;
	bitmapper_bs_iter inner_i;
	int return_value;

	bitmapper_bs_iter inner_start_pos;
	int read1_length;
	int read2_length;
	char last_two_bases[2];




	long long total_read = 0;
	long long duplicate_read = 0;
	long long unique_read = 0;
	long long abnormal_read = 0;
	long long overlap_read = 0;

	sprintf(abnormal_file_name, "%s_abnormal%d", file_name, thread_id);
	sprintf(overlap_file_name, "%s_overlap%d", file_name, thread_id);

	

	bitmapper_bs_iter tmp_genome_cuts;
	int tmp_maxDistance_pair;
	int tmp_mode;


	///bitmapper_bs_iter tmp_buffer[100];

	double load_time = 0;
	double start_load_time = 0;
	double print_time = 0;
	double start_print_time = 0;
	double update_methy_time = 0;
	double start_update_methy_time = 0;

	pthread_t inputReadsHandle;

	pthread_t *_r_threads;


	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*THREAD_COUNT);


	pthread_t *_r_threads_output;

	_r_threads_output = (pthread_t *)malloc(sizeof(pthread_t)*THREAD_COUNT);

	
	init_PE_methy_parameters(&parameters, PE_distance, cut_length);

	int j;



	for (i = 0; i < genome_cuts; i++)
	{
		/**
		fprintf(stderr, "i: %d\n", i);
		**/

		sprintf(tmp_file_name, "%s_%d.bmm", file_name, i);
		parameters.read_file = fopen(tmp_file_name, "r");

		///这里要写异常文件
		parameters.abnormal_file = fopen(abnormal_file_name, "wb");
		///这里要写多线程之间的重叠文件
		parameters.overlap_file = fopen(overlap_file_name, "wb");

		fscanf(parameters.read_file, "%llu\t%d\t%d\t%llu\t%llu\n", &tmp_genome_cuts, 
			&tmp_maxDistance_pair, &tmp_mode, &parameters.start_pos, &parameters.end_pos);

		parameters.length = parameters.end_pos - parameters.start_pos + 1;

		
		
		


		start_load_time = Get_T();

		if (!init_genome_cut_PE(parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total,
			parameters.ref_genome, parameters.length, i, parameters.extra_length, last_two_bases, &parameters.PE_de))
		{
			fprintf(stderr, "ERROR when reading genome!\n");
			exit(0);
		}

		load_time += Get_T() - start_load_time;


		start_update_methy_time = Get_T();

		init_methy_input_buffer(methylation_size/2, THREAD_COUNT,
			parameters.start_pos, parameters.end_pos, i, 1);



		input_methy_alignment_file = parameters.read_file;

		if (i == 0)
		{
			pthread_create(&inputReadsHandle, NULL, input_methy_muti_threads, NULL);

			for (j = 0; j < THREAD_COUNT; j++)
			{
				
				int *arg = (int*)malloc(sizeof(*arg));
				*arg = j;

				pthread_create(_r_threads + j, NULL, process_PE_methy_alignment_multiple_threads, (void*)arg);
			
			}

		}
		else
		{
			pthread_mutex_lock(&methy_input_buffer.input_thread_Mutex);
			pthread_cond_signal(&methy_input_buffer.input_thread_flushCond);
			pthread_mutex_unlock(&methy_input_buffer.input_thread_Mutex);

			for (j = 0; j < THREAD_COUNT; j++)
			{
				pthread_mutex_lock(&methy_input_buffer.process_Mutex[j]);
				pthread_cond_signal(&methy_input_buffer.process_stallCond[j]);
				pthread_mutex_unlock(&methy_input_buffer.process_Mutex[j]);
			}

		}

		

		pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex);
		pthread_cond_wait(&methy_input_buffer.main_thread_flushCond, &methy_input_buffer.main_thread_Mutex);
		pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex);
		
		update_methy_time += Get_T() - start_update_methy_time;
		

		////process_PE_methy_alignment(thread_id);


		///要先处理overlap文件，再处理异常文件
		///这里要关闭overlap文件, 再打开
		///请注意，这种overlap的alignment已经去过重了
		fclose(parameters.overlap_file);
		parameters.overlap_file = fopen(overlap_file_name, "rb");
		phrase_overlap_alignment(parameters.ref_genome, parameters.start_pos, L_read1, parameters.overlap_file,
			parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total, &overlap_read);
		fclose(parameters.overlap_file);
		

		///这里要关闭异常文件, 再打开
		fclose(parameters.abnormal_file);
		parameters.abnormal_file = fopen(abnormal_file_name, "rb");
		phrase_abnormal_alignment(parameters.ref_genome, parameters.start_pos, L_read1, L_read2, parameters.abnormal_file, &parameters.PE_de, 
			parameters.PE_distance,
			parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total, &unique_read, &duplicate_read);
		fclose(parameters.abnormal_file);

		
		
		


		start_print_time = Get_T();

		
		print_methylation_multiple_thread(parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total,
			parameters.start_pos, parameters.end_pos, parameters.ref_genome, last_two_bases, need_context,
			_r_threads_output, i);
		

		/**
		print_methylation(parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total,
			parameters.start_pos, parameters.end_pos, parameters.ref_genome, last_two_bases, need_context);
		**/

		print_time += Get_T() - start_print_time;


		fclose(parameters.read_file);


		last_two_bases[0] = parameters.ref_genome[parameters.length - 2];
		last_two_bases[1] = parameters.ref_genome[parameters.length - 1];

		///break;
	}


	/**
	fprintf(stderr, "******\n");
	**/

	methy_input_buffer.all_end = 1;

	pthread_mutex_lock(&methy_input_buffer.input_thread_Mutex);
	pthread_cond_signal(&methy_input_buffer.input_thread_flushCond);
	pthread_mutex_unlock(&methy_input_buffer.input_thread_Mutex);

	for (j = 0; j < THREAD_COUNT; j++)
	{
		pthread_mutex_lock(&methy_input_buffer.process_Mutex[j]);
		pthread_cond_signal(&methy_input_buffer.process_stallCond[j]);
		pthread_mutex_unlock(&methy_input_buffer.process_Mutex[j]);

		pthread_mutex_lock(&methy_input_buffer.process_Mutex_methy[j]);
		pthread_cond_signal(&methy_input_buffer.process_stallCond_methy[j]);
		pthread_mutex_unlock(&methy_input_buffer.process_Mutex_methy[j]);

	}
	/**
	fprintf(stderr, "1###########\n");
	**/


	pthread_join(inputReadsHandle, NULL);

	/**
	fprintf(stderr, "2###########\n");
	**/

	for (j = 0; j < THREAD_COUNT; j++)
	{
		///fprintf(stderr, "-j: %d\n", j);
		pthread_join(_r_threads[j], NULL);
		///这个join不注释掉程序无法退出。。
		//pthread_join(_r_threads_output[j], NULL);
	}

	/**
	for (j = 0; j < THREAD_COUNT; j++)
	{
		fprintf(stderr, "+j: %d\n", j);
		pthread_join(_r_threads_output[j], NULL);
	}
	**/


	

	///fprintf(stderr, "3###########\n");


	///至少要把这个异常文件删除了
	remove(abnormal_file_name);
	remove(overlap_file_name);

	for (i = 0; i < THREAD_COUNT; i++)
	{
		total_read += parameters.total_read[i];
		duplicate_read += parameters.duplicate_read[i];
		unique_read += parameters.unique_read[i];
		abnormal_read += parameters.abnormal_read[i];
	}


	fprintf(stderr, "total_read: %lld\n", total_read);
	fprintf(stderr, "duplicate_read: %lld\n", duplicate_read);
	fprintf(stderr, "unique_read: %lld\n", unique_read);
	fprintf(stderr, "abnormal_read: %lld\n", abnormal_read);
	fprintf(stderr, "overlap_read: %lld\n", overlap_read);

	

	fprintf(stderr, "load_time: %f\n", load_time);
	fprintf(stderr, "print_time: %f\n", print_time);
	fprintf(stderr, "update_methy_time: %f\n", update_methy_time);

	



}









void methy_extract_PE(int thread_id, char* file_name, int PE_distance, int* need_context)
{

	cut_length = _msf_refGenLength / genome_cuts;
	////if ((_msf_refGenLength * 2) % genome_cuts != 0)
	if (_msf_refGenLength % genome_cuts != 0)
	{
		cut_length = cut_length + 1;
	}


	FILE* read_file;
	FILE* abnormal_file;
	int i = 0;
	char tmp_file_name[SEQ_MAX_LENGTH];
	char abnormal_file_name[SEQ_MAX_LENGTH];
	bitmapper_bs_iter read1[SEQ_MAX_LENGTH];
	bitmapper_bs_iter read2[SEQ_MAX_LENGTH];
	bitmapper_bs_iter start_pos, end_pos;
	bitmapper_bs_iter length;
	bitmapper_bs_iter extra_length = SEQ_MAX_LENGTH * 2 + PE_distance; ///这个是为了处理边界情况


	uint16_t* methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	uint16_t* total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	uint16_t* correct_methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	uint16_t* correct_total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	char* ref_genome = (char*)malloc(sizeof(char)*(cut_length + extra_length));

	deduplicate_PE PE_de;
	///按道理来说, PE_de的长度不该有这个extra_length, 但是多加一点以免出错吧
	init_deduplicate_PE(&PE_de, PE_distance, cut_length + extra_length);

	bitmapper_bs_iter flag;
	bitmapper_bs_iter pos1;
	bitmapper_bs_iter pos2;
	int C_or_G1;
	int C_or_G2;
	int direction1;
	int direction2;
	bitmapper_bs_iter inner_i;
	int return_value;

	bitmapper_bs_iter inner_start_pos;
	int read1_length;
	int read2_length;
	char last_two_bases[2];

	


	long long total_read = 0;
	long long duplicate_read = 0;
	long long unique_read = 0;
	long long abnormal_read = 0;

	sprintf(abnormal_file_name, "%s_abnormal%d", file_name, thread_id);

	bitmapper_bs_iter tmp_genome_cuts;
	int tmp_maxDistance_pair;
	int tmp_mode;


	///bitmapper_bs_iter tmp_buffer[100];

	double load_time = 0;
	double start_load_time = 0;
	double update_methylation_time = 0;
	double start_update_methylation_time = 0;
	double print_time = 0;
	double start_print_time = 0;



	for (i = 0; i < genome_cuts; i++)
	{
		sprintf(tmp_file_name, "%s_%d.bmm", file_name, i);
		read_file = fopen(tmp_file_name, "r");

		///这里要写异常文件
		abnormal_file = fopen(abnormal_file_name, "wb");

		fscanf(read_file, "%llu\t%d\t%d\t%llu\t%llu\n", &tmp_genome_cuts, &tmp_maxDistance_pair, &tmp_mode, &start_pos, &end_pos);

		length = end_pos - start_pos + 1;



		start_load_time = Get_T();

		if (!init_genome_cut_PE(methy, total, correct_methy, correct_total, 
			ref_genome, length, i, extra_length, last_two_bases, &PE_de))
		{
			fprintf(stderr, "ERROR when reading genome!\n");
			exit(0);
		}
		load_time += Get_T() - start_load_time;

		///1是拿到一个有效的alignment
		///0是读到文件末尾了
		///-1是这个alignment重复了, 拿到了一个空的alignment
		////2是异常情况, 就是读了个5'-end位置有问题的结果
		while (return_value = get_alignment_PE(read_file, ref_genome, start_pos, 
			read1, &read1_length, &C_or_G1, &direction1, &pos1,
			read2, &read2_length, &C_or_G2, &direction2, &pos2,
			&PE_de, PE_distance))
		{



			total_read++;

			if (return_value == -1)
			{
				duplicate_read++;
				continue;
			}
			else if (return_value == 1)
			{
				unique_read++;
			}
			else   ///这是返回值为2的情况，也是5'端位置异常的情况  ///暂时还没处理
			{
				abnormal_read++;  


				///将异常结果写入临时文件
				write_abnormal_PE_alignment(
					read1, read1_length, C_or_G1, direction1, pos1,
					read2, read2_length, C_or_G2, direction2, pos2,
					abnormal_file);

				continue;

			}

			start_update_methylation_time = Get_T();
			
			update_methylation(C_or_G1, ref_genome + pos1 - start_pos, read1, read1_length,
				methy + pos1 - start_pos, total + pos1 - start_pos,
				correct_methy + pos1 - start_pos, correct_total + pos1 - start_pos);

			update_methylation(C_or_G2, ref_genome + pos2 - start_pos, read2, read2_length,
				methy + pos2 - start_pos, total + pos2 - start_pos,
				correct_methy + pos2 - start_pos, correct_total + pos2 - start_pos);
			
			update_methylation_time += Get_T() - start_update_methylation_time;


		}
		


		///这里要关闭异常文件, 再打开
		fclose(abnormal_file);
		abnormal_file = fopen(abnormal_file_name, "rb");
		
		phrase_abnormal_alignment(ref_genome, start_pos, read1, read2, abnormal_file, &PE_de, PE_distance,
			methy, total, correct_methy, correct_total, &unique_read, &duplicate_read);
		

		fclose(abnormal_file);
		
		start_print_time = Get_T();

		
		print_methylation(methy, total, correct_methy, correct_total,
			start_pos, end_pos, ref_genome, last_two_bases, need_context);
		
		print_time += Get_T() - start_print_time;


		fclose(read_file);


		last_two_bases[0] = ref_genome[length - 2];
		last_two_bases[1] = ref_genome[length - 1];


	}

	///至少要把这个异常文件删除了
	remove(abnormal_file_name);

	fprintf(stderr, "total_read: %lld\n", total_read);
	fprintf(stderr, "duplicate_read: %lld\n", duplicate_read);
	fprintf(stderr, "unique_read: %lld\n", unique_read);
	fprintf(stderr, "abnormal_read: %lld\n", abnormal_read);

	fprintf(stderr, "load_time: %f\n", load_time);
	fprintf(stderr, "update_methylation_time: %f\n", update_methylation_time);
	fprintf(stderr, "print_time: %f\n", print_time);
	
	
	


}







void methy_extract_mutiple_thread(int thread_id, char* file_name, int* need_context)
{

	cut_length = _msf_refGenLength / genome_cuts;
	////if ((_msf_refGenLength * 2) % genome_cuts != 0)
	if (_msf_refGenLength % genome_cuts != 0)
	{
		cut_length = cut_length + 1;
	}


	int i = 0;
	char tmp_file_name[SEQ_MAX_LENGTH];
	char overlap_file_name[SEQ_MAX_LENGTH];

	bitmapper_bs_iter L_read1[SEQ_MAX_LENGTH];


	bitmapper_bs_iter flag;
	bitmapper_bs_iter pos;
	int C_or_G;
	int direction;
	bitmapper_bs_iter inner_i;
	int return_value;

	bitmapper_bs_iter inner_start_pos;
	int read_length;
	char last_two_bases[2];




	long long total_read = 0;
	long long duplicate_read = 0;
	long long unique_read = 0;
	long long abnormal_read = 0;
	long long overlap_read = 0;

	sprintf(overlap_file_name, "%s_overlap%d", file_name, thread_id);



	bitmapper_bs_iter tmp_genome_cuts;
	int tmp_maxDistance_pair;
	int tmp_mode;


	///bitmapper_bs_iter tmp_buffer[100];

	double load_time = 0;
	double start_load_time = 0;
	double print_time = 0;
	double start_print_time = 0;
	double update_methy_time = 0;
	double start_update_methy_time = 0;

	pthread_t inputReadsHandle;

	pthread_t *_r_threads;


	_r_threads = (pthread_t *)malloc(sizeof(pthread_t)*THREAD_COUNT);


	pthread_t *_r_threads_output;

	_r_threads_output = (pthread_t *)malloc(sizeof(pthread_t)*THREAD_COUNT);


	init_single_methy_parameters(&parameters, cut_length);

	int j;



	for (i = 0; i < genome_cuts; i++)
	{

		sprintf(tmp_file_name, "%s_%d.bmm", file_name, i);
		parameters.read_file = fopen(tmp_file_name, "r");


		///这里要写多线程之间的重叠文件
		parameters.overlap_file = fopen(overlap_file_name, "wb");

		fscanf(parameters.read_file, "%llu\t%d\t%d\t%llu\t%llu\n", &tmp_genome_cuts,
			&tmp_maxDistance_pair, &tmp_mode, &parameters.start_pos, &parameters.end_pos);

		parameters.length = parameters.end_pos - parameters.start_pos + 1;




		start_load_time = Get_T();

		if (!init_genome_cut(parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total, 
			parameters.ref_genome, parameters.length, i, parameters.extra_length, last_two_bases))
		{
			fprintf(stderr, "ERROR when reading genome!\n");
			exit(0);
		}


		load_time += Get_T() - start_load_time;


		start_update_methy_time = Get_T();

		init_methy_input_buffer(methylation_size, THREAD_COUNT,
			parameters.start_pos, parameters.end_pos, i, 0);



		input_methy_alignment_file = parameters.read_file;

		if (i == 0)
		{
			pthread_create(&inputReadsHandle, NULL, input_methy_muti_threads_single_end, NULL);

			for (j = 0; j < THREAD_COUNT; j++)
			{

				int *arg = (int*)malloc(sizeof(*arg));
				*arg = j;

				pthread_create(_r_threads + j, NULL, process_single_methy_alignment_multiple_threads, (void*)arg);

			}

		}
		else
		{
			pthread_mutex_lock(&methy_input_buffer.input_thread_Mutex);
			pthread_cond_signal(&methy_input_buffer.input_thread_flushCond);
			pthread_mutex_unlock(&methy_input_buffer.input_thread_Mutex);

			for (j = 0; j < THREAD_COUNT; j++)
			{
				pthread_mutex_lock(&methy_input_buffer.process_Mutex[j]);
				pthread_cond_signal(&methy_input_buffer.process_stallCond[j]);
				pthread_mutex_unlock(&methy_input_buffer.process_Mutex[j]);
			}

		}



		pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex);
		pthread_cond_wait(&methy_input_buffer.main_thread_flushCond, &methy_input_buffer.main_thread_Mutex);
		pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex);

		update_methy_time += Get_T() - start_update_methy_time;


		////process_PE_methy_alignment(thread_id);


		///要先处理overlap文件，再处理异常文件
		///这里要关闭overlap文件, 再打开
		///请注意，这种overlap的alignment已经去过重了
		fclose(parameters.overlap_file);
		parameters.overlap_file = fopen(overlap_file_name, "rb");
		phrase_overlap_alignment(parameters.ref_genome, parameters.start_pos, L_read1, parameters.overlap_file,
			parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total, &overlap_read);
		fclose(parameters.overlap_file);









		start_print_time = Get_T();


		print_methylation_multiple_thread(parameters.methy, parameters.total, parameters.correct_methy, parameters.correct_total,
			parameters.start_pos, parameters.end_pos, parameters.ref_genome, last_two_bases, need_context,
			_r_threads_output, i);


		print_time += Get_T() - start_print_time;


		fclose(parameters.read_file);


		last_two_bases[0] = parameters.ref_genome[parameters.length - 2];
		last_two_bases[1] = parameters.ref_genome[parameters.length - 1];

		///break;
	}




	methy_input_buffer.all_end = 1;

	pthread_mutex_lock(&methy_input_buffer.input_thread_Mutex);
	pthread_cond_signal(&methy_input_buffer.input_thread_flushCond);
	pthread_mutex_unlock(&methy_input_buffer.input_thread_Mutex);

	for (j = 0; j < THREAD_COUNT; j++)
	{
		pthread_mutex_lock(&methy_input_buffer.process_Mutex[j]);
		pthread_cond_signal(&methy_input_buffer.process_stallCond[j]);
		pthread_mutex_unlock(&methy_input_buffer.process_Mutex[j]);

		pthread_mutex_lock(&methy_input_buffer.process_Mutex_methy[j]);
		pthread_cond_signal(&methy_input_buffer.process_stallCond_methy[j]);
		pthread_mutex_unlock(&methy_input_buffer.process_Mutex_methy[j]);

	}



	pthread_join(inputReadsHandle, NULL);



	for (j = 0; j < THREAD_COUNT; j++)
	{
		///fprintf(stderr, "-j: %d\n", j);
		pthread_join(_r_threads[j], NULL);
		///这个join不注释掉程序无法退出。。
		//pthread_join(_r_threads_output[j], NULL);
	}






	///fprintf(stderr, "3###########\n");



	remove(overlap_file_name);

	for (i = 0; i < THREAD_COUNT; i++)
	{
		total_read += parameters.total_read[i];
		duplicate_read += parameters.duplicate_read[i];
		unique_read += parameters.unique_read[i];
		////abnormal_read += parameters.abnormal_read[i];
	}


	fprintf(stderr, "total_read: %lld\n", total_read);
	fprintf(stderr, "duplicate_read: %lld\n", duplicate_read);
	fprintf(stderr, "unique_read: %lld\n", unique_read);
	////fprintf(stderr, "abnormal_read: %lld\n", abnormal_read);
	fprintf(stderr, "overlap_read: %lld\n", overlap_read);



	fprintf(stderr, "load_time: %f\n", load_time);
	fprintf(stderr, "print_time: %f\n", print_time);
	fprintf(stderr, "update_methy_time: %f\n", update_methy_time);





}






void methy_extract(int thread_id, char* file_name, int* need_context)
{

	cut_length = _msf_refGenLength / genome_cuts;
	////if ((_msf_refGenLength * 2) % genome_cuts != 0)
	if (_msf_refGenLength % genome_cuts != 0)
	{
		cut_length = cut_length + 1;
	}

	long long total_memory_size = 0;

	FILE* read_file;
	int i = 0;
	char tmp_file_name[SEQ_MAX_LENGTH];
	bitmapper_bs_iter read[SEQ_MAX_LENGTH];
	bitmapper_bs_iter start_pos, end_pos;
	bitmapper_bs_iter length;
	bitmapper_bs_iter extra_length = SEQ_MAX_LENGTH * 2; ///这个是为了处理边界情况
	///bitmapper_bs_iter extra_length = SEQ_MAX_LENGTH * 2 + maxDistance_pair; ///这个是为了处理边界情况

	///bitmapper_bs_iter tmp_buffer[100];

	uint16_t* methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	total_memory_size = total_memory_size + sizeof(uint16_t)*(cut_length + extra_length);

	uint16_t* total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	total_memory_size = total_memory_size + sizeof(uint16_t)*(cut_length + extra_length);

	uint16_t* correct_methy = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	total_memory_size = total_memory_size + sizeof(uint16_t)*(cut_length + extra_length);

	uint16_t* correct_total = (uint16_t*)malloc(sizeof(uint16_t)*(cut_length + extra_length));
	total_memory_size = total_memory_size + sizeof(uint16_t)*(cut_length + extra_length);

	char* ref_genome = (char*)malloc(sizeof(char)*(cut_length + extra_length));
	total_memory_size = total_memory_size + sizeof(char)*(cut_length + extra_length);


	bitmapper_bs_iter flag;
	bitmapper_bs_iter pos;
	int C_or_G;
	int direction;
	bitmapper_bs_iter inner_i;
	int return_value;

	bitmapper_bs_iter inner_start_pos;
	int read_length;
	char last_two_bases[2];

	bitmapper_bs_iter tmp_genome_cuts;
	int tmp_maxDistance_pair;
	int tmp_mode;

	long long total_read = 0;
	long long duplicate_read = 0;
	long long unique_read = 0;

	double tmp_pure_print_time;

	for (i = 0; i < genome_cuts; i++)
	{
		sprintf(tmp_file_name, "%s_%d.bmm", file_name, i);
		read_file = fopen(tmp_file_name, "r");

		fscanf(read_file, "%llu\t%d\t%d\t%llu\t%llu\n", &tmp_genome_cuts, &tmp_maxDistance_pair, &tmp_mode, &start_pos, &end_pos);

		length = end_pos - start_pos + 1;

		if (!init_genome_cut(methy, total, correct_methy, correct_total, ref_genome, length, i, extra_length, last_two_bases))
		{
			fprintf(stderr, "ERROR when reading genome!\n");
			exit(0);
		}


		///fprintf(stderr, "ref_genome[977002]: %llu\n", ref_genome[977002]);

		///return_value=0则读到文件末尾
		///return_value = -1读到了一个重复的alignment
		while (return_value = get_alignment(read_file, ref_genome, start_pos, read,
			&read_length, &C_or_G, &direction, &pos))
		{

			total_read++;

			if (return_value == -1)
			{
				duplicate_read++;

				continue;
			}
			else
			{
				unique_read++;
			}


			///fprintf(stderr, "C_or_G: %d, pos: %llu, read:%s\n", C_or_G, pos, read);



			///debug_pos = pos;





			update_methylation(C_or_G, ref_genome + pos - start_pos, (bitmapper_bs_iter*)read, read_length,
				methy + pos - start_pos, total + pos - start_pos,
				correct_methy + pos - start_pos, correct_total + pos - start_pos);




		}

		print_methylation(methy, total, correct_methy, correct_total,
			start_pos, end_pos, ref_genome, last_two_bases, need_context);


		fclose(read_file);


		last_two_bases[0] = ref_genome[length - 2];
		last_two_bases[1] = ref_genome[length - 1];




		//break;

	}


	///fprintf(stderr, "total_memory_size: %lld\n", total_memory_size);

	fprintf(stderr, "total_read: %lld\n", total_read);
	fprintf(stderr, "duplicate_read: %lld\n", duplicate_read);
	fprintf(stderr, "unique_read: %lld\n", unique_read);



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







inline void map_candidate_votes_mutiple_cut_end_to_end_4_sse(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m128i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
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



		BS_Reserve_Banded_BPM_4_SSE_only(t[0], t[1], t[2], t[3], p_length, read, t_length,
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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

			(*min_err) = return_sites_error[0];
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}


		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[1] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}



		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[2] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[2];

			(*min_err) = return_sites_error[2];
			(*min_err_index) = i + 2;
			min_err_site = tmp_min_err_site;
		}



		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[3] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[3];

			(*min_err) = return_sites_error[3];
			(*min_err_index) = i + 3;
			min_err_site = tmp_min_err_site;
		}




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
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (candidate_votes[i].err == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_4_SSE_only(t[0], t[1], t[2], t[3], p_length, read, t_length,
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
				min_err_site != tmp_min_err_site
				&&
				(*min_err_index) >= 0)
			{
				(*second_best_diff) = 0;

				(*min_err_index) = -2 - (*min_err_index);
				///(*min_err_index) = -2;
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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

	/**
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/
	


}

inline void map_candidate_votes_mutiple_cut_end_to_end_2_sse(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m128i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[2];
	unsigned int return_sites_error[2];

	char* t[2];


	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;


	i = 0;

	///这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	(*min_err_index) = -1;
	bitmapper_bs_iter min_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;

	while (i + 2 <= *candidate_votes_length)
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






		BS_Reserve_Banded_BPM_2_SSE_only(t[0], t[1], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i + 1].end_site = return_sites[1];


		candidate_votes[i].err = return_sites_error[0];
		candidate_votes[i + 1].err = return_sites_error[1];



		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[0] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

			(*min_err) = return_sites_error[0];
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}


		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[1] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}


		i = i + 2;

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_2_SSE_only(t[0], t[1], p_length, read, t_length,
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
				min_err_site != tmp_min_err_site
				&&
				(*min_err_index) >= 0)
			{
				(*second_best_diff) = 0;

				(*min_err_index) = -2 - (*min_err_index);
				///(*min_err_index) = -2;
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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

	/**
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/


}


inline void map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m128i* Peq_SSE, unsigned int* min_err, int* best_mapping_occ,
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





		BS_Reserve_Banded_BPM_4_SSE_only(t[0], t[1], t[2], t[3], p_length, read, t_length,
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


		BS_Reserve_Banded_BPM_4_SSE_only(t[0], t[1], t[2], t[3], p_length, read, t_length,
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





inline void map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m128i* Peq_SSE, unsigned int* min_err, int* best_mapping_occ,
	char* name)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[2];
	unsigned int return_sites_error[2];

	char* t[2];



	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;



	i = 0;

	///这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	///(*min_err_index) = -1;
	(*best_mapping_occ) = 0;
	bitmapper_bs_iter pre_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;



	while (i + 2 <= *candidate_votes_length)
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








		BS_Reserve_Banded_BPM_2_SSE_only(t[0], t[1], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);



		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i + 1].end_site = return_sites[1];


		candidate_votes[i].err = return_sites_error[0];
		candidate_votes[i + 1].err = return_sites_error[1];







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

		i = i + 2;

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


		BS_Reserve_Banded_BPM_2_SSE_only(t[0], t[1], p_length, read, t_length,
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

inline int map_candidate_votes_mutiple_end_to_end_4_sse(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m128i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
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



		BS_Reserve_Banded_BPM_4_SSE_only(t[0], t[1], t[2], t[3], p_length, read, t_length,
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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[2];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[3];

			(*min_err) = return_sites_error[3];
			(*min_err_index) = i + 3;
			min_err_site = tmp_min_err_site;
		}

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
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (candidate_votes[i].err == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_4_SSE_only(t[0], t[1], t[2], t[3], p_length, read, t_length,
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
				(*second_best_diff) = 0;

				if ((*min_err_index) >= 0)
				{
					(*min_err_index) = -2 - (*min_err_index);
				}

				if ((*min_err) == 0)
				{
					return 0;
				}
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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

	/**
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/


	return 1;

}




inline int map_candidate_votes_mutiple_end_to_end_2_sse(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m128i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	int return_sites[2];
	unsigned int return_sites_error[2];
	int mumber_of_exact_matches = 0;

	char* t[2];

	t[0] = tmp_ref;
	t[1] = t[0] + p_length + 32;


	i = 0;

	//这个数不能随便改...
	(*min_err) = ((unsigned int)-1) - 1;
	(*min_err_index) = -1;

	bitmapper_bs_iter min_err_site = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter tmp_min_err_site;

	while (i + 2 <= *candidate_votes_length)
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




		BS_Reserve_Banded_BPM_2_SSE_only(t[0], t[1], p_length, read, t_length,
			return_sites, return_sites_error, error_threshold, Peq_SSE);


		candidate_votes[i].end_site = return_sites[0];
		candidate_votes[i + 1].end_site = return_sites[1];


		candidate_votes[i].err = return_sites_error[0];
		candidate_votes[i + 1].err = return_sites_error[1];




		tmp_min_err_site = candidate_votes[i].site + candidate_votes[i].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[0] == (*min_err)
			&&
			min_err_site != tmp_min_err_site)
		{
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}

		i = i + 2;

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

			(*min_err) = candidate_votes[i].err;
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}

	}
	else if (last_length != 0)
	{


		BS_Reserve_Banded_BPM_2_SSE_only(t[0], t[1], p_length, read, t_length,
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
				(*second_best_diff) = 0;

				if ((*min_err_index) >= 0)
				{
					(*min_err_index) = -2 - (*min_err_index);
				}

				if ((*min_err) == 0)
				{
					return 0;
				}
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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
	/**
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/


	return 1;

}

















#if defined __AVX2__

inline void map_candidate_votes_mutiple_cut_end_to_end_4(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

			(*min_err) = return_sites_error[0];
			(*min_err_index) = i;
			min_err_site = tmp_min_err_site;
		}


		tmp_min_err_site = candidate_votes[i + 1].site + candidate_votes[i + 1].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[1] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

			(*min_err) = return_sites_error[1];
			(*min_err_index) = i + 1;
			min_err_site = tmp_min_err_site;
		}



		tmp_min_err_site = candidate_votes[i + 2].site + candidate_votes[i + 2].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[2] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[2];

			(*min_err) = return_sites_error[2];
			(*min_err_index) = i + 2;
			min_err_site = tmp_min_err_site;
		}



		tmp_min_err_site = candidate_votes[i + 3].site + candidate_votes[i + 3].end_site;
		///首先要最小值相同，其次位置要不同，这样的话才算有best map有多个
		if (return_sites_error[3] == (*min_err)
			&&
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[3];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

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
				min_err_site != tmp_min_err_site
				&&
				(*min_err_index) >= 0)
			{
				(*second_best_diff) = 0;

				(*min_err_index) = -2 - (*min_err_index);
				///(*min_err_index) = -2;
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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

	/**
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/


}

#endif



#if defined __AVX2__

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
#endif

#if defined __AVX2__

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

#endif




#if defined __AVX2__

inline void map_candidate_votes_mutiple_cut_end_to_end_8(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			//(*min_err_index) = -2;
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[2];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[3];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[4]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[4];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[5]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[5];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			//(*min_err_index) = -2;
		}
		else if (return_sites_error[6]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[6];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (return_sites_error[7]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[7];

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
			min_err_site != tmp_min_err_site
			&&
			(*min_err_index) >= 0)
		{
			(*second_best_diff) = 0;

			(*min_err_index) = -2 - (*min_err_index);
			///(*min_err_index) = -2;
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

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
				min_err_site != tmp_min_err_site
				&&
				(*min_err_index) >= 0)
			{
				(*second_best_diff) = 0;

				(*min_err_index) = -2 - (*min_err_index);
				///(*min_err_index) = -2;
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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

	/**
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/



}

#endif














#if defined __AVX2__

inline int map_candidate_votes_mutiple_end_to_end_8(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}


			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

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

			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[2];

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
			(*second_best_diff) = 0;


			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[3];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[4]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[4];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[5]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[5];


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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[6]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[6];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (return_sites_error[7]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[7];

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
			(*second_best_diff) = 0;
			
			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

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
				(*second_best_diff) = 0;

				if ((*min_err_index) >= 0)
				{
					(*min_err_index) = -2 - (*min_err_index);
				}

				if ((*min_err) == 0)
				{
					return 0;
				}
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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

	/**
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/





	return 1;

}



#endif




#if defined __AVX2__


inline int map_candidate_votes_mutiple_end_to_end_4(
	seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold,
	char* tmp_ref, char* read, __m256i* Peq_SSE, unsigned int* min_err, int* min_err_index, unsigned int* second_best_diff)
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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[0]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[0];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[1]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[1];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[2]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[2];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}

		}
		else if (return_sites_error[3]<(*min_err))
		{
			(*second_best_diff) = (*min_err) - return_sites_error[3];

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
			(*second_best_diff) = 0;

			if ((*min_err_index) >= 0)
			{
				(*min_err_index) = -2 - (*min_err_index);
			}

			if ((*min_err) == 0)
			{
				return 0;
			}
		}
		else if (candidate_votes[i].err < (*min_err))
		{
			(*second_best_diff) = (*min_err) - candidate_votes[i].err;

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
				(*second_best_diff) = 0;

				if ((*min_err_index) >= 0)
				{
					(*min_err_index) = -2 - (*min_err_index);
				}

				if ((*min_err) == 0)
				{
					return 0;
				}
			}
			else if (return_sites_error[inner_i]<(*min_err))
			{
				(*second_best_diff) = (*min_err) - return_sites_error[inner_i];

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
	if((*second_best_diff) > error_threshold)
	{
		(*second_best_diff) = - 2;
	}
	**/


	/**
	if (mumber_of_exact_matches >= 2)
	{
	(*min_err_index) = -2;
	return 0;
	}
	**/




	return 1;

}

#endif
















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




inline void directly_output_read1_return_buffer_unmapped_PE(char* name, char* read, char* r_read, char* qulity,
	int read_length, Output_buffer_sub_block* sub_block, bam_phrase* bam_groups, int flag)
{

	if (bam_output == 0)
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



		if (flag == 1)
		{
			sub_block->buffer[sub_block->length] = '7';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '7';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			sub_block->buffer[sub_block->length] = '1';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '4';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '1';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}


		///fprintf(schema_out_fp, "*\t0\t0\t*\t*\t0\t0\t");
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
		sub_block->buffer[sub_block->length] = '*';
		sub_block->length++;
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

		if (flag == 1)
		{
			output_to_buffer_char_length(sub_block, read, read_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			output_to_buffer_char_length(sub_block, r_read, read_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}


		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;


	}
	else
	{


		///这个一定一定要有
		///directly_output_read1_return_buffer里面一定要有置0, 但是在最后并不转换成bam
		///directly_output_read2_return_buffer这里不能置0，但是在最后需要转成bam
		///这样是为了避免在bam结果中，一对paired-end read的结果被分开了
		sub_block->length = 0;

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

		if (flag == 1)
		{
			sub_block->buffer[sub_block->length] = '7';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '7';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			sub_block->buffer[sub_block->length] = '1';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '4';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '1';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}


		///fprintf(schema_out_fp, "*\t0\t0\t*\t*\t0\t0\t");
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
		sub_block->buffer[sub_block->length] = '*';
		sub_block->length++;
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

		if (flag == 1)
		{
			output_to_buffer_char_length(sub_block, read, read_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			output_to_buffer_char_length(sub_block, r_read, read_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}

		///这里应该是'\0', 不是'\t'也不是'\n'
		output_to_buffer_char_length(sub_block, qulity, read_length);
		sub_block->buffer[sub_block->length] = '\0';

		convert_string_to_bam(sub_block->buffer, sub_block->length, bam_groups);


	}

}






inline void directly_output_read1_return_buffer(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length,
	Output_buffer_sub_block* sub_block, int paired_end_distance,
	bam_phrase* bam_groups, int output_mask, int mapq)
{

	if (bam_output == 0)
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


		result->flag = result->flag | output_mask;

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

		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);
	

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


			///output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

			///output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

				///output_to_buffer_int(sub_block, matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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


				///output_to_buffer_int(sub_block, another_matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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


	}
	else
	{


		///这个一定一定要有
		///directly_output_read1_return_buffer里面一定要有置0, 但是在最后并不转换成bam
		///directly_output_read2_return_buffer这里不能置0，但是在最后需要转成bam
		///这样是为了避免在bam结果中，一对paired-end read的结果被分开了
		sub_block->length = 0;



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


		result->flag = result->flag | output_mask;

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

		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);

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


			///output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

			///output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

				///output_to_buffer_int(sub_block, matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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


				///output_to_buffer_int(sub_block, another_matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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
		/**
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/


		sub_block->buffer[sub_block->length] = '\0';


		convert_string_to_bam(sub_block->buffer, sub_block->length, bam_groups);


	}

}





inline void output_paired_methy_buffer(char* name,
	char* read1, char* r_read1, map_result* result1, bitmapper_bs_iter r_length1,
	char* read2, char* r_read2, map_result* result2, bitmapper_bs_iter r_length2,
	Pair_Methylation* methy, pair_distance_count* distances)
{
	if ((*methy).current_size < (*methy).total_size)
	{

		long long end_5_pos1, end_5_pos2;

		/**read1的信息**/
		/*************这里要改**************/

		///注意r_length一定要+1, 这个是为了放/0的
		///Check_Space((*methy).R[0].reads[(*methy).current_size], (*methy).R[0].r_size[(*methy).current_size], r_length1 + 1);

		Check_Space_3_Bit((*methy).R[0].reads_3_bit[(*methy).current_size], 
			(*methy).R[0].r_size_3_bit[(*methy).current_size], r_length1)
		(*methy).R[0].r_real_length[(*methy).current_size] = r_length1;
		/*************这里要改**************/

		(*methy).R[0].r_length[(*methy).current_size] = r_length1;
		

		result1->site = result1->site + _ih_refGenName[result1->chrome_id].start_location - 1;



		if (result1->flag == 16)   ///read1比到反向互补链上了
		{
			result1->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ1;

			/*************这里要改**************/
			////correct_read_using_cigar(result1->cigar, r_read1, (*methy).R[0].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(result1->cigar, r_read1, (*methy).R[0].reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/

			end_5_pos1 = result1->site + r_length1;
		}
		else  ///正向链
		{

			result1->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ1;

			/*************这里要改**************/
			///correct_read_using_cigar(result1->cigar, read1, (*methy).R[0].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(result1->cigar, read1, (*methy).R[0].reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/


			end_5_pos1 = result1->site;
		}


		(*methy).R[0].r_length[(*methy).current_size] = result1->flag;

		(*methy).R[0].sites[(*methy).current_size] = result1->site;
		///(*methy).R[0].sites[(*methy).current_size] = end_5_pos1;


		/**read1的信息**/




		/**read2的信息**/

		/*************这里要改**************/
		///注意r_length一定要+1, 这个是为了放/0的
		///Check_Space((*methy).R[1].reads[(*methy).current_size], (*methy).R[1].r_size[(*methy).current_size], r_length2 + 1);
		Check_Space_3_Bit(
			(*methy).R[1].reads_3_bit[(*methy).current_size], 
			(*methy).R[1].r_size_3_bit[(*methy).current_size], 
			r_length2);
		(*methy).R[1].r_real_length[(*methy).current_size] = r_length2;
		/*************这里要改**************/

		(*methy).R[1].r_length[(*methy).current_size] = r_length2;


		result2->site = result2->site + _ih_refGenName[result2->chrome_id].start_location - 1;

		if (result2->flag == 16)   ///代表read2的反向互补链比对到了参考组的反向互补上, read2本身实际上比到了参考组正向链上
		{
			result2->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ2;

			/*************这里要改**************/
			///correct_read_using_cigar(result2->cigar, r_read2, (*methy).R[1].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit_over_lap_length_r2(
				result2->cigar, r_read2, (*methy).R[1].reads_3_bit[(*methy).current_size],
				result1->site, r_length1, result2->site, r_length2);
			/*************这里要改**************/

			
			end_5_pos2 = result2->site;
		}
		else
		{
			result2->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ2;

			/*************这里要改**************/
			///correct_read_using_cigar(result2->cigar, read2, (*methy).R[1].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit_over_lap_length_r2(
				result2->cigar, read2, (*methy).R[1].reads_3_bit[(*methy).current_size],
				result1->site, r_length1, result2->site, r_length2);
			/*************这里要改**************/


			end_5_pos2 = result2->site + r_length2;
		}


		/*************这里要改**************/
		///over_lap_length(result1->site, r_length1, result2->site, r_length2, (*methy).R[1].reads[(*methy).current_size]);
		/*************这里要改**************/


		(*methy).R[1].r_length[(*methy).current_size] = result2->flag;

		(*methy).R[1].sites[(*methy).current_size] = result2->site;
		///(*methy).R[1].sites[(*methy).current_size] = end_5_pos2;


		/**read2的信息**/

		int block_ID = get_cut_id_pair(methy, get_key_pair_methylation(methy, (*methy).current_size));

		///把小的那个当做key
		(*methy).cut_index[block_ID]++;

		(*methy).current_size++;




		long long end5_distance = end_5_pos2 - end_5_pos1;


		if (end5_distance < 0)
		{
			end5_distance = end5_distance * -1;
		}


		/**
		if (end5_distance != calculate_TLEN(result1->site, r_length1, result2->site, r_length2) && result1->site != result2->site)
		{
			fprintf(stderr, "end5_distance: %lld, TLEN: %lld\n", end5_distance, calculate_TLEN(result1->site, r_length1, result2->site, r_length2));
			fprintf(stderr, "site1: %lld, r_length1: %lld, flag1: %d\n", result1->site, r_length1, result1->flag);
			fprintf(stderr, "site2: %lld, r_length2: %lld, flag2: %d\n\n", result2->site, r_length2, result2->flag);
		}
		**/
	
		/**
		if (result1->flag == 99 && result1->site > result2->site)
		{
			fprintf(stderr, "+%s\n", name);
			fprintf(stderr, "site1: %lld, r_length1: %lld, flag1: %d\n", result1->site, r_length1, result1->flag);
			fprintf(stderr, "site2: %lld, r_length2: %lld, flag2: %d\n\n", result2->site, r_length2, result2->flag);
		}


		if (result1->flag == 83 && result1->site < result2->site)
		{
			fprintf(stderr, "-%s\n", name);
			fprintf(stderr, "site1: %lld, r_length1: %lld, flag1: %d\n", result1->site, r_length1, result1->flag);
			fprintf(stderr, "site2: %lld, r_length2: %lld, flag2: %d\n\n", result2->site, r_length2, result2->flag);
		}
		**/

		/**
		if (result1->flag == 99 && end_5_pos1 > end_5_pos2)
		{
			fprintf(stderr, "+%s\n", name);
			fprintf(stderr, "end_5_pos1: %lld, flag1: %d, end_5_pos2: %lld, flag2: %d\n\n", end_5_pos1, result1->flag, 
				end_5_pos2, result2->flag);
		}


		if (result1->flag == 83 && end_5_pos1 < end_5_pos2)
		{
			fprintf(stderr, "-%s\n", name);
			fprintf(stderr, "end_5_pos1: %lld, flag1: %d, end_5_pos2: %lld, flag2: %d\n\n", end_5_pos1, result1->flag,
				end_5_pos2, result2->flag);
		}
		**/


		///distances->count[block_ID][end5_distance - distances->minDistance_pair]++;
		distances->count[block_ID][end5_distance]++;


	}


	if ((*methy).current_size == (*methy).total_size)
	{
		assign_cuts_pair(methy);

		output_methylation_pair(methy);

		clear_methylation_pair(methy);
	}

}






inline void output_paired_methy_buffer_mutiple_threads(
	char* read1, char* r_read1, map_result* result1, bitmapper_bs_iter r_length1,
	char* read2, char* r_read2, map_result* result2, bitmapper_bs_iter r_length2,
	Pair_Methylation* methy, pair_distance_count* distances)
{
	if ((*methy).current_size < (*methy).total_size)
	{

		long long end_5_pos1, end_5_pos2;

		/**read1的信息**/
		///注意r_length一定要+1, 这个是为了放/0的
		/*************这里要改**************/
		///Check_Space((*methy).R[0].reads[(*methy).current_size], (*methy).R[0].r_size[(*methy).current_size], r_length1 + 1);
		Check_Space_3_Bit((*methy).R[0].reads_3_bit[(*methy).current_size],
			(*methy).R[0].r_size_3_bit[(*methy).current_size], r_length1)
			(*methy).R[0].r_real_length[(*methy).current_size] = r_length1;
		(*methy).R[0].r_real_length[(*methy).current_size] = r_length1;
		/*************这里要改**************/


		(*methy).R[0].r_length[(*methy).current_size] = r_length1;
		


		result1->site = result1->site + _ih_refGenName[result1->chrome_id].start_location - 1;



		if (result1->flag == 16)   ///read1比到反向互补链上了
		{
			result1->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ1;

			/*************这里要改**************/
			///correct_read_using_cigar(result1->cigar, r_read1, (*methy).R[0].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(result1->cigar, r_read1, (*methy).R[0].reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/


			end_5_pos1 = result1->site + r_length1;
		}
		else  ///正向链
		{

			result1->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ1;

			/*************这里要改**************/
			///correct_read_using_cigar(result1->cigar, read1, (*methy).R[0].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(result1->cigar, read1, (*methy).R[0].reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/


			end_5_pos1 = result1->site;
		}


		(*methy).R[0].r_length[(*methy).current_size] = result1->flag;

		(*methy).R[0].sites[(*methy).current_size] = result1->site;

		///(*methy).R1.cut_index[get_cut_id(&(methy->R1), result1->site)]++;


		/**read1的信息**/




		/**read2的信息**/
		/*************这里要改**************/
		///注意r_length一定要+1, 这个是为了放/0的
		///Check_Space((*methy).R[1].reads[(*methy).current_size], (*methy).R[1].r_size[(*methy).current_size], r_length2 + 1);
		Check_Space_3_Bit(
			(*methy).R[1].reads_3_bit[(*methy).current_size],
			(*methy).R[1].r_size_3_bit[(*methy).current_size],
			r_length2);
		(*methy).R[1].r_real_length[(*methy).current_size] = r_length2;
		/*************这里要改**************/

		(*methy).R[1].r_length[(*methy).current_size] = r_length2;


		result2->site = result2->site + _ih_refGenName[result2->chrome_id].start_location - 1;

		if (result2->flag == 16)   ///代表read2的反向互补链比对到了参考组的反向互补上, read2本身实际上比到了参考组正向链上
		{
			result2->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ2;

			/*************这里要改**************/
			///correct_read_using_cigar(result2->cigar, r_read2, (*methy).R[1].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit_over_lap_length_r2(
				result2->cigar, r_read2, (*methy).R[1].reads_3_bit[(*methy).current_size],
				result1->site, r_length1, result2->site, r_length2);
			/*************这里要改**************/


			end_5_pos2 = result2->site;
		}
		else
		{
			result2->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ2;

			/*************这里要改**************/
			////correct_read_using_cigar(result2->cigar, read2, (*methy).R[1].reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit_over_lap_length_r2(
				result2->cigar, read2, (*methy).R[1].reads_3_bit[(*methy).current_size],
				result1->site, r_length1, result2->site, r_length2);
			/*************这里要改**************/


			end_5_pos2 = result2->site + r_length2;
		}


		/*************这里要改**************/
		///over_lap_length(result1->site, r_length1, result2->site, r_length2, (*methy).R[1].reads[(*methy).current_size]);
		/*************这里要改**************/

		(*methy).R[1].r_length[(*methy).current_size] = result2->flag;

		(*methy).R[1].sites[(*methy).current_size] = result2->site;

		///(*methy).R2.cut_index[get_cut_id(&(methy->R2), result2->site)]++;

		/**read2的信息**/

		int block_ID = get_cut_id_pair(methy, get_key_pair_methylation(methy, (*methy).current_size));

		///把小的那个当做key
		(*methy).cut_index[block_ID]++;

		(*methy).current_size++;




		long long end5_distance = end_5_pos2 - end_5_pos1;


		if (end5_distance < 0)
		{
			end5_distance = end5_distance * -1;
		}



		///distances->count[block_ID][end5_distance - distances->minDistance_pair]++;
		distances->count[block_ID][end5_distance]++;


	}


	if ((*methy).current_size == (*methy).total_size)
	{
		assign_cuts_pair(methy);

		push_methy_to_buffer_pair(methy);

		///这个不用清理，因为再push 的时候清理过了
		///clear_methylation_pair(methy);
		
	}

}




///flag == 1是read1
///flag == 2是read2
inline void directly_output_unmapped_PE(char* name, char* read, char* r_read, char* qulity,
	int read_length, bam_output_cell* cell, Output_buffer_sub_block* sub_block, int flag)
{

	if (bam_output == 0)
	{
		if (name[0] == '@')
		{
			fprintf(schema_out_fp, "%s\t", name + 1);
		}
		else
		{
			fprintf(schema_out_fp, "%s\t", name);
		}

		if (flag == 1)
		{
			fprintf(schema_out_fp, "77\t");
		}
		else
		{
			fprintf(schema_out_fp, "141\t");
		}

		fprintf(schema_out_fp, "*\t0\t0\t*\t*\t0\t0\t");

		if (flag == 1)
		{
			fprintf(schema_out_fp, "%s\t", read);
		}
		else
		{
			fprintf(schema_out_fp, "%s\t", r_read);
		}

		fprintf(schema_out_fp, "%s\n", qulity);
	}
	else
	{

		///这个一定一定要有
		sub_block->length = 0;


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

		if (flag == 1)
		{
			sub_block->buffer[sub_block->length] = '7';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '7';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			sub_block->buffer[sub_block->length] = '1';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '4';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '1';
			sub_block->length++;
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}


		///fprintf(schema_out_fp, "*\t0\t0\t*\t*\t0\t0\t");
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
		sub_block->buffer[sub_block->length] = '*';
		sub_block->length++;
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

		if (flag == 1)
		{
			output_to_buffer_char_length(sub_block, read, read_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			output_to_buffer_char_length(sub_block, r_read, read_length);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}

		///这里应该是'\0', 不是'\t'也不是'\n'
		output_to_buffer_char_length(sub_block, qulity, read_length);
		sub_block->buffer[sub_block->length] = '\0';

		write_alignment_directly(sub_block->buffer, sub_block->length, cell);


	}






}




inline void directly_output_read1(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length,
	int paired_end_distance, bam_output_cell* cell,
	Output_buffer_sub_block* sub_block, int output_mask, int mapq)
{

	if (bam_output == 0)
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

		///这里要改
		result->flag = result->flag | output_mask;


		fprintf(schema_out_fp, "%d\t", result->flag);

		fprintf(schema_out_fp, "%s\t", _ih_refGenName[result->chrome_id]._rg_chrome_name);


		fprintf(schema_out_fp, "%llu\t", result->site);

		///fprintf(schema_out_fp, "255\t");
		fprintf(schema_out_fp, "%d\t", mapq);

		

		fprintf(schema_out_fp, "%s\t", result->cigar);





		if (another_result->site > result->site)
		{
			fprintf(schema_out_fp, "=\t%llu\t%llu\t",
				another_result->site, paired_end_distance);

		}
		else if (another_result->site < result->site)
		{
			fprintf(schema_out_fp, "=\t%llu\t-%llu\t",
				another_result->site, paired_end_distance);


		}
		else
		{
			if (matched_length >= another_matched_length)
			{
				fprintf(schema_out_fp, "=\t%llu\t%llu\t",
					another_result->site, paired_end_distance);
			}
			else
			{
				fprintf(schema_out_fp, "=\t%llu\t%llu\t",
					another_result->site, paired_end_distance);
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
	else
	{

		///这个一定一定要有
		sub_block->length = 0;







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

		///这里要改
		result->flag = result->flag | output_mask;

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

		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);

		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;



		output_to_buffer_char_no_length(sub_block, result->cigar);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;





		if (another_result->site > result->site)
		{
			sub_block->buffer[sub_block->length] = '=';
			sub_block->length++;

			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

			output_to_buffer_int(sub_block, another_result->site);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;


			///output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;

		}
		else if (another_result->site < result->site)
		{
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

			///output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
			//上面扩容时有32的余量，所以下面\t不需要检查了
			sub_block->buffer[sub_block->length] = '\t';
			sub_block->length++;
		}
		else
		{
			if (matched_length >= another_matched_length)
			{

				sub_block->buffer[sub_block->length] = '=';
				sub_block->length++;

				sub_block->buffer[sub_block->length] = '\t';
				sub_block->length++;

				output_to_buffer_int(sub_block, another_result->site);
				//上面扩容时有32的余量，所以下面\t不需要检查了
				sub_block->buffer[sub_block->length] = '\t';
				sub_block->length++;

				///output_to_buffer_int(sub_block, matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
				//上面扩容时有32的余量，所以下面\t不需要检查了
				sub_block->buffer[sub_block->length] = '\t';
				sub_block->length++;
			}
			else
			{

				sub_block->buffer[sub_block->length] = '=';
				sub_block->length++;

				sub_block->buffer[sub_block->length] = '\t';
				sub_block->length++;

				output_to_buffer_int(sub_block, another_result->site);
				//上面扩容时有32的余量，所以下面\t不需要检查了
				sub_block->buffer[sub_block->length] = '\t';
				sub_block->length++;


				///output_to_buffer_int(sub_block, another_matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
				//上面扩容时有32的余量，所以下面\t不需要检查了
				sub_block->buffer[sub_block->length] = '\t';
				sub_block->length++;

			}

		}





		if (result->flag & RC_MATE_MAPPED)
		{

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

		///注意这里和多线程原始的算法不一样
		/**
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/

		sub_block->buffer[sub_block->length] = '\0';

		write_alignment_directly(sub_block->buffer, sub_block->length, cell);


	}

	




}









inline void directly_output_read2_return_buffer(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length,
	Output_buffer_sub_block* sub_block, int paired_end_distance,
	bam_phrase* bam_groups, int output_mask, int mapq)
{

	if (bam_output == 0)
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

		result->flag = result->flag | output_mask;

		

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

		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);
	

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


			///output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

			///output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

				///output_to_buffer_int(sub_block, matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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


				///output_to_buffer_int(sub_block, another_matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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
	else
	{



		///这个一定一定要有
		sub_block->length = 0;


		///这个和正向read是相反的
		if (result->flag == 0)
		{
			result->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ2;
		}
		else
		{
			result->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ2;
		}

		result->flag = result->flag | output_mask;

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
		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);

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


			///output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

			///output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

				///output_to_buffer_int(sub_block, matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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


				///output_to_buffer_int(sub_block, another_matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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
		/**
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/

		sub_block->buffer[sub_block->length] = '\0';


		convert_string_to_bam(sub_block->buffer, sub_block->length, bam_groups);

	}
}








inline void directly_output_read2(char* name, char* read, char* r_read, char* qulity,
	map_result* result, map_result* another_result, int read_length, int matched_length, int another_matched_length,
	int paired_end_distance, bam_output_cell* cell,
	Output_buffer_sub_block* sub_block, int output_mask, int mapq)
{

	if (bam_output == 0)
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

		///这里要改
		result->flag = result->flag | output_mask;


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

		///fprintf(schema_out_fp, "255\t");
		fprintf(schema_out_fp, "%d\t", mapq);


		fprintf(schema_out_fp, "%s\t", result->cigar);

		///fprintf(schema_out_fp, "*\t0\t0\t");
		///fprintf(schema_out_fp, "=\t%llu\t%lld\t", another_result->site, (long long)(another_result->site) - (long long)(result->site));




		if (another_result->site > result->site)
		{
			fprintf(schema_out_fp, "=\t%llu\t%llu\t",
				another_result->site, paired_end_distance);

		}
		else if (another_result->site < result->site)
		{
			fprintf(schema_out_fp, "=\t%llu\t-%llu\t",
				another_result->site, paired_end_distance);
		}
		else
		{
			if (matched_length >= another_matched_length)
			{
				fprintf(schema_out_fp, "=\t%llu\t-%llu\t",
					another_result->site, paired_end_distance);
			}
			else
			{
				fprintf(schema_out_fp, "=\t%llu\t-%llu\t",
					another_result->site, paired_end_distance);
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
	else
	{


		///这个一定一定要有
		sub_block->length = 0;


















		///这个和正向read是相反的
		if (result->flag == 0)
		{
			result->flag = IS_PAIRED | MATE_MAPPED | RC_MAPPED | READ2;
		}
		else
		{
			result->flag = IS_PAIRED | MATE_MAPPED | RC_MATE_MAPPED | READ2;
		}

		///这里要改
		result->flag = result->flag | output_mask;


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

		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);
	

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


			///output_to_buffer_int(sub_block, another_result->site - result->site + another_matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

			///output_to_buffer_int(sub_block, result->site - another_result->site + matched_length);
			output_to_buffer_int(sub_block, paired_end_distance);
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

				///output_to_buffer_int(sub_block, matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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


				///output_to_buffer_int(sub_block, another_matched_length);
				output_to_buffer_int(sub_block, paired_end_distance);
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


		///注意这里和多线程原始的算法不一样
		/**
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/
		sub_block->buffer[sub_block->length] = '\0';

		write_alignment_directly(sub_block->buffer, sub_block->length, cell);
	}

	




}





inline void output_sam_end_to_end(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length,
	bam_output_cell* cell,
	Output_buffer_sub_block* sub_block,
	int* map_among_references, int output_mask, int mapq)
{
	(*map_among_references) = 0;

	if (bam_output == 0)
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


		///bitmapper_bs_iter r_length = end_site - start_site + 1;
		///end_site - start_site + 1这个是read在基因组上的alignment长度
		///map_location此刻是1-based
		///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) < _ih_refGenName[now_ref_name]._rg_chrome_length )
		if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
		{

			(*map_among_references) = 1;
			return;
		}


		if (name[0] == '@')
		{
			fprintf(schema_out_fp, "%s\t", name + 1);
		}
		else
		{
			fprintf(schema_out_fp, "%s\t", name);
		}


		fprintf(schema_out_fp, "%d\t", flag | output_mask);

		fprintf(schema_out_fp, "%s\t", _ih_refGenName[now_ref_name]._rg_chrome_name);


		fprintf(schema_out_fp, "%llu\t", map_location);

		fprintf(schema_out_fp, "%d\t", mapq);

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


		fprintf(schema_out_fp, "NM:i:%d\n", err);
	}
	else
	{


		///这个一定一定要有
		sub_block->length = 0;












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



		///bitmapper_bs_iter r_length = end_site - start_site + 1;
		///end_site - start_site + 1这个是read在基因组上的alignment长度
		///map_location此刻是1-based
		///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) < _ih_refGenName[now_ref_name]._rg_chrome_length )
		if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
		{
			(*map_among_references) = 1;
			return;
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






		///result_string << flag << "\t";
		output_to_buffer_int(sub_block, flag | output_mask);
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
		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);
	

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

		///注意这里和多线程原始的算法不一样
		/**
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/
		sub_block->buffer[sub_block->length] = '\0';

		write_alignment_directly(sub_block->buffer, sub_block->length, cell);

	}

	



}





inline void output_methy_end_to_end(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length,
	Methylation* methy,
	int* map_among_references)
{
	int flag;

	(*map_among_references) = 0;

	if ((*methy).current_size < (*methy).total_size)
	{

		///只要site >= _msf_refGenLength，就能判断出是正向还是反向

		///矫正之后的read长度就是这个
		bitmapper_bs_iter r_length = end_site - start_site + 1;

		/*************这里要改**************/
		///注意r_length一定要+1, 这个是为了放/0的
		///Check_Space((*methy).reads[(*methy).current_size], (*methy).r_size[(*methy).current_size], r_length+1);
		Check_Space_3_Bit((*methy).reads_3_bit[(*methy).current_size], (*methy).r_size_3_bit[(*methy).current_size],
			r_length);
		(*methy).r_real_length[(*methy).current_size] = r_length;
		/*************这里要改**************/

		(*methy).r_length[(*methy).current_size] = r_length;

		if (site >= _msf_refGenLength)
		{

			site = site + end_site;

			site = _msf_refGenLength * 2 - site - 1;


			/**************************这个是加的东西*******************************/
			int now_ref_name;

			for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
			{
				if ((site >= _ih_refGenName[now_ref_name].start_location) && (site <= _ih_refGenName[now_ref_name].end_location))
				{

					///bitmapper_bs_iter r_length = end_site - start_site + 1;
					///end_site - start_site + 1这个是read在基因组上的alignment长度
					if (site + end_site - start_site > _ih_refGenName[now_ref_name].end_location)
					{
						(*map_among_references) = 1;
						return;
					}


					break;
				}
					
			}
			/**************************这个是加的东西*******************************/





			/*************这里要改**************/
			///correct_read_using_cigar(best_cigar, r_read, (*methy).reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(best_cigar, r_read, (*methy).reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/
			

			flag = 16;

		}
		else
		{
			site = site + start_site;



			/**************************这个是加的东西*******************************/
			int now_ref_name;

			for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
			{
				if ((site >= _ih_refGenName[now_ref_name].start_location) && (site <= _ih_refGenName[now_ref_name].end_location))
				{

					///bitmapper_bs_iter r_length = end_site - start_site + 1;
					///end_site - start_site + 1这个是read在基因组上的alignment长度
					if (site + end_site - start_site > _ih_refGenName[now_ref_name].end_location)
					{
						(*map_among_references) = 1;
						return;
					}


					break;
				}

			}
			/**************************这个是加的东西*******************************/




			/*************这里要改**************/
			///correct_read_using_cigar(best_cigar, read, (*methy).reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(best_cigar, read, (*methy).reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/

			flag = 0;
		}



		(*methy).r_length[(*methy).current_size] = flag;

		(*methy).sites[(*methy).current_size] = site;

		(*methy).cut_index[get_cut_id(methy, site)]++;

		/**
		fprintf(stderr, "cut_length: %d\n", (*methy).cut_length);
		fprintf(stderr, "get_cut_id(methy, site): %d\n", get_cut_id(methy, site));
		**/

		(*methy).current_size++;

	}

	if ((*methy).current_size == (*methy).total_size)
	{
		assign_cuts(methy);

		output_methylation(methy);

		clear_methylation(methy);
	}




}





inline void output_methy_end_to_end_output_buffer(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length, Methylation* methy,
	int* map_among_references)
{

	(*map_among_references) = 0;


	int flag;

	if ((*methy).current_size < (*methy).total_size)
	{

		///只要site >= _msf_refGenLength，就能判断出是正向还是反向

		///矫正之后的read长度就是这个
		bitmapper_bs_iter r_length = end_site - start_site + 1;





		/*************这里要改**************/
		///注意r_length一定要+1, 这个是为了放/0的
		///Check_Space((*methy).reads[(*methy).current_size], (*methy).r_size[(*methy).current_size], r_length + 1);
		Check_Space_3_Bit((*methy).reads_3_bit[(*methy).current_size], (*methy).r_size_3_bit[(*methy).current_size],
			r_length);
		(*methy).r_real_length[(*methy).current_size] = r_length;

		/*************这里要改**************/













		(*methy).r_length[(*methy).current_size] = r_length;

		if (site >= _msf_refGenLength)
		{

			site = site + end_site;

			site = _msf_refGenLength * 2 - site - 1;






			/**************************这个是加的东西*******************************/
			int now_ref_name;

			for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
			{
				if ((site >= _ih_refGenName[now_ref_name].start_location) && (site <= _ih_refGenName[now_ref_name].end_location))
				{

					///bitmapper_bs_iter r_length = end_site - start_site + 1;
					///end_site - start_site + 1这个是read在基因组上的alignment长度
					if (site + end_site - start_site > _ih_refGenName[now_ref_name].end_location)
					{
						(*map_among_references) = 1;
						return;
					}


					break;
				}

			}
			/**************************这个是加的东西*******************************/








			/*************这里要改**************/
			///correct_read_using_cigar(best_cigar, r_read, (*methy).reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(best_cigar, r_read, (*methy).reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/

			flag = 16;

		}
		else
		{
			site = site + start_site;





			/**************************这个是加的东西*******************************/
			int now_ref_name;

			for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
			{
				if ((site >= _ih_refGenName[now_ref_name].start_location) && (site <= _ih_refGenName[now_ref_name].end_location))
				{

					///bitmapper_bs_iter r_length = end_site - start_site + 1;
					///end_site - start_site + 1这个是read在基因组上的alignment长度
					if (site + end_site - start_site > _ih_refGenName[now_ref_name].end_location)
					{
						(*map_among_references) = 1;
						return;
					}


					break;
				}

			}
			/**************************这个是加的东西*******************************/






			/*************这里要改**************/
			///correct_read_using_cigar(best_cigar, read, (*methy).reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(best_cigar, read, (*methy).reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/

			flag = 0;
		}

		(*methy).r_length[(*methy).current_size] = flag;

		(*methy).sites[(*methy).current_size] = site;

		(*methy).cut_index[get_cut_id(methy, site)]++;

		(*methy).current_size++;

	}

	if ((*methy).current_size == (*methy).total_size)
	{
		assign_cuts(methy);

		push_methy_to_buffer(methy);

		///这个不用清理，因为再push 的时候清理过了
		///clear_methylation(methy);

	}




}



inline void output_sam_end_to_end_output_buffer(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length, Output_buffer_sub_block* sub_block,
	bam_phrase* bam_groups, int* map_among_references, int output_mode, int mapq)
{
	(*map_among_references) = 0;

	if (bam_output == 0)
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


		///bitmapper_bs_iter r_length = end_site - start_site + 1;
		///end_site - start_site + 1这个是read在基因组上的alignment长度
		///map_location此刻是1-based
		///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) < _ih_refGenName[now_ref_name]._rg_chrome_length )
		if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
		{
			(*map_among_references) = 1;
			return;
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


		///result_string << flag << "\t";
		output_to_buffer_int(sub_block, flag | output_mode);
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

		/**
		///result_string << "255\t";
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);
		//上面扩容时有32的余量，所以下面\t不需要检查了
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
	else
	{
		

			///这个一定一定要有
			sub_block->length = 0;


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



			///bitmapper_bs_iter r_length = end_site - start_site + 1;
			///end_site - start_site + 1这个是read在基因组上的alignment长度
			///map_location此刻是1-based
			///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) < _ih_refGenName[now_ref_name]._rg_chrome_length )
			if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
			{

				(*map_among_references) = 1;
				return;
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






			///result_string << flag << "\t";
			output_to_buffer_int(sub_block, flag | output_mode);
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
			/**
			sub_block->buffer[sub_block->length] = '2';
			sub_block->length++;

			sub_block->buffer[sub_block->length] = '5';
			sub_block->length++;

			sub_block->buffer[sub_block->length] = '5';
			sub_block->length++;
			**/
			output_to_buffer_int(sub_block, mapq);
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

			///注意这里和多线程原始的算法不一样
			
			//上面扩容时有32的余量，所以下面\t不需要检查了
			//sub_block->buffer[sub_block->length] = '\n';
			//sub_block->length++;
			
			sub_block->buffer[sub_block->length] = '\0';


			convert_string_to_bam(sub_block->buffer, sub_block->length, bam_groups);


	}

	
}








inline void output_sam_end_to_end_output_buffer_unmapped(char* name, char* read, char* r_read, char* qulity,
	int read_length, Output_buffer_sub_block* sub_block, bam_phrase* bam_groups)
{

	if (bam_output == 0)
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



		///fprintf(schema_out_fp, "4\t*\t0\t0\t*\t*\t0\t0\t");
		sub_block->buffer[sub_block->length] = '4';
		sub_block->length++;
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
		sub_block->buffer[sub_block->length] = '*';
		sub_block->length++;
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



		///result_string << read << "\t";
		output_to_buffer_char_length(sub_block, read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		///result_string << qulity << "\t";
		output_to_buffer_char_length(sub_block, qulity, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;




	}
	else
	{


		///这个一定一定要有
		sub_block->length = 0;




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







		///fprintf(schema_out_fp, "4\t*\t0\t0\t*\t*\t0\t0\t");
		sub_block->buffer[sub_block->length] = '4';
		sub_block->length++;
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
		sub_block->buffer[sub_block->length] = '*';
		sub_block->length++;
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



		///result_string << read << "\t";
		output_to_buffer_char_length(sub_block, read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		///result_string << qulity << "\t";
		output_to_buffer_char_length(sub_block, qulity, read_length);


		sub_block->buffer[sub_block->length] = '\0';


		convert_string_to_bam(sub_block->buffer, sub_block->length, bam_groups);


	}


}











inline void output_sam_end_to_end_pbat_output_buffer(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length, Output_buffer_sub_block* sub_block,
	bam_phrase* bam_groups, int* map_among_references, int output_mode, int mapq)
{

	(*map_among_references) = 0;

	if (bam_output == 0)
	{

		bitmapper_bs_iter map_location = site;

		int flag;

		if (map_location >= _msf_refGenLength)
		{

			map_location = map_location + end_site;

			map_location = _msf_refGenLength * 2 - map_location - 1;

			///flag = 0;
			flag = 16;

		}
		else
		{
			map_location = map_location + start_site;

			///flag = 16;
			flag = 0;
		}


		int now_ref_name = 0;


		for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
		{
			if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
				break;
		}


		map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;




		///bitmapper_bs_iter r_length = end_site - start_site + 1;
		///end_site - start_site + 1这个是read在基因组上的alignment长度
		///map_location此刻是1-based
		///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) < _ih_refGenName[now_ref_name]._rg_chrome_length )
		if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
		{

			(*map_among_references) = 1;
			return;
		}



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


		output_to_buffer_int(sub_block, flag | output_mode);
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

		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);
	

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




		///if (flag == 16)
		if (flag == 0)
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
	else
	{



		///这个一定一定要有
		sub_block->length = 0;


		


		bitmapper_bs_iter map_location = site;

		int flag;

		if (map_location >= _msf_refGenLength)
		{

			map_location = map_location + end_site;

			map_location = _msf_refGenLength * 2 - map_location - 1;

			///flag = 0;
			flag = 16;

		}
		else
		{
			map_location = map_location + start_site;

			///flag = 16;
			flag = 0;
		}


		int now_ref_name = 0;


		for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
		{
			if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
				break;
		}


		map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;



		///bitmapper_bs_iter r_length = end_site - start_site + 1;
		///end_site - start_site + 1这个是read在基因组上的alignment长度
		///map_location此刻是1-based
		///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) < _ih_refGenName[now_ref_name]._rg_chrome_length )
		if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
		{

			(*map_among_references) = 1;
			return;
		}


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


		output_to_buffer_int(sub_block, flag | output_mode);
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

		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);
	

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




		///if (flag == 16)
		if (flag == 0)
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
		/**
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/

		sub_block->buffer[sub_block->length] = '\0';

		convert_string_to_bam(sub_block->buffer, sub_block->length, bam_groups);


	}


}










inline void output_sam_end_to_end_pbat(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length,
	bam_output_cell* cell,
	Output_buffer_sub_block* sub_block, int* map_among_references, int output_mode, int mapq)
{

	(*map_among_references) = 0;


	if (bam_output == 0)
	{
		///fprintf(stderr, "******************\n");

		bitmapper_bs_iter map_location = site;

		int flag;

		if (map_location >= _msf_refGenLength)
		{

			map_location = map_location + end_site;

			map_location = _msf_refGenLength * 2 - map_location - 1;

			flag = 16;
			///flag = 0;

		}
		else
		{
			map_location = map_location + start_site;

			flag = 0;
			///flag = 16;
		}


		int now_ref_name = 0;


		for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
		{
			if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
				break;
		}


		map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;




		///bitmapper_bs_iter r_length = end_site - start_site + 1;
		///end_site - start_site + 1这个是read在基因组上的alignment长度
		///map_location此刻是1-based
		///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) > _ih_refGenName[now_ref_name]._rg_chrome_length )
		if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
		{

			(*map_among_references) = 1;
			return;
		}











		if (name[0] == '@')
		{
			fprintf(schema_out_fp, "%s\t", name + 1);
		}
		else
		{
			fprintf(schema_out_fp, "%s\t", name);
		}


		fprintf(schema_out_fp, "%d\t", flag | output_mode);

		fprintf(schema_out_fp, "%s\t", _ih_refGenName[now_ref_name]._rg_chrome_name);


		fprintf(schema_out_fp, "%llu\t", map_location);

		///fprintf(schema_out_fp, "255\t");
		fprintf(schema_out_fp, "%d\t", mapq);

		fprintf(schema_out_fp, "%s\t", best_cigar);

		fprintf(schema_out_fp, "*\t0\t0\t");



		///if (flag == 16)
		if (flag == 0)
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
	else
	{

		///这个一定一定要有
		sub_block->length = 0;








		bitmapper_bs_iter map_location = site;

		int flag;

		if (map_location >= _msf_refGenLength)
		{

			map_location = map_location + end_site;

			map_location = _msf_refGenLength * 2 - map_location - 1;

			///flag = 0;
			flag = 16;

		}
		else
		{
			map_location = map_location + start_site;

			///flag = 16;
			flag = 0;
		}


		int now_ref_name = 0;


		for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
		{
			if ((map_location >= _ih_refGenName[now_ref_name].start_location) && (map_location <= _ih_refGenName[now_ref_name].end_location))
				break;
		}


		map_location = map_location + 1 - _ih_refGenName[now_ref_name].start_location;



		///bitmapper_bs_iter r_length = end_site - start_site + 1;
		///end_site - start_site + 1这个是read在基因组上的alignment长度
		///map_location此刻是1-based
		///所以原来应该是if ((map_location -1) + (end_site - start_site + 1) > _ih_refGenName[now_ref_name]._rg_chrome_length )
		if (map_location + end_site - start_site > _ih_refGenName[now_ref_name]._rg_chrome_length)
		{

			(*map_among_references) = 1;
			return;
		}




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


		output_to_buffer_int(sub_block, flag | output_mode);
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


		/**
		sub_block->buffer[sub_block->length] = '2';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;

		sub_block->buffer[sub_block->length] = '5';
		sub_block->length++;
		**/
		output_to_buffer_int(sub_block, mapq);

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




		///if (flag == 16)
		if (flag == 0)
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
		/**
		///注意这里和多线程原始的算法不一样
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/
		sub_block->buffer[sub_block->length] = '\0';

		write_alignment_directly(sub_block->buffer, sub_block->length, cell);
	}


}



inline void output_methy_end_to_end_pbat(char* name, char* read, char* r_read, char* qulity,
	bitmapper_bs_iter site,
	bitmapper_bs_iter end_site,
	bitmapper_bs_iter start_site,
	int err, char* best_cigar, int read_length,
	Methylation* methy,
	int* map_among_references)
{


	(*map_among_references) = 0;

	int flag;

	if ((*methy).current_size < (*methy).total_size)
	{

		///只要site >= _msf_refGenLength，就能判断出是正向还是反向

		///矫正之后的read长度就是这个
		bitmapper_bs_iter r_length = end_site - start_site + 1;

		/*************这里要改**************/
		///注意r_length一定要+1, 这个是为了放/0的
		///Check_Space((*methy).reads[(*methy).current_size], (*methy).r_size[(*methy).current_size], r_length + 1);
		Check_Space_3_Bit((*methy).reads_3_bit[(*methy).current_size], (*methy).r_size_3_bit[(*methy).current_size],
			r_length);
		(*methy).r_real_length[(*methy).current_size] = r_length;
		/*************这里要改**************/

		(*methy).r_length[(*methy).current_size] = r_length;

		if (site >= _msf_refGenLength)
		{

			site = site + end_site;

			site = _msf_refGenLength * 2 - site - 1;







			/**************************这个是加的东西*******************************/
			int now_ref_name;

			for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
			{
				if ((site >= _ih_refGenName[now_ref_name].start_location) && (site <= _ih_refGenName[now_ref_name].end_location))
				{

					///bitmapper_bs_iter r_length = end_site - start_site + 1;
					///end_site - start_site + 1这个是read在基因组上的alignment长度
					if (site + end_site - start_site > _ih_refGenName[now_ref_name].end_location)
					{
						(*map_among_references) = 1;
						return;
					}


					break;
				}

			}
			/**************************这个是加的东西*******************************/










			/*************这里要改**************/
			///correct_read_using_cigar(best_cigar, r_read, (*methy).reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(best_cigar, r_read, (*methy).reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/

			flag = 16;

		}
		else
		{
			site = site + start_site;





			/**************************这个是加的东西*******************************/
			int now_ref_name;

			for (now_ref_name = 0; now_ref_name < refChromeCont; ++now_ref_name)
			{
				if ((site >= _ih_refGenName[now_ref_name].start_location) && (site <= _ih_refGenName[now_ref_name].end_location))
				{

					///bitmapper_bs_iter r_length = end_site - start_site + 1;
					///end_site - start_site + 1这个是read在基因组上的alignment长度
					if (site + end_site - start_site > _ih_refGenName[now_ref_name].end_location)
					{
						(*map_among_references) = 1;
						return;
					}


					break;
				}

			}
			/**************************这个是加的东西*******************************/












			/*************这里要改**************/
			////correct_read_using_cigar(best_cigar, read, (*methy).reads[(*methy).current_size]);
			correct_read_using_cigar_3_bit(best_cigar, read, (*methy).reads_3_bit[(*methy).current_size]);
			/*************这里要改**************/

			flag = 0;
		}

		(*methy).r_length[(*methy).current_size] = flag;

		(*methy).sites[(*methy).current_size] = site;

		(*methy).cut_index[get_cut_id(methy, site)]++;

		(*methy).current_size++;

	}

	if ((*methy).current_size == (*methy).total_size)
	{
		assign_cuts(methy);

		output_methylation(methy);

		clear_methylation(methy);
	}



}








inline void calculate_cigar_end_to_end_back(
	bitmapper_bs_iter site,
	char* path, int path_length,
	char* cigar,
	int t_length,
	int* start_site,
	bitmapper_bs_iter* end_site,
	bitmapper_bs_iter errthold,
	char *pattern, int p_length,
	char *text)
{

	int i = 0;
	char pre_ciga = 100;
	int pre_ciga_length = 0;




	cigar[0] = '\0';

	pre_ciga = 100;
	pre_ciga_length = 0;

	///调整cigar, 其实就是把cigar两边的I转成M
	for (i = 0; i < path_length; i++)
	{
		if (map_cigar[path[i]] != 'I')
		{
			break;
		}
		else
		{
			path[i] = 1;
			///(*start_site)--;
			(*end_site)++;
		}
	}


	for (i = path_length - 1; i >= 0; i--)
	{
		if (map_cigar[path[i]] != 'I')
		{
			break;
		}
		else
		{
			path[i] = 1;
			///(*end_site)++;
			(*start_site)--;
		}
	}
	///调整cigar, 其实就是把cigar两边的I转成M



	int no_gap_site = t_length - 1 + errthold;



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


inline int mismatch_score(char *pattern, int p_length,
	char *text, int t_length, int errthold)
{

	int i = 0;
	int tmp_err = 0;


	int start_site = errthold;


	for (i = 0; i < t_length; i++)
	{
		if (pattern[i + start_site] == 'N')
		{
			tmp_err = tmp_err + 1;
		}
		else if (text[i] != pattern[i + start_site])
		{

			if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
			{
				tmp_err = tmp_err + 6;
			}

		}
	}

	return tmp_err;

}



inline void calculate_cigar_end_to_end(
	bitmapper_bs_iter site,
	char* path, int path_length,
	char* cigar,
	int t_length,
	int* start_site,
	bitmapper_bs_iter* end_site,
	bitmapper_bs_iter errthold,
	char *pattern, int p_length,
	char *text,
	unsigned int* error)
{



	



	


	int i = 0;
	char pre_ciga = 100;
	int pre_ciga_length = 0;




	cigar[0] = '\0';

	pre_ciga = 100;
	pre_ciga_length = 0;

	///调整cigar, 其实就是把cigar两边的I转成M
	for (i = 0; i < path_length; i++)
	{
		if (map_cigar[path[i]] != 'I')
		{
			break;
		}
		else
		{
			path[i] = 1;
			///(*start_site)--;
			(*end_site)++;
		}
	}


	for (i = path_length - 1; i >= 0; i--)
	{
		if (map_cigar[path[i]] != 'I')
		{
			break;
		}
		else
		{
			path[i] = 1;
			///(*end_site)++;
			(*start_site)--;
		}
	}
	///调整cigar, 其实就是把cigar两边的I转成M



	int score = 0;
	int gap_length = 0;
	int gap_open = 5;
	int gap_extend = 3;
	int np = 1;
	int mp = 6;


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

					if (pre_ciga != 'M')
					{
						gap_length = gap_length + pre_ciga_length;
						score = score + gap_open + gap_extend*pre_ciga_length;
					}
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

			if (pre_ciga != 'M')
			{
				gap_length = gap_length + pre_ciga_length;
				score = score + gap_open + gap_extend*pre_ciga_length;
			}
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

					if (pre_ciga != 'M')
					{
						gap_length = gap_length + pre_ciga_length;
						score = score + gap_open + gap_extend*pre_ciga_length;
					}
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

			if (pre_ciga != 'M')
			{
				gap_length = gap_length + pre_ciga_length;
				score = score + gap_open + gap_extend*pre_ciga_length;
			}
		}
	}


	
	score = score + ((*error) - gap_length)*mp;

	
	///trim_Ns_new(pattern, p_length, text, t_length, errthold, error, start_site, end_site, *end_site, score, mp, np, cigar);
	

	
	trim_Ns_new_all(pattern, p_length, text, t_length, errthold, error, start_site, end_site,
		*end_site, score, mp, np, cigar);
	

}







inline int calculate_best_map_cigar_end_to_end_return(bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar, map_result* result, int* matched_length, int need_reverse_quality)
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


		/**
		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(result->err),
			cigar, result->end_site, &start_site, path, &path_length, matrix_bit, &result->end_site) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			///double start = Get_T();
			calculate_cigar_end_to_end(result->origin_site, path, path_length, best_cigar, t_length, &start_site, &result->end_site, 
				error_threshold, tmp_ref, p_length, read, &(result->err));
			///total = total + Get_T() - start;
		}
		**/
		///need modification
		fast_recalculate_bs_Cigar(tmp_ref, p_length, read, t_length, error_threshold, result->end_site,
		result->err, &start_site, &result->end_site, &result->err, &result->score, best_cigar, 
		result->origin_site < _msf_refGenLength, mat, mat_diff, GapOpenPenalty, GapExtensionPenalty, 
		MistMatchPenaltyMax, MistMatchPenaltyMin, N_Penalty, qulity, need_reverse_quality, Q_base);


	}
	else
	{

		///start_site = 0;
		///need modification
		result->score = 0;

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
	char* best_cigar, Output_buffer_sub_block* sub_block, Methylation* methy,
	bam_phrase* bam_groups, int* map_among_references, int output_mode, unsigned int second_best_diff)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;


	///need modification
	int score;

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





		/**
		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit, &candidate_votes[0].end_site) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length, &start_site, &candidate_votes[0].end_site, 
				error_threshold, tmp_ref, p_length, read, &(candidate_votes[0].err));
		}
		**/

		///need modification
		fast_recalculate_bs_Cigar(tmp_ref, p_length, read, t_length, error_threshold, candidate_votes[0].end_site,
		candidate_votes[0].err, &start_site, &candidate_votes[0].end_site, &candidate_votes[0].err, &score, best_cigar, 
		candidate_votes[0].site < _msf_refGenLength, mat, mat_diff, GapOpenPenalty, GapExtensionPenalty, 
		MistMatchPenaltyMax, MistMatchPenaltyMin, N_Penalty, qulity, 0, Q_base);

	}
	else
	{
		///need modification
		score = 0;

		///start_site = 0;
		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}


	///need modification
	int mapq = MAP_Calculation(second_best_diff, error_threshold, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
	/**
	fprintf(stderr, "mapq: %d, second_best_diff: %llu, candidate_votes[0].err: %d, score: %d, error_threshold: %d, gapoe: %d, mis: %d\n", 
	mapq, second_best_diff, candidate_votes[0].err, score, error_threshold, GapOpenPenalty+GapExtensionPenalty, MistMatchPenaltyMax);
	**/


	if (output_methy == 1)
	{
		output_methy_end_to_end_output_buffer
			(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length, methy, map_among_references);
	}
	else
	{
		output_sam_end_to_end_output_buffer
			(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length, sub_block, bam_groups, map_among_references, output_mode, mapq);
	}

	




	return 0;




}







inline int calculate_best_map_cigar_end_to_end_pbat_output_buffer(seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar, Output_buffer_sub_block* sub_block, Methylation* methy,
	bam_phrase* bam_groups, int* map_among_references, int output_mode, unsigned int second_best_diff)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;



	int start_site;
	int path_length;

	///need modification
	int score;




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





		/**
		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit, &candidate_votes[0].end_site) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length, &start_site, &candidate_votes[0].end_site, 
				error_threshold, tmp_ref, p_length, read, &(candidate_votes[0].err));
		}
		**/
		///need modification
		///quality need to be reverse
		fast_recalculate_bs_Cigar(tmp_ref, p_length, read, t_length, error_threshold, candidate_votes[0].end_site,
		candidate_votes[0].err, &start_site, &candidate_votes[0].end_site, &candidate_votes[0].err, &score, best_cigar, 
		candidate_votes[0].site < _msf_refGenLength, mat, mat_diff, GapOpenPenalty, GapExtensionPenalty, 
		MistMatchPenaltyMax, MistMatchPenaltyMin, N_Penalty, qulity, 1, Q_base);

	}
	else
	{

		///start_site = 0;
		///need modification
		score = 0;

		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}

	///need modification
	int mapq = MAP_Calculation(second_best_diff, error_threshold, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);


	if (output_methy == 1)
	{
		output_methy_end_to_end_output_buffer
			(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length, methy, map_among_references);
	}
	else
	{
		output_sam_end_to_end_pbat_output_buffer(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length, sub_block, bam_groups, map_among_references, output_mode,mapq);

	}

	


	return 0;




}





inline int calculate_best_map_cigar_end_to_end(seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar, Methylation* methy, 
	bam_output_cell* cell,
	Output_buffer_sub_block* current_sub_buffer,
	int *map_among_references, int output_mask, unsigned int second_best_diff)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;
	
	///need modification
	int score;


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


		/**
		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit, &candidate_votes[0].end_site) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length, &start_site, &candidate_votes[0].end_site, 
				error_threshold, tmp_ref, p_length, read, &(candidate_votes[0].err));
		}
		**/
	
		///need modification
		fast_recalculate_bs_Cigar(tmp_ref, p_length, read, t_length, error_threshold, candidate_votes[0].end_site,
		candidate_votes[0].err, &start_site, &candidate_votes[0].end_site, &candidate_votes[0].err, &score, best_cigar, 
		candidate_votes[0].site < _msf_refGenLength, mat, mat_diff, GapOpenPenalty, GapExtensionPenalty, 
		MistMatchPenaltyMax, MistMatchPenaltyMin, N_Penalty, qulity, 0, Q_base);


	}
	else
	{
		///start_site = 0;
		///need modification
		score = 0;
		
		
		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);
	}



	///need modification
	int mapq = MAP_Calculation(second_best_diff, error_threshold, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
	/**
	fprintf(stderr, "mapq: %d, second_best_diff: %llu, candidate_votes[0].err: %d, score: %d, error_threshold: %d, gapoe: %d, mis: %d\n", 
	mapq, second_best_diff, candidate_votes[0].err, score, error_threshold, GapOpenPenalty+GapExtensionPenalty, MistMatchPenaltyMax);
	**/


	if (output_methy == 1)
	{
		output_methy_end_to_end(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length, 
			methy, map_among_references);
	}
	else
	{
		output_sam_end_to_end(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length,
			cell,
			current_sub_buffer,
			map_among_references, output_mask, mapq);
	}
	


	return 0;




}






inline int calculate_best_map_cigar_end_to_end_pbat(seed_votes* candidate_votes, bitmapper_bs_iter* candidate_votes_length,
	bitmapper_bs_iter read_length, bitmapper_bs_iter error_threshold, char* tmp_ref,
	char* read, char* r_read, char* qulity, char* cigar, char* path, uint16_t* matrix, Word* matrix_bit, char* name,
	char* best_cigar, Methylation* methy,
	bam_output_cell* cell,
	Output_buffer_sub_block* current_sub_buffer, int* map_among_references, int output_mode,
	unsigned int second_best_diff)
{
	bitmapper_bs_iter t_length = read_length;
	bitmapper_bs_iter p_length = read_length + 2 * error_threshold;
	bitmapper_bs_iter site;



	int start_site;
	int path_length;
	///need modification
	int score;




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



		/**
		if (fast_bs_Calculate_Cigar(tmp_ref, p_length, read, t_length, error_threshold, &(candidate_votes[0].err),
			cigar, candidate_votes[0].end_site, &start_site, path, &path_length, matrix_bit, &candidate_votes[0].end_site) == 6)
		{
			sprintf(best_cigar, "%dM", t_length);
		}
		else
		{
			///double start = Get_T();
			calculate_cigar_end_to_end(candidate_votes[0].site, path, path_length, best_cigar, t_length, &start_site, &candidate_votes[0].end_site, 
				error_threshold, tmp_ref, p_length, read, &(candidate_votes[0].err));
			///total = total + Get_T() - start;
		}
		**/
		///need modification
		///quality need to be reverse
		fast_recalculate_bs_Cigar(tmp_ref, p_length, read, t_length, error_threshold, candidate_votes[0].end_site,
		candidate_votes[0].err, &start_site, &candidate_votes[0].end_site, &candidate_votes[0].err, &score, best_cigar, 
		candidate_votes[0].site < _msf_refGenLength, mat, mat_diff, GapOpenPenalty, GapExtensionPenalty, 
		MistMatchPenaltyMax, MistMatchPenaltyMin, N_Penalty, qulity, 1, Q_base);

	}
	else
	{

		///start_site = 0;
		///need modification
		score = 0;

		start_site = candidate_votes[0].end_site - t_length + 1;
		sprintf(best_cigar, "%dM", t_length);

	}


	///need modification
	int mapq = MAP_Calculation(second_best_diff, error_threshold, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);


	if (output_methy == 1)
	{
		output_methy_end_to_end_pbat(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length, methy, map_among_references);
	}
	else
	{
		output_sam_end_to_end_pbat(name, read, r_read, qulity,
			candidate_votes[0].site,
			candidate_votes[0].end_site,
			start_site,
			candidate_votes[0].err, best_cigar,
			read_length,
			cell,
			current_sub_buffer, map_among_references, output_mode, mapq);
	}

	



	return 0;




}








inline int try_process_unique_mismatch_end_to_end(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, int* get_error,
	Methylation* methy,
	bam_output_cell* cell,
	Output_buffer_sub_block* current_sub_buffer,
	int* map_among_references)
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

		if (output_methy == 1)
		{
			output_methy_end_to_end(name, read, r_read, qulity,
				locates[0],
				read_length - 1,
				0,
				0,
				cigar, read_length, methy, map_among_references);
		}
		else
		{
			output_sam_end_to_end(name, read, r_read, qulity,
				locates[0],
				read_length - 1,
				0,
				0,
				cigar, read_length,
				cell,
				current_sub_buffer,
				map_among_references, 0, 42);
		}

		

		return 1;

	}

	return 0;



}



inline int try_process_unique_mismatch_end_to_end_pbat(
		bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
		bitmapper_bs_iter* locates,
		bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, int* get_error, 
		Methylation* methy,
		bam_output_cell* cell,
		Output_buffer_sub_block* current_sub_buffer,
		int *map_among_references)
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

			if (output_methy == 1)
			{
				output_methy_end_to_end_pbat(name, read, r_read, qulity,
					locates[0],
					read_length - 1,
					0,
					0,
					cigar, read_length, methy, map_among_references);
			}
			else
			{
				output_sam_end_to_end_pbat(name, read, r_read, qulity,
					locates[0],
					read_length - 1,
					0,
					0,
					cigar, read_length,
					cell,
					current_sub_buffer, map_among_references,0,42);
			}

			


			return 1;

		}

		return 0;



}






inline int try_process_unique_mismatch_end_to_end_output_buffer(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, 
	Output_buffer_sub_block* sub_block, int* get_error,
	Methylation* methy,
	bam_phrase* bam_groups,
	int* map_among_references)
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

		if (output_methy == 1)
		{
			output_methy_end_to_end_output_buffer(name, read, r_read, qulity,
				locates[0],
				read_length - 1,
				0,
				0,
				cigar, read_length, methy, map_among_references);
		}
		else
		{
			output_sam_end_to_end_output_buffer(name, read, r_read, qulity,
				locates[0],
				read_length - 1,
				0,
				0,
				cigar, read_length, sub_block, bam_groups, map_among_references, 0, 42);
		}

		


		return 1;

	}

	return 0;



}













inline int try_process_unique_mismatch_end_to_end_pbat_get_output_buffer(
		bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
		bitmapper_bs_iter* locates,
		bitmapper_bs_iter top, char* tmp_ref, bitmapper_bs_iter* matched_length, int first_C_site, 
		Output_buffer_sub_block* sub_block, int* get_error, Methylation* methy,
		bam_phrase* bam_groups, int* map_among_references)
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

			if (output_methy == 1)
			{
				output_methy_end_to_end_output_buffer(name, read, r_read, qulity,
					locates[0],
					read_length - 1,
					0,
					0,
					cigar, read_length, methy, map_among_references);
			}
			else
			{
				output_sam_end_to_end_pbat_output_buffer(name, read, r_read, qulity,
					locates[0],
					read_length - 1,
					0,
					0,
					cigar, read_length, sub_block, bam_groups, map_among_references,0,42);
			}

			


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
	int inner_minDistance_pair,
	unsigned int* second_best_diff
	)
{
	int mapping_pair = 0;
	///int best_sum_err = 2 * error_threshold + 1;
	int best_sum_err = 4 * error_threshold + 2;
	int inner_i, inner_j;
	long long distance;
	long long current_sum_err;
	long long best_pair_1_index, best_pair_2_index;
	long long first_index_i;
	long long second_index_i;
	///need modification
	///change this value would result in bug
	long long second_best_err = best_sum_err * 2;
	///need modification
	(*second_best_diff) = 0;

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
							if (distance <= maxDistance_pair + current_sum_err
								&&
								distance >= minDistance_pair - current_sum_err)**/
							{


								if (current_sum_err < best_sum_err)
								{
									///need modification
									second_best_err = best_sum_err;

									best_sum_err = current_sum_err;
									best_pair_1_index = inner_i;
									best_pair_2_index = inner_j;
									mapping_pair = 1;
								}
								else if (current_sum_err == best_sum_err)
								{
									///need modification
									second_best_err = best_sum_err;

									mapping_pair++;

									if (best_sum_err == 0)
									{
										///need modification
										(*second_best_diff) = 0;

										(*best_pair_1_index_return) = best_pair_1_index;
										(*best_pair_2_index_return) = best_pair_2_index;
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
							if (distance <= maxDistance_pair + current_sum_err
								&&
								distance >= minDistance_pair - current_sum_err)**/
							{

								if (current_sum_err < best_sum_err)
								{
									///need modification
									second_best_err = best_sum_err;

									best_sum_err = current_sum_err;
									best_pair_1_index = inner_i;
									best_pair_2_index = inner_j;
									mapping_pair = 1;
								}
								else if (current_sum_err == best_sum_err)
								{
									///need modification
									second_best_err = best_sum_err;

									mapping_pair++;

									if (best_sum_err == 0)
									{
										///need modification
										(*second_best_diff) = 0;

										(*best_pair_1_index_return) = best_pair_1_index;
										(*best_pair_2_index_return) = best_pair_2_index;
										return mapping_pair;
									}

								}
							}
						}



					}








				}
			}

		




	}

	///need modification
	if(mapping_pair != 0)
	{
		///need modification
		(*second_best_diff) = second_best_err - best_sum_err;
	}
	


	(*best_pair_1_index_return) = best_pair_1_index;
	(*best_pair_2_index_return) = best_pair_2_index;

	return mapping_pair;
}




inline void debug_new_faster_verify_pairs_second_best_diff(int best_mapp_occ1, int best_mapp_occ2, int error_threshold,
	seed_votes* result1_array, seed_votes* result2_array,
	long long* best_pair_1_index_return, long long* best_pair_2_index_return,
	int inner_maxDistance_pair,
	int inner_minDistance_pair,
	unsigned int second_best_diff
	)
{
	int mapping_pair = 0;
	int best_sum_err = 10 * error_threshold + 2;
	int inner_i, inner_j;
	long long distance;
	long long current_sum_err;
	long long best_pair_1_index, best_pair_2_index;
	///need modification
	///change this value would result in bug
	long long second_best_err = best_sum_err * 2;


	if (best_mapp_occ1 > 0 && best_mapp_occ2 > 0)
	{

		///外层循环是result1_array
		for (inner_i = 0; inner_i < best_mapp_occ1; inner_i++)
		{
			///内层循环是result2_array
			for (inner_j = 0; inner_j < best_mapp_occ2; inner_j++)
			{
				///在内层循环里，result1_array[inner_i].site是不变的
				///变的是result2_array[inner_j].site
				if (result1_array[inner_i].site > result2_array[inner_j].site)
				{
					distance = result1_array[inner_i].site - result2_array[inner_j].site;
				}
				else
				{
					distance = result2_array[inner_j].site - result1_array[inner_i].site;
				}

				if (distance<=inner_maxDistance_pair && distance >= inner_minDistance_pair)
				{
					mapping_pair++;
					current_sum_err = result1_array[inner_i].err + result2_array[inner_j].err;
					if (current_sum_err <= best_sum_err)
					{
						second_best_err = best_sum_err;
						best_sum_err = current_sum_err;
					}

				}
			}
		}
	}
	///no alignment
	if(mapping_pair == 0)
	{
		if(second_best_diff != 0)
		{
			fprintf(stderr, "error 1: second_best_diff: %d\n", second_best_diff);
		}
	}
	else
	{
		int diff = second_best_err - best_sum_err;
		if(diff != second_best_diff)
		{
			if(second_best_diff <= 2*error_threshold || diff <= 2*error_threshold)
			{
				fprintf(stderr, "***error 2: diff: %d, second_best_diff: %d, large_error_threshold: %d\n", 
			diff, second_best_diff, error_threshold);
			}
		}
	}
	
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
#if defined __AVX2__
	__m256i* Peq_SSE,
#else
	__m128i* Peq_SSE,
#endif
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
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
		}
		else
		{
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
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
#if defined __AVX2__
	__m256i* Peq_SSE,
#else
	__m128i* Peq_SSE,
#endif
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
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
		}
		else
		{
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
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
#if defined __AVX2__
	__m256i* Peq_SSE,
#else
	__m128i* Peq_SSE,
#endif
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
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
		}
		else
		{
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
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
#if defined __AVX2__
	__m256i* Peq_SSE,
#else
	__m128i* Peq_SSE,
#endif
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
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
		}
		else
		{
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
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
#if defined __AVX2__
	__m256i* Peq_SSE,
#else
	__m128i* Peq_SSE,
#endif
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
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
		}
		else
		{
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
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
#if defined __AVX2__
	__m256i* Peq_SSE,
#else
	__m128i* Peq_SSE,
#endif
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
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
		}
		else
		{
#if defined __AVX2__
			map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#else
			map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(candidates_votes1, &candidates_votes_length1, (*current_read).length,
				error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ1, (*current_read).name);
#endif
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
#if defined __AVX2__
	__m256i* Peq_SSE)
#else
	__m128i* Peq_SSE)
#endif
{
	unsigned int min_err;
	int best_mapp_occ = 0;

	if (error_threshold <= 15)
	{
#if defined __AVX2__
		map_candidate_votes_mutiple_cut_end_to_end_8_for_paired_end(candidates_votes, &candidates_votes_length, (*current_read).length,
			error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ, (*current_read).name);
#else
		map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end_sse(candidates_votes, &candidates_votes_length, (*current_read).length,
			error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ, (*current_read).name);
#endif
	}
	else
	{
#if defined __AVX2__
		map_candidate_votes_mutiple_cut_end_to_end_4_for_paired_end(candidates_votes, &candidates_votes_length, (*current_read).length,
			error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ, (*current_read).name);
#else
		map_candidate_votes_mutiple_cut_end_to_end_2_for_paired_end_sse(candidates_votes, &candidates_votes_length, (*current_read).length,
			error_threshold, tmp_ref, (*current_read).seq, Peq_SSE, &min_err, &best_mapp_occ, (*current_read).name);
#endif
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

	int my_methylation_size = methylation_size;

	Pair_Methylation methy;

	if (output_methy == 1)
	{
		init_pair_methylation(&methy, my_methylation_size / 2);

		init_pe_distance(&PE_distance[thread_id], minDistance_pair, maxDistance_pair);
	}



	


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






#if defined __AVX2__

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

#else

	__m128i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

#endif











	long long total_number_of_hit = 0;


	fprintf(stderr, "Welcome to BitMapperBS!\n");


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

	bitmapper_bs_iter read1_site;
	bitmapper_bs_iter read1_length;

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
	unsigned int second_best_diff;






	bam_output_cell cell;
	Output_buffer_sub_block current_sub_buffer;

	if (bam_output)
	{
		init_bam_output_cell(&cell);
		init_buffer_sub_block(&current_sub_buffer);
	}

	
	int max_length;



	long long total_bases = 0;
	long long error_bases = 0;



	///这里要改
	bitmapper_bs_iter pre_unique_matched_read = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter pre_ambious_matched_read = (bitmapper_bs_iter)-1;
	int output_mask = 0;

	//正向模式
	i = 0;
	while (1)
	{


		/********************************************输出不匹配结果用的********************************************/
		///这里要改
		///这说明上一个read没有匹配
		///初始pre_ambious_matched_read = (bitmapper_bs_iter)-1, ambious_matched_read = 0
		///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
		if (pre_ambious_matched_read == ambious_matched_read
			&&
			pre_unique_matched_read == unique_matched_read)
		{
			if (unmapped_out == 1)
			{
				directly_output_unmapped_PE(current_read1.name, current_read1.seq, current_read1.rseq, current_read1.qual,
					current_read1.length, &cell, &current_sub_buffer, 1);
				directly_output_unmapped_PE(current_read2.name, current_read2.seq, current_read2.rseq, current_read2.qual,
					current_read2.length, &cell, &current_sub_buffer, 2);
			}
		}
		pre_ambious_matched_read = ambious_matched_read;
		pre_unique_matched_read = unique_matched_read;
		/********************************************输出不匹配结果用的********************************************/




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
		///need verification
		second_best_diff = best_sum_err;


		/*****************************有变化*********************************/
		max_length = current_read1.length;
		if (current_read2.length > max_length)
		{
			max_length = current_read2.length;
		}


		inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
		inner_minDistance_pair = minDistance_pair - large_error_threshold * 2 - max_length;
		/*****************************有变化*********************************/



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
				inner_maxDistance_pair, inner_minDistance_pair, &second_best_diff);



		/**
		debug_new_faster_verify_pairs_second_best_diff(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
				candidates_votes1, candidates_votes2,
				&best_pair_1_index, &best_pair_2_index,
				inner_maxDistance_pair, inner_minDistance_pair,second_best_diff);
		
		if(mapping_pair == 1 && second_best_diff == 0)
		{
			fprintf(stderr, "1 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
		}

		if(mapping_pair > 1 && second_best_diff != 0)
		{
			fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
		}

		if(mapping_pair == 0 && second_best_diff != 0)
		{
			fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
		}
		**/
		


	end_i:
		///这里要改
		if (
			(mapping_pair == 1) 
			||
			(ambiguous_out == 1 &&mapping_pair > 1)
			)
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
					&result1, &matched_length1, 0);



			}
			else
			{
				result1.score = 0;
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
					&result2, &matched_length2, 1);
			}
			else
			{
				result2.score = 0;
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

			paired_end_distance = calculate_TLEN(result1.site, matched_length1, result2.site, matched_length2);

			if ((paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
				&&
				(result1.site + matched_length1 <=_ih_refGenName[result1.chrome_id]._rg_chrome_length + 1)
				&&
				(result2.site + matched_length2 <=_ih_refGenName[result2.chrome_id]._rg_chrome_length + 1))
			{

				///这里要改
				if (mapping_pair == 1)
				{
					unique_matched_read++;
				}
				else
				{
					ambious_matched_read++;
				}

				

				total_bases = total_bases + current_read1.length + current_read2.length;
				error_bases = error_bases + result1.err + result2.err;

				///need modification
				int mapq = MAP_Calculation(second_best_diff, error_threshold1 + error_threshold2, 
				result1.score + result2.score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
				/**
				fprintf(stderr, "mapq: %d, second_best_diff: %d, error_threshold1: %d, error_threshold2: %d\n", 
				mapq, second_best_diff, error_threshold1, error_threshold2);
				fprintf(stderr, "result1.score: %d, result2.score: %d, GapOpenPenalty: %d, GapExtensionPenalty: %d, MistMatchPenaltyMax: %d\n", 
				result1.score, result2.score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
				**/

				if (output_methy == 1)
				{

					output_paired_methy_buffer(current_read1.name,
						current_read1.seq, current_read1.rseq, &result1, matched_length1,
						current_read2.seq, current_read2.rseq, &result2, matched_length2,
						&methy, &PE_distance[thread_id]);
				
				}
				else
				{
					/**
					///这里要改
					if (mapping_pair == 1)
					{
						output_mask = 0;
					}
					else
					{
						output_mask = 0x100;
					}
					**/

					directly_output_read1(current_read1.name, current_read1.seq, current_read1.rseq, current_read1.qual,
						&result1, &result2, current_read1.length,
						matched_length1, matched_length2, paired_end_distance,
						&cell,
						&current_sub_buffer, output_mask, mapq);

					directly_output_read2(current_read2.name, current_read2.seq, current_read2.rseq, current_read2.qual,
						&result2, &result1, current_read2.length,
						matched_length2, matched_length1, paired_end_distance,
						&cell,
						&current_sub_buffer, output_mask, mapq);

				}

			}






		}
		else if (mapping_pair > 1)
		{
			ambious_matched_read++;
		}

	}


	if (output_methy == 1 && methy.current_size != 0)
	{
		/**
		assign_cuts(&methy);
		output_methylation(&methy);
		clear_methylation(&methy);
		**/
		assign_cuts_pair(&methy);
		output_methylation_pair(&methy);
		clear_methylation_pair(&methy);

	}


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;


	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;


	/**
	fprintf(stderr, "debug_1: %lld\n", debug_1);
	fprintf(stderr, "debug_2: %lld\n", debug_2);
	**/

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
	/**
	FILE* debug_f = fopen("sb_debug.txt", "w");
	**/

	int my_methylation_size = methylation_size;

	Pair_Methylation methy;
	/**
	fprintf(debug_f, "methy.current_size: %d\n", methy.current_size);
	fflush(debug_f);

	fprintf(debug_f, "output_methy: %d\n", output_methy);
	fflush(debug_f);
	**/

	if (output_methy == 1)
	{
		init_pair_methylation(&methy, my_methylation_size / 2);

		init_pe_distance(&PE_distance[thread_id], minDistance_pair, maxDistance_pair);
	}



	


	bitmapper_bs_iter read1_site;
	bitmapper_bs_iter read1_length;


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


#if defined __AVX2__

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

#else

	__m128i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

#endif



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






	bam_output_cell cell;
	Output_buffer_sub_block current_sub_buffer;

	if (bam_output)
	{
		init_bam_output_cell(&cell);
		init_buffer_sub_block(&current_sub_buffer);
	}


	int max_length;


	long long total_bases = 0;
	long long error_bases = 0;


	///这里要改
	bitmapper_bs_iter pre_unique_matched_read = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter pre_ambious_matched_read = (bitmapper_bs_iter)-1;
	int output_mask = 0;
	unsigned int second_best_diff;

	

	//正向模式
	i = 0;
	while (1)
	{



		/********************************************输出不匹配结果用的********************************************/
		///这里要改
		///这说明上一个read没有匹配
		///初始pre_ambious_matched_read = (bitmapper_bs_iter)-1, ambious_matched_read = 0
		///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
		if (pre_ambious_matched_read == ambious_matched_read
			&&
			pre_unique_matched_read == unique_matched_read)
		{
			if (unmapped_out == 1)
			{
				directly_output_unmapped_PE(current_read1.name, current_read1.seq, current_read1.rseq, current_read1.qual,
					current_read1.length, &cell, &current_sub_buffer, 1);
				directly_output_unmapped_PE(current_read2.name, current_read2.seq, current_read2.rseq, current_read2.qual,
					current_read2.length, &cell, &current_sub_buffer, 2);
			}
		}
		pre_ambious_matched_read = ambious_matched_read;
		pre_unique_matched_read = unique_matched_read;
		/********************************************输出不匹配结果用的********************************************/




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
		///need verification
		second_best_diff = best_sum_err;


		/*****************************有变化*********************************/
		max_length = current_read1.length;
		if (current_read2.length > max_length)
		{
			max_length = current_read2.length;
		}


		inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
		inner_minDistance_pair = minDistance_pair - large_error_threshold * 2 - max_length;
		/*****************************有变化*********************************/





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
		///in current version, unique_best_read1 and unique_best_read2 are always -1
		///so the following codes are not useful
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
				inner_maxDistance_pair, inner_minDistance_pair, &second_best_diff);
		}
		else
		{
			///second_best_diff > large_error_threshold means there is only one alignment
			second_best_diff = large_error_threshold * 2 + 1;
		}



		/**
		debug_new_faster_verify_pairs_second_best_diff(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
				candidates_votes1, candidates_votes2,
				&best_pair_1_index, &best_pair_2_index,
				inner_maxDistance_pair, inner_minDistance_pair,second_best_diff);
		
		if(mapping_pair == 1 && second_best_diff == 0)
		{
			fprintf(stderr, "1 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
		}

		if(mapping_pair > 1 && second_best_diff != 0)
		{
			fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
		}

		if(mapping_pair == 0 && second_best_diff != 0)
		{
			fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
		}
		**/
		
		

	end_i:
		///这里要改
		if (
			(mapping_pair == 1)
			||
			(ambiguous_out == 1 && mapping_pair > 1)
			)
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
					&result1, &matched_length1, 0);

				

			}
			else
			{
				result1.score = 0;
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
					&result2, &matched_length2, 1);
			}
			else
			{
				result2.score = 0;
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

			paired_end_distance = calculate_TLEN(result1.site, matched_length1, result2.site, matched_length2);


			///if (paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
		if ((paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
			&&
			(result1.site + matched_length1 <= _ih_refGenName[result1.chrome_id]._rg_chrome_length + 1)
			&&
			(result2.site + matched_length2 <= _ih_refGenName[result2.chrome_id]._rg_chrome_length + 1))
			{

				///这里要改
				if (mapping_pair == 1)
				{
					unique_matched_read++;
				}
				else
				{
					ambious_matched_read++;
				}


				total_bases = total_bases + current_read1.length + current_read2.length;
				error_bases = error_bases + result1.err + result2.err;

				///need modification
				int mapq = MAP_Calculation(second_best_diff, error_threshold1 + error_threshold2, 
				result1.score + result2.score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);

				if (output_methy == 1)
				{

					output_paired_methy_buffer(current_read1.name,
						current_read1.seq, current_read1.rseq, &result1, matched_length1,
						current_read2.seq, current_read2.rseq, &result2, matched_length2,
						&methy, &PE_distance[thread_id]);
				}
				else
				{

					/**
					///这里要改
					if (mapping_pair == 1)
					{
						output_mask = 0;
					}
					else
					{
						output_mask = 0x100;
					}
					**/


					directly_output_read1(current_read1.name, current_read1.seq, current_read1.rseq, current_read1.qual,
						&result1, &result2, current_read1.length,
						matched_length1, matched_length2, paired_end_distance,
						&cell,
						&current_sub_buffer, output_mask, mapq);

					directly_output_read2(current_read2.name, current_read2.seq, current_read2.rseq, current_read2.qual,
						&result2, &result1, current_read2.length,
						matched_length2, matched_length1, paired_end_distance,
						&cell,
						&current_sub_buffer, output_mask, mapq);
				}
				
				
				
			}

			
			


		}
		else if (mapping_pair > 1)
		{
			ambious_matched_read++;
		}

	

		i++;
	}

	/**
	fprintf(debug_f, "1********************************\n");
	fflush(debug_f);

	fprintf(debug_f, "methy.current_size: %d\n", methy.current_size);
	fflush(debug_f);
	**/


	if (output_methy == 1 && methy.current_size != 0)
	{
		assign_cuts_pair(&methy);
		output_methylation_pair(&methy);
		clear_methylation_pair(&methy);

	}
	/**
	fprintf(debug_f, "2********************************\n");
	fflush(debug_f);
	**/

	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;
	/**
	fprintf(debug_f, "3********************************\n");
	fflush(debug_f);
	**/

	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;
	/**
	fprintf(debug_f, "4********************************\n");
	fflush(debug_f);
	**/
	/**
	fprintf(stderr, "debug_1: %lld\n", debug_1);
	fprintf(stderr, "debug_2: %lld\n", debug_2);
	**/

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

	bitmapper_bs_iter read1_site;
	bitmapper_bs_iter read1_length;


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




#if defined __AVX2__

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

#else

	__m128i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

#endif

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


	Pair_Methylation methy;
	unsigned int second_best_diff;
	



	if (output_methy == 1)
	{
		int my_methylation_size = methylation_size;
		
		if (output_methy == 1)
		{
			///注意这里要除以2
			init_pair_methylation(&methy, my_methylation_size / 2);
		}


		init_pe_distance(&PE_distance[thread_id], minDistance_pair, maxDistance_pair);
	}
	else
	{
		init_buffer_sub_block(&current_sub_buffer);
	}

	
	bam_phrase bam_buffer;

	if (bam_output == 1)
	{
		init_bam_buffer(&bam_buffer);
	}



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


	int max_length;



	long long total_bases = 0;
	long long error_bases = 0;

	///这里要改
	bitmapper_bs_iter pre_unique_matched_read = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter pre_ambious_matched_read = (bitmapper_bs_iter)-1;
	int output_mask = 0;


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

		////这里要改
		pre_unique_matched_read = (bitmapper_bs_iter)-1;
		pre_ambious_matched_read = (bitmapper_bs_iter)-1;

		////这里要改
		while (1)
		{


			/********************************************输出不匹配结果用的********************************************/
			///这里要改
			///这说明上一个read没有匹配
			///初始pre_ambious_matched_read = (bitmapper_bs_iter)-1, ambious_matched_read = 0
			///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
			if (pre_ambious_matched_read == ambious_matched_read
				&&
				pre_unique_matched_read == unique_matched_read)
			{
				if (unmapped_out == 1)
				{
					directly_output_read1_return_buffer_unmapped_PE(
						read_batch1[i - 1].name, read_batch1[i - 1].seq, read_batch1[i - 1].rseq, read_batch1[i - 1].qual,
						read_batch1[i - 1].length, &current_sub_buffer, &bam_buffer, 1);
					directly_output_read1_return_buffer_unmapped_PE(
						read_batch2[i - 1].name, read_batch2[i - 1].seq, read_batch2[i - 1].rseq, read_batch2[i - 1].qual,
						read_batch2[i - 1].length, &current_sub_buffer, &bam_buffer, 2);
				}
			}
			pre_ambious_matched_read = ambious_matched_read;
			pre_unique_matched_read = unique_matched_read;
			/********************************************输出不匹配结果用的********************************************/

			////这里要改
			if (i >= obtain_reads_num)
			{
				break;
			}







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
			///need verification
			second_best_diff = best_sum_err;


			/*****************************有变化*********************************/
			max_length = read_batch1[i].length;
			if (read_batch2[i].length > max_length)
			{
				max_length = read_batch2[i].length;
			}


			inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
			inner_minDistance_pair = minDistance_pair - large_error_threshold * 2 - max_length;
			/*****************************有变化*********************************/


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
					inner_maxDistance_pair, inner_minDistance_pair, &second_best_diff);

			
			/**
			debug_new_faster_verify_pairs_second_best_diff(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
				candidates_votes1, candidates_votes2,
				&best_pair_1_index, &best_pair_2_index,
				inner_maxDistance_pair, inner_minDistance_pair,second_best_diff);
		
			if(mapping_pair == 1 && second_best_diff == 0)
			{
				fprintf(stderr, "1 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
			}

			if(mapping_pair > 1 && second_best_diff != 0)
			{
				fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
			}

			if(mapping_pair == 0 && second_best_diff != 0)
			{
				fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
			}
			**/



		end_i:
			///这里要改
			if (
				(mapping_pair == 1)
				||
				(ambiguous_out == 1 && mapping_pair > 1)
				)
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
						&result1, &matched_length1, 0);



				}
				else
				{
					result1.score = 0;
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
						&result2, &matched_length2, 1);
				}
				else
				{
					result2.score = 0;
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


				/**
				if (result2.site >= result1.site)
				{
					paired_end_distance = result2.site - result1.site + matched_length2;

				}
				else
				{
					paired_end_distance = result1.site - result2.site + matched_length1;
				}
				**/

				paired_end_distance = calculate_TLEN(result1.site, matched_length1, result2.site, matched_length2);


				///if (paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
				if ((paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
					&&
					(result1.site + matched_length1 <= _ih_refGenName[result1.chrome_id]._rg_chrome_length + 1)
					&&
					(result2.site + matched_length2 <= _ih_refGenName[result2.chrome_id]._rg_chrome_length + 1))
				{
					///这里要改
					if (mapping_pair == 1)
					{
						unique_matched_read++;
					}
					else
					{
						ambious_matched_read++;
					}


					total_bases = total_bases + read_batch1[i].length + read_batch2[i].length;
					error_bases = error_bases + result1.err + result2.err;

					///need modification
					int mapq = MAP_Calculation(second_best_diff, error_threshold1 + error_threshold2, 
					result1.score + result2.score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
					
					if (output_methy == 1)
					{
						output_paired_methy_buffer_mutiple_threads(
							read_batch1[i].seq, read_batch1[i].rseq, &result1, matched_length1,
							read_batch2[i].seq, read_batch2[i].rseq, &result2, matched_length2,
							&methy, &PE_distance[thread_id]);

					}
					else
					{
						/**
						///这里要改
						if (mapping_pair == 1)
						{
							output_mask = 0;
						}
						else
						{
							output_mask = 0x100;
						}
						**/

						directly_output_read1_return_buffer(read_batch1[i].name, read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual,
							&result1, &result2, read_batch1[i].length,
							matched_length1, matched_length2, &current_sub_buffer, paired_end_distance, &bam_buffer, output_mask, mapq);

						directly_output_read2_return_buffer(read_batch2[i].name, read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual,
							&result2, &result1, read_batch2[i].length,
							matched_length2, matched_length1, &current_sub_buffer, paired_end_distance, &bam_buffer, output_mask, mapq);
					}

					

				}

				

			}
			else if (mapping_pair > 1)
			{
				ambious_matched_read++;

				
			}


			i++;
		}

		///if (output_methy == 0)
		if (output_methy == 0 && bam_output == 0)
		{

			current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

			push_results_to_buffer(&current_sub_buffer);

			current_sub_buffer.length = 0;
		}

	}

	if (output_methy == 1 && methy.current_size != 0)
	{
		assign_cuts_pair(&methy);

		push_methy_to_buffer_pair(&methy);

		///这个不用清理，因为再push 的时候清理过了
		///clear_methylation(&methy);

	}



	if (bam_output == 1)
	{

		flush_bam_buffer(&bam_buffer);
	}


	if (bam_output == 1)
	{
		finish_bam_output_buffer();

	}
	else
	{
		finish_output_buffer_pair();
	}

	


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;


	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;

	fprintf(stderr, "thread %d completed!\n", thread_id);


}



///输出成对匹配
///和一端unique一端不匹配的结果
void* Map_Pair_Seq_split(void* arg)
{

	int thread_id = *((int *)arg);

	FILE* output_file = get_Ouput_Dec();

	bwt_locate_queue get_queue;

	init_locate_queue_muti_thread(&get_queue);



	bitmapper_bs_iter read1_site;
	bitmapper_bs_iter read1_length;


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



#if defined __AVX2__

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

#else

	__m128i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

#endif





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




	Read* read_batch1;
	Read* read_batch2;

	Read_buffer_pe_sub_block curr_sub_block;

	init_single_sub_block_pe(&curr_sub_block);

	Output_buffer_sub_block current_sub_buffer;
	

	Pair_Methylation methy;



	if (output_methy == 1)
	{
		int my_methylation_size = methylation_size;

		if (output_methy == 1)
		{
			///注意这里要除以2
			init_pair_methylation(&methy, my_methylation_size / 2);
		}


		init_pe_distance(&PE_distance[thread_id], minDistance_pair, maxDistance_pair);
	}
	else
	{
		init_buffer_sub_block(&current_sub_buffer);
	}



	bam_phrase bam_buffer;

	if (bam_output == 1)
	{
		init_bam_buffer(&bam_buffer);
	}


	


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





	long long total_bases = 0;
	long long error_bases = 0;


	///这里要改
	bitmapper_bs_iter pre_unique_matched_read = (bitmapper_bs_iter)-1;
	bitmapper_bs_iter pre_ambious_matched_read = (bitmapper_bs_iter)-1;
	int output_mask = 0;
	unsigned int second_best_diff;


	int max_length;

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

		////这里要改
		pre_unique_matched_read = (bitmapper_bs_iter)-1;
		pre_ambious_matched_read = (bitmapper_bs_iter)-1;


		////这里要改
		while (1)
		{




			/********************************************输出不匹配结果用的********************************************/
			///这里要改
			///这说明上一个read没有匹配
			///初始pre_ambious_matched_read = (bitmapper_bs_iter)-1, ambious_matched_read = 0
			///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
			if (pre_ambious_matched_read == ambious_matched_read
				&&
				pre_unique_matched_read == unique_matched_read)
			{
				if (unmapped_out == 1)
				{
					directly_output_read1_return_buffer_unmapped_PE(
						read_batch1[i - 1].name, read_batch1[i - 1].seq, read_batch1[i - 1].rseq, read_batch1[i - 1].qual,
						read_batch1[i - 1].length, &current_sub_buffer, &bam_buffer, 1);
					directly_output_read1_return_buffer_unmapped_PE(
						read_batch2[i - 1].name, read_batch2[i - 1].seq, read_batch2[i - 1].rseq, read_batch2[i - 1].qual,
						read_batch2[i - 1].length, &current_sub_buffer, &bam_buffer, 2);
				}
			}
			pre_ambious_matched_read = ambious_matched_read;
			pre_unique_matched_read = unique_matched_read;
			/********************************************输出不匹配结果用的********************************************/




			////这里要改
			if (i >= obtain_reads_num)
			{
				break;
			}











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
			///need verification
			second_best_diff = best_sum_err;

			

			


			/*****************************有变化*********************************/
			max_length = read_batch1[i].length;
			if (read_batch2[i].length > max_length)
			{
				max_length = read_batch2[i].length;
			}


			inner_maxDistance_pair = maxDistance_pair + large_error_threshold * 2;
			inner_minDistance_pair = minDistance_pair - large_error_threshold * 2 - max_length;
			/*****************************有变化*********************************/






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
			///in current version, unique_best_read1 and unique_best_read2 are always -1
			///so the following codes are not useful
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
					inner_maxDistance_pair, inner_minDistance_pair, &second_best_diff);
			}
			else
			{
				///second_best_diff > large_error_threshold means there is only one alignment
				second_best_diff = large_error_threshold * 2 + 1;
			}

			
			/**
			debug_new_faster_verify_pairs_second_best_diff(best_mapp_occ1, best_mapp_occ2, large_error_threshold,
				candidates_votes1, candidates_votes2,
				&best_pair_1_index, &best_pair_2_index,
				inner_maxDistance_pair, inner_minDistance_pair,second_best_diff);
		
			if(mapping_pair == 1 && second_best_diff == 0)
			{
				fprintf(stderr, "1 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
			}

			if(mapping_pair > 1 && second_best_diff != 0)
			{
				fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
			}

			if(mapping_pair == 0 && second_best_diff != 0)
			{
				fprintf(stderr, "2 second_best_diff: %d, large_error_threshold: %d\n", second_best_diff, large_error_threshold);
			}
			**/


		end_i:
			///这里要改
			if (
				(mapping_pair == 1)
				||
				(ambiguous_out == 1 && mapping_pair > 1)
				)
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
						&result1, &matched_length1, 0);



				}
				else
				{
					result1.score = 0;
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
						&result2, &matched_length2, 1);
				}
				else
				{
					result2.score = 0;
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


				paired_end_distance = calculate_TLEN(result1.site, matched_length1, result2.site, matched_length2);


				///if (paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
				if ((paired_end_distance <= maxDistance_pair && paired_end_distance >= minDistance_pair)
					&&
					(result1.site + matched_length1 <= _ih_refGenName[result1.chrome_id]._rg_chrome_length + 1)
					&&
					(result2.site + matched_length2 <= _ih_refGenName[result2.chrome_id]._rg_chrome_length + 1))
				{
					///这里要改
					if (mapping_pair == 1)
					{
						unique_matched_read++;
					}
					else
					{
						ambious_matched_read++;
					}


					total_bases = total_bases + read_batch1[i].length + read_batch2[i].length;
					error_bases = error_bases + result1.err + result2.err;

					///need modification
					int mapq = MAP_Calculation(second_best_diff, error_threshold1 + error_threshold2, 
					result1.score + result2.score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);


					if (output_methy == 1)
					{

						output_paired_methy_buffer_mutiple_threads(
							read_batch1[i].seq, read_batch1[i].rseq, &result1, matched_length1,
							read_batch2[i].seq, read_batch2[i].rseq, &result2, matched_length2,
							&methy, &PE_distance[thread_id]);
					}
					else
					{
						/**
						///这里要改
						if (mapping_pair == 1)
						{
							output_mask = 0;
						}
						else
						{
							output_mask = 0x100;
						}
						**/


						directly_output_read1_return_buffer(read_batch1[i].name, read_batch1[i].seq, read_batch1[i].rseq, read_batch1[i].qual,
							&result1, &result2, read_batch1[i].length,
							matched_length1, matched_length2, &current_sub_buffer, paired_end_distance, &bam_buffer, output_mask, mapq);

						directly_output_read2_return_buffer(read_batch2[i].name, read_batch2[i].seq, read_batch2[i].rseq, read_batch2[i].qual,
							&result2, &result1, read_batch2[i].length,
							matched_length2, matched_length1, &current_sub_buffer, paired_end_distance, &bam_buffer, output_mask, mapq);
					}

					

				}



			}
			else if (mapping_pair > 1)
			{
				ambious_matched_read++;
			}


			i++;
		}


		if (output_methy == 0 && bam_output == 0)
		{
			current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

			push_results_to_buffer(&current_sub_buffer);

			current_sub_buffer.length = 0;
		}

	}



	if (output_methy == 1 && methy.current_size != 0)
	{
		assign_cuts_pair(&methy);

		push_methy_to_buffer_pair(&methy);

		///这个不用清理，因为再push 的时候清理过了
		///clear_methylation(&methy);

	}


	if (bam_output == 1)
	{
		flush_bam_buffer(&bam_buffer);
	}


	if (bam_output == 1)
	{
		finish_bam_output_buffer();

	}
	else
	{
		finish_output_buffer_pair();
	}
	

	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = ambious_matched_read;

	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;

	fprintf(stderr, "thread %d completed!\n", thread_id);


}


inline void output_sam_unmapped(char* name, char* read, char* r_read, char* qulity, int read_length,
	bam_output_cell* cell, Output_buffer_sub_block* sub_block)
{

	if (bam_output == 0)
	{



		if (name[0] == '@')
		{
			fprintf(schema_out_fp, "%s\t", name + 1);
		}
		else
		{
			fprintf(schema_out_fp, "%s\t", name);
		}



		fprintf(schema_out_fp, "4\t*\t0\t0\t*\t*\t0\t0\t");


		fprintf(schema_out_fp, "%s\t", read);

		fprintf(schema_out_fp, "%s\n", qulity);

	}
	else
	{


		///这个一定一定要有
		sub_block->length = 0;




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


		///fprintf(schema_out_fp, "4\t*\t0\t0\t*\t*\t0\t0\t");
		sub_block->buffer[sub_block->length] = '4';
		sub_block->length++;
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
		sub_block->buffer[sub_block->length] = '*';
		sub_block->length++;
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



		///result_string << read << "\t";
		output_to_buffer_char_length(sub_block, read, read_length);
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\t';
		sub_block->length++;

		///result_string << qulity << "\t";
		output_to_buffer_char_length(sub_block, qulity, read_length);


		///注意这里和多线程原始的算法不一样
		/**
		//上面扩容时有32的余量，所以下面\t不需要检查了
		sub_block->buffer[sub_block->length] = '\n';
		sub_block->length++;
		**/
		sub_block->buffer[sub_block->length] = '\0';

		write_alignment_directly(sub_block->buffer, sub_block->length, cell);

	}


}



inline int output_ambiguous_exact_map(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top,
	bitmapper_bs_iter number_of_hits,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter* matched_length,
	bam_output_cell* cell,
	Output_buffer_sub_block* current_sub_buffer,
	int* map_among_references)
{
	bitmapper_bs_iter tmp_SA_length = 0;

	///限制遍历上限
	if (number_of_hits > max_seed_matches)
	{
		number_of_hits = max_seed_matches;
	}

	///这个cigar一定是精确匹配的
	sprintf(cigar, "%dM", read_length);


	bitmapper_bs_iter i = 0;

	for (i = 0; i < number_of_hits; i++)
	{

		tmp_SA_length = 0;

		locate_one_position_direct(locates, top + i, &tmp_SA_length, total_SA_length, *matched_length, 0);

		///in this case, mapq should be 1
		output_sam_end_to_end(name, read, r_read, qulity,
			locates[0],
			read_length - 1,
			0,
			0,
			cigar, read_length,
			cell,
			current_sub_buffer,
			map_among_references, 0, 1);

		///如果找到了一个匹配，就结束
		if ((*map_among_references) == 0)
		{
			return 1;
		}
	}
}



inline int debug_1_mismatch_score(bitmapper_bs_iter map_location, char* ref, char* read, int readLen, 
	char* qulity, int quality_base,
	int MistMatchPenaltyMax, int MistMatchPenaltyMin, int N_Penalty, int need_r_quality, int mismatch_site, int check_score)
{
	int score = 0;
	int error = 0;
	int i = 0;

	if(map_location < _msf_refGenLength)
	{
		get_actuall_genome(ref, map_location, readLen);
	}
	else
	{
		get_actuall_rc_genome(ref, map_location - _msf_refGenLength, readLen);
	}

	for (i = 0; i < readLen; i++)
	{
		if (ref[i] != read[i])
		{
			if (!(read[i] == 'T' && ref[i] == 'C'))
			{
				if(i != mismatch_site)
				{
					fprintf(stderr, "mismatch_site: %d, i: %d\n", mismatch_site, i);
				}

				error++;

				if(read[i] == 'N' || ref[i] == 'N')
				{
					score -= N_Penalty;
				}
				else
				{
					if(!need_r_quality)
					{
						score -= 
						MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, quality_base, qulity[i]);
						///fprintf(stderr, "i: %d, qulity[i]: %d, penalty: %d\n", i, qulity[i], MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, quality_base, qulity[i]));
					}
					else
					{
						score -= 
						MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, quality_base, qulity[readLen - i - 1]);
					}
				}
			}


			
		}

	}


	if(error != 1)
	{
		fprintf(stderr, "(*error): %d, (*score): %d\n", error, score);
	}

	if(check_score != score)
	{
		fprintf(stderr, "check_score: %d, score: %d\n", check_score, score);
	}


}



				





int Map_Single_Seq_end_to_end(int thread_id)
{
	///int my_methylation_size = methylation_size / THREAD_COUNT;
	int my_methylation_size = methylation_size;

	Methylation methy;
	

	if (output_methy == 1)
	{
		init_methylation(&methy, my_methylation_size);
	}



	

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

	fprintf(stderr, "Welcome to BitMapperBS!\n");

	#if defined __AVX2__

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

	#else

	__m128i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

	#endif

	


	int is_mutiple_map = 0;

	int map_among_references;

	int C_site;


	int return_flag;

	double startTime = Get_T();

	int file_flag;

	int get_error;

	//long long debug_0_mismatch = 0;
	///long long debug_1_mismatch = 0;
	//long long debug_seed_match_length = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;
	///int processed = 0;









	bam_output_cell cell;
	Output_buffer_sub_block current_sub_buffer;

	if (bam_output)
	{
		init_bam_output_cell(&cell);
		init_buffer_sub_block(&current_sub_buffer);
	}
	






	long long total_bases = 0;
	long long error_bases = 0;





	///这里要改
	bitmapper_bs_iter pre_matched_read = (bitmapper_bs_iter)-1;
	int output_mask;
	long long ambiguous_index;


	///need modification
	unsigned int second_best_diff;
	int one_mismatch_site;
	/**
	kswr_t qry;
	init_qry_total(&qry);
	**/
	


	//正向模式
	i = 0;
	while (1)
	{


		///这里要改
		/********************************************输出不匹配结果用的********************************************/
		///这里要改
		///这说明上一个read没有匹配
		///初始pre_matched_read = (bitmapper_bs_iter)-1, matched_read = 0
		///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
		if (pre_matched_read == matched_read)
		{
			/**
			if(second_best_diff != 0 && map_among_references == 0)
			{
				fprintf(stderr, "2: second_best_diff: %d, candidate_length: %d, name: %s, current_read.length: %d\n", 
				second_best_diff, candidate_length, current_read.name, current_read.length);
			}
			**/
			
			if (unmapped_out == 1)
			{
				output_sam_unmapped(current_read.name, current_read.seq, current_read.rseq, current_read.qual, 
					current_read.length, &cell, &current_sub_buffer);
			}
		}
		pre_matched_read = matched_read;
		/********************************************输出不匹配结果用的********************************************/








		file_flag = inputReads_single_directly(&current_read);
		///如果等于0, 估计说明read读完了吧
		if (file_flag == 0)
		{
			break;
		}


		enq_i++;


		///等于3应该是有过多的N吧
		///这个功能应该是不存在的
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

		map_among_references = 0;

		second_best_diff = 0;

		


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

			
			if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end(
				current_read.length, current_read.seq, current_read.rseq, current_read.qual, cigar, current_read.name,
				locates,
				top, tmp_ref, &match_length, current_read.length - C_site - 1, &get_error,
				&methy,
				&cell,
				&current_sub_buffer,
				&map_among_references))
			{

				if (map_among_references == 0)
				{
					unique_matched_read++;
					matched_read++;

					total_bases = total_bases + current_read.length;

				}

				
				continue;
			}

			one_mismatch_site = match_length;


			///只有第一个seed才有可能出现这种情况
			if (match_length == current_read.length
				&&
				number_of_hits > 1)
			{
				is_mutiple_map = 1;

				///等于-1说明read里面没有C
				if (C_site == -1)
				{
					i++;
					matched_read++;

					///这里要改
					if (ambiguous_out == 1)
					{

						output_ambiguous_exact_map(
							current_read.length,
							current_read.seq,
							current_read.rseq,
							current_read.qual,
							cigar,
							current_read.name,
							locates,
							top,
							number_of_hits,
							max_seed_matches,
							&match_length,
							&cell,
							&current_sub_buffer,
							&map_among_references);

						if (map_among_references != 0)
						{
							matched_read--;
						}

					}
					
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




		
		///两种特殊情况为error为0和error为1，这两种情况在上面特殊处理就好了
		///除了这两种情况, 均在下面处理
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

		///这也是个特殊情况, 比较复杂....
		///candidate_length == 2 && candidates[0] == candidates[1] 这是两个种子都指向了同一个位置
		///1-mismatch
		if (extra_seed_flag == 0 
			&& 
			(candidate_length == 1 || 
			(candidate_length == 2 && candidates[0] == candidates[1])
			))
		{
			sprintf(cigar, "%dM", current_read.length);


			if (output_methy == 1)
			{
				output_methy_end_to_end(
					current_read.name,
					current_read.seq,
					current_read.rseq,
					current_read.qual,
					candidates[0],
					current_read.length - 1,
					0,
					1,
					cigar,
					current_read.length,
					&methy,
					&map_among_references);
			}
			else
			{

				///need modification
				int score = 0;

				if(current_read.seq[one_mismatch_site] == 'N')
				{
					score -= N_Penalty;
				}
				else
				{
					score -= MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, Q_base, current_read.qual[one_mismatch_site]);
				}

				int mapq = MAP_Calculation((unsigned int)-1, error_threshold1, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
				/**
				debug_1_mismatch_score(candidates[0], tmp_ref, current_read.seq, current_read.length, 
					current_read.qual, Q_base, MistMatchPenaltyMax, MistMatchPenaltyMin, N_Penalty, 0, one_mismatch_site, score);
				**/
				
				

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
					current_read.length,
					&cell,
					&current_sub_buffer,
					&map_among_references, 0, mapq);
			}




			if (map_among_references == 0)
			{

				unique_matched_read++;
				matched_read++;

				total_bases = total_bases + current_read.length;
				error_bases = error_bases + 1;
				
			}

			///debug_1_mismatch++;
		}
		else if (candidate_length != 0)   ///这里才是正常情况
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
				#if defined __AVX2__
					map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#else
					map_candidate_votes_mutiple_cut_end_to_end_4_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#endif
					
				}
				else
				{
				#if defined __AVX2__
					map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#else
					map_candidate_votes_mutiple_cut_end_to_end_2_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#endif
				}
					
				
			}
			else
			{

				if (error_threshold1 <= 15)
				{
				#if defined __AVX2__
					map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#else
					map_candidate_votes_mutiple_end_to_end_4_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#endif
					
				}
				else
				{
				#if defined __AVX2__
					map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#else
					map_candidate_votes_mutiple_end_to_end_2_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
				#endif
				}
				
				
			}


			if (min_err_index >= 0)
			{
				/**
				if(second_best_diff <= 0)
				{
					fprintf(stderr, "0: second_best_diff: %d\n", second_best_diff);
				}
				**/
				
				
				calculate_best_map_cigar_end_to_end
					(candidates_votes, &min_candidates_votes_length,
					current_read.length, error_threshold1, tmp_ref,
					current_read.seq, current_read.rseq, current_read.qual, cigar, path, matrix, matrix_bit, current_read.name,
					best_cigar, &methy,
					&cell,
					&current_sub_buffer,
					&map_among_references,0,second_best_diff);

				///fprintf(stderr, "map_among_references: %d\n", map_among_references);


				if (map_among_references == 0)
				{
					matched_read++;
					unique_matched_read++;

					total_bases = total_bases + current_read.length;
					error_bases = error_bases + candidates_votes[0].err;
				}


				
			}
			else if (min_err_index != -1)
			{
				/**
				if(second_best_diff != 0)
				{
					fprintf(stderr, "1: second_best_diff: %d\n", second_best_diff);
				}
				**/
				

				matched_read++;

				///这里要改
				if (ambiguous_out == 1)
				{
					ambiguous_index = -2 - min_err_index;
					candidates_votes[0].err = candidates_votes[ambiguous_index].err;
					candidates_votes[0].end_site = candidates_votes[ambiguous_index].end_site;
					candidates_votes[0].site = candidates_votes[ambiguous_index].site;

					calculate_best_map_cigar_end_to_end
						(candidates_votes, &min_candidates_votes_length,
						current_read.length, error_threshold1, tmp_ref,
						current_read.seq, current_read.rseq, current_read.qual, cigar, path, matrix, matrix_bit, current_read.name,
						best_cigar, &methy,
						&cell,
						&current_sub_buffer,
						&map_among_references, 0,second_best_diff);

					if (map_among_references != 0)
					{
						matched_read--;
					}
				}
			}

		}
	}

	
	if (output_methy == 1 && methy.current_size != 0)
	{
		assign_cuts(&methy);
		output_methylation(&methy);
		clear_methylation(&methy);

	}

	/**
	for (i = 0; i < methy_out.genome_cut; i++)
	{
		fprintf(stderr, "i: %d, %llu\n", i, methy_out.cut_index[i]);
	}
	**/
	

	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;

	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;

	/**
	fprintf(stderr, "debug_0_mismatch: %lld\n", debug_0_mismatch);
	fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);
	fprintf(stderr, "debug_seed_match_length: %lld\n", debug_seed_match_length);
	**/
	
	///fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);

	return 1;
}



inline int output_ambiguous_exact_map_pbat(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top,
	bitmapper_bs_iter number_of_hits,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter* matched_length,
	bam_output_cell* cell,
	Output_buffer_sub_block* current_sub_buffer,
	int* map_among_references)
{
	bitmapper_bs_iter tmp_SA_length = 0;

	///限制遍历上限
	if (number_of_hits > max_seed_matches)
	{
		number_of_hits = max_seed_matches;
	}

	///这个cigar一定是精确匹配的
	sprintf(cigar, "%dM", read_length);


	bitmapper_bs_iter i = 0;

	for (i = 0; i < number_of_hits; i++)
	{

		tmp_SA_length = 0;

		locate_one_position_direct(locates, top + i, &tmp_SA_length, total_SA_length, *matched_length, 0);

		///in this case, mapq should be 1
		output_sam_end_to_end_pbat(name, read, r_read, qulity,
			locates[0],
			read_length - 1,
			0,
			0,
			cigar, read_length,
			cell,
			current_sub_buffer,
			map_among_references, 0, 1);

		///如果找到了一个匹配，就结束
		if ((*map_among_references) == 0)
		{
			return 1;
		}
	}
}



int Map_Single_Seq_end_to_end_pbat(int thread_id)
{
	///int my_methylation_size = methylation_size / THREAD_COUNT;
	int my_methylation_size = methylation_size;

	Methylation methy;


	if (output_methy == 1)
	{
		init_methylation(&methy, my_methylation_size);
	}





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
	///bitmapper_bs_iter total_best_mapping_site = 0;
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







#if defined __AVX2__

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

#else

	__m128i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

#endif






	long long total_number_of_hit = 0;




	int is_mutiple_map = 0;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	///long long direct_cut = 0;
	///long long number_of_short_seeds = 0;
	///long long number_of_short_seeds_matches = 0;
	int C_site;

	///double total_c_time = 0;
	///double start_c_time;


	double startTime = Get_T();

	int file_flag;

	int get_error;
	///long long debug_1_mismatch = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;




	int map_among_references;



	bam_output_cell cell;
	Output_buffer_sub_block current_sub_buffer;

	if (bam_output)
	{
		init_bam_output_cell(&cell);
		init_buffer_sub_block(&current_sub_buffer);
	}


	long long total_bases = 0;
	long long error_bases = 0;



	///这里要改
	bitmapper_bs_iter pre_matched_read = (bitmapper_bs_iter)-1;
	int output_mask;
	long long ambiguous_index;


	///need modification
	unsigned int second_best_diff;
	int one_mismatch_site;


	//正向模式
	i = 0;
	///while (getReads_single(&current_read))
	while (1)
	{


		///这里要改
		/********************************************输出不匹配结果用的********************************************/
		///这里要改
		///这说明上一个read没有匹配
		///初始pre_matched_read = (bitmapper_bs_iter)-1, matched_read = 0
		///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
		if (pre_matched_read == matched_read)
		{
			/**
			if(second_best_diff != 0 && map_among_references == 0)
			{
				fprintf(stderr, "2: second_best_diff: %d, candidate_length: %d, name: %s, current_read.length: %d\n", 
				second_best_diff, candidate_length, current_read.name, current_read.length);
			}
			**/

			if (unmapped_out == 1)
			{
				///这个seq和rseq要逆过来
				output_sam_unmapped(current_read.name, current_read.rseq, current_read.seq, current_read.qual,
					current_read.length, &cell, &current_sub_buffer);
			}
		}
		pre_matched_read = matched_read;
		/********************************************输出不匹配结果用的********************************************/



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

		map_among_references = 0;

		second_best_diff = 0;



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
				top, tmp_ref, &match_length, current_read.length - C_site - 1, &get_error, &methy,
				&cell,
				&current_sub_buffer, &map_among_references))
			{
				if (map_among_references == 0)
				{
					unique_matched_read++;
					matched_read++;

					total_bases = total_bases + current_read.length;
				}

				i++;
				continue;


			}


			one_mismatch_site = match_length;

			











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
					

					///这里要改
					if (ambiguous_out == 1)
					{

						output_ambiguous_exact_map_pbat(
							current_read.length,
							current_read.seq,
							current_read.rseq,
							current_read.qual,
							cigar,
							current_read.name,
							locates,
							top,
							number_of_hits,
							max_seed_matches,
							&match_length,
							&cell,
							&current_sub_buffer,
							&map_among_references);

						if (map_among_references != 0)
						{
							matched_read--;
						}

					}



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


			if (output_methy == 1)
			{
				output_methy_end_to_end_pbat
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
					current_read.length,
					&methy, &map_among_references);
			}
			else
			{

				///need modification
				int score = 0;
				
				if(current_read.seq[one_mismatch_site] == 'N')
				{
					score -= N_Penalty;
				}
				else
				{
					///seq saves the reverse complement strand of read
					///while qual saves the forward quality
					///and the calculation is based on seq
					score -= MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, Q_base, 
					current_read.qual[current_read.length - one_mismatch_site - 1]);
				}

				int mapq = MAP_Calculation((unsigned int)-1, error_threshold1, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);



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
					current_read.length,
					&cell,
					&current_sub_buffer, &map_among_references,0,mapq);

			}

			if (map_among_references == 0)
			{
				unique_matched_read++;
				matched_read++;

				total_bases = total_bases + current_read.length;
				error_bases = error_bases + 1;
			}

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
					#if defined __AVX2__
					map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#else
					map_candidate_votes_mutiple_cut_end_to_end_4_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#endif


				}
				else
				{
					#if defined __AVX2__
					map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#else
					map_candidate_votes_mutiple_cut_end_to_end_2_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#endif
				}


			}
			else
			{

				if (error_threshold1 <= 15)
				{
					#if defined __AVX2__
					map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#else
					map_candidate_votes_mutiple_end_to_end_4_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#endif

				}
				else
				{
					#if defined __AVX2__
					map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#else
					map_candidate_votes_mutiple_end_to_end_2_sse(candidates_votes, &candidates_votes_length, current_read.length,
						error_threshold1, tmp_ref, current_read.seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
					#endif
				}


			}





			min_candidates_votes_length = 0;

			if (min_err_index >= 0)
			{

				/**
				if(second_best_diff <= 0)
				{
					fprintf(stderr, "0: second_best_diff: %d\n", second_best_diff);
				}
				**/
				
				
				min_candidates_votes_length = 0;

				calculate_best_map_cigar_end_to_end_pbat
					(candidates_votes, &min_candidates_votes_length,
					current_read.length, error_threshold1, tmp_ref,
					current_read.seq, current_read.rseq, current_read.qual, cigar, path, matrix, matrix_bit, current_read.name,
					best_cigar, &methy, &cell,
					&current_sub_buffer, &map_among_references, 0, second_best_diff);

				if (map_among_references == 0)
				{
					unique_matched_read++;
					matched_read++;


					total_bases = total_bases + current_read.length;
					error_bases = error_bases + candidates_votes[0].err;
				}
				


			}
			else if (min_err_index == -1)
			{
				unmatched_read++;
			}
			else
			{

				/**
				if(second_best_diff != 0)
				{
					fprintf(stderr, "1: second_best_diff: %d\n", second_best_diff);
				}
				**/
				
				
				


				matched_read++;


				///这里要改
				if (ambiguous_out == 1)
				{
					ambiguous_index = -2 - min_err_index;
					candidates_votes[0].err = candidates_votes[ambiguous_index].err;
					candidates_votes[0].end_site = candidates_votes[ambiguous_index].end_site;
					candidates_votes[0].site = candidates_votes[ambiguous_index].site;

					calculate_best_map_cigar_end_to_end_pbat
						(candidates_votes, &min_candidates_votes_length,
						current_read.length, error_threshold1, tmp_ref,
						current_read.seq, current_read.rseq, current_read.qual, cigar, path, matrix, matrix_bit, current_read.name,
						best_cigar, &methy,
						&cell,
						&current_sub_buffer,
						&map_among_references, 0, second_best_diff);

					if (map_among_references != 0)
					{
						matched_read--;
					}
				}
			}


		}
		else
		{
			empty_read++;
		}


		i++;
	}

	if (output_methy == 1 && methy.current_size != 0)
	{
		assign_cuts(&methy);
		output_methylation(&methy);
		clear_methylation(&methy);

	}


	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;
	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;


	////fprintf(stderr, "debug_1_mismatch: %lld\n", debug_1_mismatch);


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

	if (bam_output == 0)
	{

		if (output_methy == 0)
		{
			init_output_buffer(THREAD_COUNT);
		}
		else
		{
			///init_output_methy_buffer(THREAD_COUNT);
			init_output_methy_buffer_pair(THREAD_COUNT);
		}
	}
	else
	{
		init_multiple_buffer(THREAD_COUNT);
	}

	

	pthread_t inputReadsHandle;
	init_Pair_Seq_input_buffer(THREAD_COUNT);

	///if (read_format == FASTQ)
	{
		pthread_create(&inputReadsHandle, NULL, input_pe_reads_muti_threads, NULL);
	}



	pthread_t outputResultSinkHandle;

	if (bam_output == 0)
	{

		if (output_methy == 0)
		{
			pthread_create(&outputResultSinkHandle, NULL, pop_buffer, NULL);
		}
		else
		{
			pthread_create(&outputResultSinkHandle, NULL, pop_methy_buffer_pair, NULL);
		}
	}
	else
	{
		pthread_create(&outputResultSinkHandle, NULL, pop_buffer_bam, NULL);
	}





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

	if (read_format == FASTQ)
	{
		pthread_join(inputReadsHandle, NULL);
	}

	for (i = 0; i<THREAD_COUNT; i++)
		pthread_join(_r_threads[i], NULL);

	pthread_join(outputResultSinkHandle, NULL);

	

	free(_r_threads);

	
}


int Map_Single_Seq_muti_thread(int thread_id)
{

	if (bam_output == 0)
	{
		if (output_methy == 0)
		{
			init_output_buffer(THREAD_COUNT);
		}
		else
		{
			init_output_methy_buffer(THREAD_COUNT);
		}
	}
	else
	{
		init_multiple_buffer(THREAD_COUNT);
	}

	
	
	pthread_t inputReadsHandle;

	

	init_Single_Seq_input_buffer(THREAD_COUNT);

	///if (read_format == FASTQ)
	{
		
		pthread_create(&inputReadsHandle, NULL, input_single_reads_muti_threads, NULL);
	}

	



	pthread_t outputResultSinkHandle;

	if (bam_output == 0)
	{

		if (output_methy == 0)
		{
			pthread_create(&outputResultSinkHandle, NULL, pop_buffer, NULL);
		}
		else
		{
			pthread_create(&outputResultSinkHandle, NULL, pop_methy_buffer, NULL);
		}
	}
	else
	{
		pthread_create(&outputResultSinkHandle, NULL, pop_buffer_bam, NULL);
	}

	

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

	if (read_format == FASTQ)
	{
		pthread_join(inputReadsHandle, NULL);
	}


	for (i = 0; i<THREAD_COUNT; i++)
		pthread_join(_r_threads[i], NULL);


	pthread_join(outputResultSinkHandle, NULL);
	



	free(_r_threads);


}


int Map_Single_Seq_pbat_muti_thread(int thread_id)
{

	if (bam_output == 0)
	{
		if (output_methy == 0)
		{
			init_output_buffer(THREAD_COUNT);
		}
		else
		{
			init_output_methy_buffer(THREAD_COUNT);
		}
	}
	else
	{
		init_multiple_buffer(THREAD_COUNT);
	}


	pthread_t inputReadsHandle;

	init_Single_Seq_input_buffer(THREAD_COUNT);

	///if (read_format == FASTQ)
	{
		pthread_create(&inputReadsHandle, NULL, input_single_reads_muti_threads_pbat, NULL);
	}


	pthread_t outputResultSinkHandle;

	if (bam_output == 0)
	{
		if (output_methy == 0)
		{
			pthread_create(&outputResultSinkHandle, NULL, pop_buffer, NULL);
		}
		else
		{
			pthread_create(&outputResultSinkHandle, NULL, pop_methy_buffer, NULL);
		}
	}
	else
	{
		pthread_create(&outputResultSinkHandle, NULL, pop_buffer_bam, NULL);
	}




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

	if (read_format == FASTQ)
	{
		pthread_join(inputReadsHandle, NULL);
	}


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








inline int output_ambiguous_exact_map_output_buffer(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top,
	bitmapper_bs_iter number_of_hits,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter* matched_length,
	Output_buffer_sub_block* sub_block,
	bam_phrase* bam_groups,
	int* map_among_references, int thread_id)
{



	bitmapper_bs_iter tmp_SA_length = 0;

	///限制遍历上限
	if (number_of_hits > max_seed_matches)
	{
		number_of_hits = max_seed_matches;
	}

	///这个cigar一定是精确匹配的
	sprintf(cigar, "%dM", read_length);


	bitmapper_bs_iter i = 0;

	for (i = 0; i < number_of_hits; i++)
	{

		tmp_SA_length = 0;

		locate_one_position_direct(locates, top + i, &tmp_SA_length, total_SA_length, *matched_length, 0);

		
		///in this case, mapq should be 1
		output_sam_end_to_end_output_buffer(name, read, r_read, qulity,
			locates[0],
			read_length - 1,
			0,
			0,
			cigar, read_length,
			sub_block,
			bam_groups,
			map_among_references, 0, 1);
		

		///如果找到了一个匹配，就结束
		if ((*map_among_references) == 0)
		{
			return 1;
		}
	}
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
	///bitmapper_bs_iter total_best_mapping_site = 0;

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







#if defined __AVX2__

	__m256i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm256_setzero_si256();
	}

#else

	__m128i Peq_SSE[256];

	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

#endif






	long long total_number_of_hit = 0;


	int is_mutiple_map = 0;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	///long long direct_cut = 0;
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
	Methylation methy;

	
	if (output_methy == 1)
	{
		int my_methylation_size = methylation_size;
		init_methylation(&methy, my_methylation_size);
	}
	else
	{
		init_buffer_sub_block(&current_sub_buffer);
	}



	///int inner_i;

	double startTime = Get_T();

	int file_flag = 1;
	int obtain_reads_num;


	int get_error;

	///long long debug_1_mismatch = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;



	bam_phrase bam_buffer;

	if (bam_output == 1)
	{
		init_bam_buffer(&bam_buffer);
	}

	int map_among_references;

	long long total_bases = 0;
	long long error_bases = 0;





	///这里要改
	bitmapper_bs_iter pre_matched_read = (bitmapper_bs_iter)-1;
	int output_mask;
	long long ambiguous_index;


	unsigned int second_best_diff;
	int one_mismatch_site;






	file_flag = 1;
	while (file_flag != 0)
	{


		file_flag = get_single_reads_mul_thread(&curr_sub_block);
		obtain_reads_num = curr_sub_block.sub_block_read_number;
		read_batch = curr_sub_block.read;


		
		
		enq_i = enq_i + obtain_reads_num;

		i = 0;

		///这里要改
		pre_matched_read = (bitmapper_bs_iter)-1;
		///这里要改
		while (1)
		{



			/********************************************输出不匹配结果用的********************************************/
			///这里要改
			///这说明上一个read没有匹配
			///初始pre_ambious_matched_read = (bitmapper_bs_iter)-1, ambious_matched_read = 0
			///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
			if (pre_matched_read == matched_read)
			{
				if (unmapped_out == 1)
				{
					output_sam_end_to_end_output_buffer_unmapped(
						read_batch[i - 1].name, read_batch[i - 1].seq, read_batch[i - 1].rseq, read_batch[i - 1].qual,
						read_batch[i - 1].length, &current_sub_buffer, &bam_buffer);
				}
			}
			pre_matched_read = matched_read;
			/********************************************输出不匹配结果用的********************************************/


			////这里要改
			if (i >= obtain_reads_num)
			{
				break;
			}





			
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

			map_among_references = 0;

			second_best_diff = 0;


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








		
				if (number_of_hits == 1 && try_process_unique_mismatch_end_to_end_output_buffer(
					read_batch[i].length, read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual, cigar, read_batch[i].name,
					locates,
					top, tmp_ref, &match_length, read_batch[i].length - C_site - 1, &current_sub_buffer, &get_error,
					&methy, &bam_buffer, &map_among_references))
				{
					///direct_match++;

					if (map_among_references == 0)
					{
						unique_matched_read++;
						matched_read++;

						total_bases = total_bases + read_batch[i].length;
					}

					
					i++;
					continue;


				}


				one_mismatch_site = match_length;


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
						


						
						///这里要改
						if (ambiguous_out == 1)
						{
							output_ambiguous_exact_map_output_buffer(
								read_batch[i].length,
								read_batch[i].seq,
								read_batch[i].rseq,
								read_batch[i].qual,
								cigar,
								read_batch[i].name,
								locates,
								top,
								number_of_hits,
								max_seed_matches,
								&match_length,
								&current_sub_buffer, 
								&bam_buffer,
								&map_among_references, thread_id);
						

							if (map_among_references != 0)
							{
								matched_read--;
							}

						}
						
						///注意i一定要放在output_ambiguous_exact_map_output_buffer之后啊
					    ///血的教训，不放会报错啊啊啊
						i++;
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

				if (output_methy == 1)
				{
					output_methy_end_to_end_output_buffer(
						read_batch[i].name,
						read_batch[i].seq,
						read_batch[i].rseq,
						read_batch[i].qual,
						candidates[0],
						read_batch[i].length - 1,
						0,
						1,
						cigar,
						read_batch[i].length, &methy, &map_among_references);
				}
				else
				{

					///need modification
					int score = 0;

					if(read_batch[i].seq[one_mismatch_site] == 'N')
					{
						score -= N_Penalty;
					}
					else
					{
						score -= MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, Q_base, read_batch[i].qual[one_mismatch_site]);
					}

					int mapq = MAP_Calculation((unsigned int)-1, error_threshold1, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
					/**
					debug_1_mismatch_score(candidates[0], tmp_ref, current_read.seq, current_read.length, 
						current_read.qual, Q_base, MistMatchPenaltyMax, MistMatchPenaltyMin, N_Penalty, 0, one_mismatch_site, score);
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
						read_batch[i].length, &current_sub_buffer, &bam_buffer, &map_among_references, 0, mapq);
				}

				

				if (map_among_references == 0)
				{
					unique_matched_read++;
					matched_read++;

					total_bases = total_bases + read_batch[i].length;
					error_bases = error_bases + 1;
				}

				//debug_1_mismatch++;
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
						#if defined __AVX2__
						map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_cut_end_to_end_4_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif
					}
					else
					{
						#if defined __AVX2__
						map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_cut_end_to_end_2_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif
					}


				}
				else
				{
					total_is_mutiple_map++;


					if (error_threshold1 <= 15)
					{
						#if defined __AVX2__
						map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_end_to_end_4_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif

					}
					else
					{
						#if defined __AVX2__
						map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_end_to_end_2_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif
					}


				}


				min_candidates_votes_length = 0;

				if (min_err_index >= 0)
				{
					
					min_candidates_votes_length = 0;


					calculate_best_map_cigar_end_to_end_output_buffer
						(candidates_votes, &min_candidates_votes_length,
						read_batch[i].length, error_threshold1, tmp_ref,
						read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual,
						cigar, path, matrix, matrix_bit, read_batch[i].name,
						best_cigar, &current_sub_buffer, &methy, &bam_buffer, &map_among_references,0,
						second_best_diff);

					if (map_among_references == 0)
					{
						matched_read++;
						unique_matched_read++;

						total_bases = total_bases + read_batch[i].length;
						error_bases = error_bases + candidates_votes[0].err;
					}
					


				}
				else if (min_err_index == -1)
				{
					unmatched_read++;
				}
				else
				{
					matched_read++;


					///这里要改
					if (ambiguous_out == 1)
					{
						ambiguous_index = -2 - min_err_index;
						candidates_votes[0].err = candidates_votes[ambiguous_index].err;
						candidates_votes[0].end_site = candidates_votes[ambiguous_index].end_site;
						candidates_votes[0].site = candidates_votes[ambiguous_index].site;

						calculate_best_map_cigar_end_to_end_output_buffer
							(candidates_votes, &min_candidates_votes_length,
							read_batch[i].length, error_threshold1, tmp_ref,
							read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual,
							cigar, path, matrix, matrix_bit, read_batch[i].name,
							best_cigar, &current_sub_buffer, &methy, &bam_buffer, &map_among_references, 0,
							second_best_diff);

						if (map_among_references != 0)
						{
							matched_read--;
						}
					}



				}







			}
			else
			{
				empty_read++;
			}

			i++;

		}


		if (output_methy == 0 && bam_output == 0)
		{

			current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

			push_results_to_buffer(&current_sub_buffer);

			current_sub_buffer.length = 0;
		}

	}


	if (output_methy == 1 && methy.current_size != 0)
	{
		assign_cuts(&methy);

		push_methy_to_buffer(&methy);

		///这个不用清理，因为再push 的时候清理过了
		///clear_methylation(&methy);

	}


	if (bam_output == 1)
	{
		flush_bam_buffer(&bam_buffer);
	}


	if (bam_output == 1)
	{
		finish_bam_output_buffer();

	}
	else
	{
		finish_output_buffer();
	}

	



	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;

	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;



	fprintf(stderr, "thread %d completed!\n", thread_id);


}




inline int output_ambiguous_exact_map_output_buffer_pbat(
	bitmapper_bs_iter read_length, char* read, char* r_read, char* qulity, char* cigar, char* name,
	bitmapper_bs_iter* locates,
	bitmapper_bs_iter top,
	bitmapper_bs_iter number_of_hits,
	bitmapper_bs_iter max_seed_matches,
	bitmapper_bs_iter* matched_length,
	Output_buffer_sub_block* sub_block,
	bam_phrase* bam_groups,
	int* map_among_references, int thread_id)
{



	bitmapper_bs_iter tmp_SA_length = 0;

	///限制遍历上限
	if (number_of_hits > max_seed_matches)
	{
		number_of_hits = max_seed_matches;
	}

	///这个cigar一定是精确匹配的
	sprintf(cigar, "%dM", read_length);


	bitmapper_bs_iter i = 0;

	for (i = 0; i < number_of_hits; i++)
	{

		tmp_SA_length = 0;

		locate_one_position_direct(locates, top + i, &tmp_SA_length, total_SA_length, *matched_length, 0);



		output_sam_end_to_end_pbat_output_buffer(name, read, r_read, qulity,
			locates[0],
			read_length - 1,
			0,
			0,
			cigar, read_length,
			sub_block,
			bam_groups,
			map_among_references, 0, 1);


		///如果找到了一个匹配，就结束
		if ((*map_among_references) == 0)
		{
			return 1;
		}
	}
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
	///bitmapper_bs_iter total_best_mapping_site = 0;

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




	#if defined __AVX2__

		__m256i Peq_SSE[256];

		for (i = 0; i < 256; i++)
		{
			Peq_SSE[i] = _mm256_setzero_si256();
		}

	#else

		__m128i Peq_SSE[256];

		for (i = 0; i < 256; i++)
		{
			Peq_SSE[i] = _mm_setzero_si128();
		}

	#endif




	long long total_number_of_hit = 0;


	int is_mutiple_map = 0;

	long long total_is_mutiple_map = 0;
	long long total_cutted = 0;
	///long long direct_cut = 0;
	///long long direct_match = 0;
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
	Methylation methy;


	if (output_methy == 1)
	{
		int my_methylation_size = methylation_size;
		init_methylation(&methy, my_methylation_size);
	}
	else
	{
		init_buffer_sub_block(&current_sub_buffer);
	}






	///int inner_i;

	double startTime = Get_T();

	int file_flag = 1;
	int obtain_reads_num;


	int get_error;
	///long long debug_1_mismatch = 0;
	long long first_seed_match_length = 0;
	long long second_seed_length = 0;
	int extra_seed_flag = 1;





	bam_phrase bam_buffer;

	if (bam_output == 1)
	{
		init_bam_buffer(&bam_buffer);
	}


	int map_among_references;



	long long total_bases = 0;
	long long error_bases = 0;



	///这里要改
	bitmapper_bs_iter pre_matched_read = (bitmapper_bs_iter)-1;
	int output_mask;
	long long ambiguous_index;

	///need modification
	unsigned int second_best_diff;
	int one_mismatch_site;



	file_flag = 1;
	while (file_flag != 0)
	{


		file_flag = get_single_reads_mul_thread_pbat(&curr_sub_block);
		obtain_reads_num = curr_sub_block.sub_block_read_number;
		read_batch = curr_sub_block.read;


		enq_i = enq_i + obtain_reads_num;

		i = 0;


		///这里要改
		pre_matched_read = (bitmapper_bs_iter)-1;
		///这里要改
		while (1)
		{




			/********************************************输出不匹配结果用的********************************************/
			///这里要改
			///这说明上一个read没有匹配
			///初始pre_ambious_matched_read = (bitmapper_bs_iter)-1, ambious_matched_read = 0
			///这样就避免了初始的问题, 因为初始情况下不该输出任何东西
			if (pre_matched_read == matched_read)
			{
				if (unmapped_out == 1)
				{
					///seq和rseq要逆过来
					output_sam_end_to_end_output_buffer_unmapped(
						read_batch[i - 1].name, read_batch[i - 1].rseq, read_batch[i - 1].seq, read_batch[i - 1].qual,
						read_batch[i - 1].length, &current_sub_buffer, &bam_buffer);
				}
			}
			pre_matched_read = matched_read;
			/********************************************输出不匹配结果用的********************************************/


			////这里要改
			if (i >= obtain_reads_num)
			{
				break;
			}




















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

			map_among_references = 0;

			second_best_diff = 0;



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
					top, tmp_ref, &match_length, read_batch[i].length - C_site - 1, &current_sub_buffer, &get_error, &methy,
					&bam_buffer, &map_among_references))
				{

					if (map_among_references == 0)
					{
						unique_matched_read++;
						matched_read++;

						total_bases = total_bases + read_batch[i].length;
					}

				
					i++;
					continue;


				}

				one_mismatch_site = match_length;





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




						///这里要改
						if (ambiguous_out == 1)
						{
							output_ambiguous_exact_map_output_buffer_pbat(
								read_batch[i].length,
								read_batch[i].seq,
								read_batch[i].rseq,
								read_batch[i].qual,
								cigar,
								read_batch[i].name,
								locates,
								top,
								number_of_hits,
								max_seed_matches,
								&match_length,
								&current_sub_buffer,
								&bam_buffer,
								&map_among_references, thread_id);


							if (map_among_references != 0)
							{
								matched_read--;
							}

						}






						i++;
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


				if (output_methy == 1)
				{
					output_methy_end_to_end_output_buffer
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
						&methy, &map_among_references);
				}
				else
				{


					///need modification
					int score = 0;
					
					if(read_batch[i].seq[one_mismatch_site] == 'N')
					{
						score -= N_Penalty;
					}
					else
					{
						///seq saves the reverse complement strand of read
						///while qual saves the forward quality
						///and the calculation is based on seq
						score -= MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, Q_base, 
						read_batch[i].qual[read_batch[i].length - one_mismatch_site - 1]);
					}

					int mapq = MAP_Calculation((unsigned int)-1, error_threshold1, score, GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);


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
						&current_sub_buffer,
						&bam_buffer, &map_among_references,0,mapq);
				}

				if (map_among_references == 0)
				{

					unique_matched_read++;
					matched_read++;

					total_bases = total_bases + read_batch[i].length;
					error_bases = error_bases + 1;
				}

				///debug_1_mismatch++;
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
						#if defined __AVX2__
						map_candidate_votes_mutiple_cut_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_cut_end_to_end_4_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif


					}
					else
					{
						#if defined __AVX2__
						map_candidate_votes_mutiple_cut_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_cut_end_to_end_2_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif
					}


				}
				else
				{
					total_is_mutiple_map++;


					if (error_threshold1 <= 15)
					{
						#if defined __AVX2__
						map_candidate_votes_mutiple_end_to_end_8(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_end_to_end_4_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif

					}
					else
					{
						#if defined __AVX2__
						map_candidate_votes_mutiple_end_to_end_4(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#else
						map_candidate_votes_mutiple_end_to_end_2_sse(candidates_votes, &candidates_votes_length, read_batch[i].length,
							error_threshold1, tmp_ref, read_batch[i].seq, Peq_SSE, &min_err, &min_err_index, &second_best_diff);
						#endif
					}


				}


				min_candidates_votes_length = 0;

				if (min_err_index >= 0)
				{

					///total_best_mapping_site++;
					
					min_candidates_votes_length = 0;

					calculate_best_map_cigar_end_to_end_pbat_output_buffer
						(candidates_votes, &min_candidates_votes_length,
						read_batch[i].length, error_threshold1, tmp_ref,
						read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual,
						cigar, path, matrix, matrix_bit, read_batch[i].name,
						best_cigar, &current_sub_buffer, &methy,
						&bam_buffer, &map_among_references,0, second_best_diff);

					if (map_among_references == 0)
					{
						matched_read++;
						unique_matched_read++;

						total_bases = total_bases + read_batch[i].length;
						error_bases = error_bases + candidates_votes[0].err;
					}

					


				}
				else if (min_err_index == -1)
				{
					unmatched_read++;
				}
				else
				{
					matched_read++;



					///这里要改
					if (ambiguous_out == 1)
					{
						ambiguous_index = -2 - min_err_index;
						candidates_votes[0].err = candidates_votes[ambiguous_index].err;
						candidates_votes[0].end_site = candidates_votes[ambiguous_index].end_site;
						candidates_votes[0].site = candidates_votes[ambiguous_index].site;

						calculate_best_map_cigar_end_to_end_pbat_output_buffer
							(candidates_votes, &min_candidates_votes_length,
							read_batch[i].length, error_threshold1, tmp_ref,
							read_batch[i].seq, read_batch[i].rseq, read_batch[i].qual,
							cigar, path, matrix, matrix_bit, read_batch[i].name,
							best_cigar, &current_sub_buffer, &methy, &bam_buffer, &map_among_references, 0, second_best_diff);

						if (map_among_references != 0)
						{
							matched_read--;
						}
					}

				}







			}
			else
			{
				empty_read++;
			}

			i++;

		}


		if (output_methy == 0 && bam_output == 0)
		{
			current_sub_buffer.buffer[current_sub_buffer.length] = '\0';

			push_results_to_buffer(&current_sub_buffer);

			current_sub_buffer.length = 0;
		}

	}

	if (output_methy == 1 && methy.current_size != 0)
	{
		assign_cuts(&methy);

		push_methy_to_buffer(&methy);

		///这个不用清理，因为再push 的时候清理过了
		///clear_methylation(&methy);

	}






	if (bam_output == 1)
	{
		flush_bam_buffer(&bam_buffer);
	}


	if (bam_output == 1)
	{
		finish_bam_output_buffer();

	}
	else
	{
		finish_output_buffer();
	}











	double totalMappingTime;
	totalMappingTime = Get_T() - startTime;

	completedSeqCnt[thread_id] = enq_i;
	unique_mapped_read[thread_id] = unique_matched_read;
	ambiguous_mapped_read[thread_id] = matched_read - unique_matched_read;

	mapped_bases[thread_id] = total_bases;
	error_mapped_bases[thread_id] = error_bases;

	fprintf(stderr, "thread %d completed!\n", thread_id);

}



void init_methy_input_buffer(long long each_sub_buffer_size, int number_of_threads,
	bitmapper_bs_iter total_start, bitmapper_bs_iter total_end, int is_update, int is_paired)
{
	int i = 0;

	bitmapper_bs_iter total_length = total_end - total_start + 1;

	///int my_methylation_size = methylation_size / 2;
	int my_methylation_size = methylation_size;

	methy_input_buffer.number_of_intervals = number_of_threads;

	methy_input_buffer.each_interval_length = total_length / number_of_threads + (total_length % number_of_threads != 0);

	methy_input_buffer.end = 0;

	methy_input_buffer.all_end = 0;

	methy_input_buffer.all_completed = 0;

	methy_input_buffer.total_start = total_start;
	methy_input_buffer.total_end = total_end;
	

	////如果不是update，那就要初始化，就要分配空间
	if (is_update == 0)
	{
		methy_input_buffer.intervals
			= (Pair_Methylation*)malloc(sizeof(Pair_Methylation)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.each_buffer_interval_start
			= (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.each_buffer_interval_end
			= (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*methy_input_buffer.number_of_intervals);


		methy_input_buffer.M_methy
			= (uint16_t**)malloc(sizeof(uint16_t*)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.M_total
			= (uint16_t**)malloc(sizeof(uint16_t*)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.M_correct_methy
			= (uint16_t**)malloc(sizeof(uint16_t*)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.M_correct_total
			= (uint16_t**)malloc(sizeof(uint16_t*)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.M_ref_genome
			= (char**)malloc(sizeof(char*)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.M_start_pos
			= (long long*)malloc(sizeof(long long)*methy_input_buffer.number_of_intervals);
		methy_input_buffer.M_end_pos
			= (long long*)malloc(sizeof(long long)*methy_input_buffer.number_of_intervals);

		methy_input_buffer.need_context[0] = CpG;
		methy_input_buffer.need_context[1] = CHG;
		methy_input_buffer.need_context[2] = CHH;

		if (CpG)
		{
			methy_input_buffer.CpG_buffer
				= (Output_buffer_sub_block*)malloc(sizeof(Output_buffer_sub_block)*methy_input_buffer.number_of_intervals);
		}

		if (CHG)
		{
			methy_input_buffer.CHG_buffer
				= (Output_buffer_sub_block*)malloc(sizeof(Output_buffer_sub_block)*methy_input_buffer.number_of_intervals);
		}

		if (CHH)
		{
			methy_input_buffer.CHH_buffer
				= (Output_buffer_sub_block*)malloc(sizeof(Output_buffer_sub_block)*methy_input_buffer.number_of_intervals);
		}
		


		methy_input_buffer.Mutex
			= (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t) * methy_input_buffer.number_of_intervals);

		methy_input_buffer.stallCond 
			= (pthread_cond_t*)malloc(sizeof(pthread_cond_t) * methy_input_buffer.number_of_intervals);

		methy_input_buffer.flushCond
			= (pthread_cond_t*)malloc(sizeof(pthread_cond_t) * methy_input_buffer.number_of_intervals);


		methy_input_buffer.process_Mutex
			= (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)* methy_input_buffer.number_of_intervals);

		methy_input_buffer.process_stallCond
			= (pthread_cond_t*)malloc(sizeof(pthread_cond_t)* methy_input_buffer.number_of_intervals);

		
		methy_input_buffer.process_Mutex_methy
			= (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)* methy_input_buffer.number_of_intervals);

		methy_input_buffer.process_stallCond_methy
			= (pthread_cond_t*)malloc(sizeof(pthread_cond_t)* methy_input_buffer.number_of_intervals);




		for (i = 0; i < methy_input_buffer.number_of_intervals; i++)
		{

			if (is_paired)
			{
				init_pair_methylation(&methy_input_buffer.intervals[i], my_methylation_size);
			}
			else
			{
				init_single_methylation(&methy_input_buffer.intervals[i], my_methylation_size);
			}

			
			if (CpG)
			{
				init_buffer_sub_block(&methy_input_buffer.CpG_buffer[i]);
			}
			if (CHG)
			{
				init_buffer_sub_block(&methy_input_buffer.CHG_buffer[i]);
			}
			if (CHH)
			{
				init_buffer_sub_block(&methy_input_buffer.CHH_buffer[i]);
			}
		}
	}
	else  ///如果不是初次，就不必要分配空间
	{
		for (i = 0; i < methy_input_buffer.number_of_intervals; i++)
		{

			methy_input_buffer.intervals[i].current_size = 0;

		}

	}




	bitmapper_bs_iter tmp_start = total_start;
	for (i = 0; i < methy_input_buffer.number_of_intervals; i++)
	{
		methy_input_buffer.each_buffer_interval_start[i] = tmp_start;
		if (methy_input_buffer.each_buffer_interval_start[i] > total_end)
		{
			methy_input_buffer.each_buffer_interval_start[i] = total_end;
		}


		methy_input_buffer.each_buffer_interval_end[i]
			= methy_input_buffer.each_buffer_interval_start[i] + methy_input_buffer.each_interval_length - 1;
		if (methy_input_buffer.each_buffer_interval_end[i] > total_end)
		{
			methy_input_buffer.each_buffer_interval_end[i] = total_end;
		}

		tmp_start = tmp_start + methy_input_buffer.each_interval_length;
	}
}






void* input_methy_muti_threads(void*)
{

	FILE* read_file;

	uint16_t flag1;
	bitmapper_bs_iter pos1;
	int seq_length1;
	bitmapper_bs_iter read1[SEQ_MAX_LENGTH];

	uint16_t flag2;
	bitmapper_bs_iter pos2;
	int seq_length2; 
	bitmapper_bs_iter read2[SEQ_MAX_LENGTH];

	bitmapper_bs_iter key_pos;
	bitmapper_bs_iter sub_interval_ID;

	

	buffer_3_bit tmp_read1;
	buffer_3_bit tmp_read2;

	Init_buffer_3_bit(tmp_read1);
	Init_buffer_3_bit(tmp_read2);

	buffer_3_bit tmp_swap;



	Pair_Methylation* current_interval;

	int round_ID = 0;

	while (methy_input_buffer.all_end == 0)
	{

		read_file = input_methy_alignment_file;

		while (1)
		{
			if (!input_PE_alignment(read_file,
				&flag1, &pos1, &seq_length1, read1,
				&flag2, &pos2, &seq_length2, read2))
			{
				break;
			}

			Check_Space_3_Bit(tmp_read1.buffer, tmp_read1.size, seq_length1);
			Check_Space_3_Bit(tmp_read2.buffer, tmp_read2.size, seq_length2);
			memcpy(tmp_read1.buffer, read1, tmp_read1.size * sizeof(bitmapper_bs_iter));
			memcpy(tmp_read2.buffer, read2, tmp_read2.size * sizeof(bitmapper_bs_iter));


			key_pos = pos1<=pos2? pos1:pos2;
			key_pos = key_pos - methy_input_buffer.total_start;
			sub_interval_ID = key_pos / methy_input_buffer.each_interval_length;

			current_interval = methy_input_buffer.intervals + sub_interval_ID;



			pthread_mutex_lock(&methy_input_buffer.Mutex[sub_interval_ID]);
			///如果满了，则该线程本身要wait，并通知消费者线程消费数据
			while (current_interval->current_size
					== current_interval->total_size)
			{

				

				///按道理这个信号量似乎不用发
				///因为队列不可能一边满一边空
				pthread_cond_signal(&methy_input_buffer.stallCond[sub_interval_ID]);


				pthread_cond_wait(&methy_input_buffer.flushCond[sub_interval_ID], 
					&methy_input_buffer.Mutex[sub_interval_ID]);
			}

			/**************************read1***********************************/
			tmp_swap.buffer = current_interval->R[0].reads_3_bit[current_interval->current_size];
			tmp_swap.size = current_interval->R[0].r_size_3_bit[current_interval->current_size];

			current_interval->R[0].reads_3_bit[current_interval->current_size] = tmp_read1.buffer;
			current_interval->R[0].r_size_3_bit[current_interval->current_size] = tmp_read1.size;

			tmp_read1 = tmp_swap;


			current_interval->R[0].r_length[current_interval->current_size] = flag1;
			current_interval->R[0].sites[current_interval->current_size] = pos1;
			current_interval->R[0].r_real_length[current_interval->current_size] = seq_length1;
			/**************************read1***********************************/
			 
			/**************************read2***********************************/

			tmp_swap.buffer = current_interval->R[1].reads_3_bit[current_interval->current_size];
			tmp_swap.size = current_interval->R[1].r_size_3_bit[current_interval->current_size];

			current_interval->R[1].reads_3_bit[current_interval->current_size] = tmp_read2.buffer;
			current_interval->R[1].r_size_3_bit[current_interval->current_size] = tmp_read2.size;

			tmp_read2 = tmp_swap;


			current_interval->R[1].r_length[current_interval->current_size] = flag2;
			current_interval->R[1].sites[current_interval->current_size] = pos2;
			current_interval->R[1].r_real_length[current_interval->current_size] = seq_length2;
			/**************************read2***********************************/
			current_interval->current_size++;


			pthread_cond_signal(&methy_input_buffer.stallCond[sub_interval_ID]);
			pthread_mutex_unlock(&methy_input_buffer.Mutex[sub_interval_ID]);

		}


		methy_input_buffer.end = 1;
		int i;
		for (i = 0; i < methy_input_buffer.number_of_intervals; i++)
		{
			pthread_cond_signal(&methy_input_buffer.stallCond[i]);
		}


		pthread_mutex_lock(&methy_input_buffer.all_completed_Mutex);
		if (methy_input_buffer.all_completed == THREAD_COUNT)
		{
			pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex);
			pthread_cond_signal(&methy_input_buffer.main_thread_flushCond);
			pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex);
		}
		methy_input_buffer.all_completed++;
		pthread_mutex_unlock(&methy_input_buffer.all_completed_Mutex);


		pthread_mutex_lock(&methy_input_buffer.input_thread_Mutex);
		pthread_cond_wait(&methy_input_buffer.input_thread_flushCond,&methy_input_buffer.input_thread_Mutex);
		pthread_mutex_unlock(&methy_input_buffer.input_thread_Mutex);

		round_ID++;

	}


	
}






void* input_methy_muti_threads_single_end(void*)
{

	FILE* read_file;

	uint16_t flag;
	bitmapper_bs_iter pos;
	int seq_length;
	bitmapper_bs_iter read[SEQ_MAX_LENGTH];

	

	bitmapper_bs_iter key_pos;
	bitmapper_bs_iter sub_interval_ID;



	buffer_3_bit tmp_read;

	Init_buffer_3_bit(tmp_read);

	buffer_3_bit tmp_swap;



	Pair_Methylation* current_interval;

	int round_ID = 0;

	while (methy_input_buffer.all_end == 0)
	{

		read_file = input_methy_alignment_file;

		while (1)
		{
			

			if (!input_single_alignment(read_file,
				&flag, &pos, &seq_length, read))
			{
				break;
			}

			Check_Space_3_Bit(tmp_read.buffer, tmp_read.size, seq_length);
			memcpy(tmp_read.buffer, read, tmp_read.size * sizeof(bitmapper_bs_iter));

			/**
			if (pos == 33955283)
			{
				pos = pos - methy_input_buffer.total_start;
				sub_interval_ID = pos / methy_input_buffer.each_interval_length;
				current_interval = methy_input_buffer.intervals + sub_interval_ID;

				fprintf(stderr, "total_start: %llu\n", methy_input_buffer.total_start);
				fprintf(stderr, "each_interval_length: %llu\n", methy_input_buffer.each_interval_length);
				fprintf(stderr, "sub_interval_ID: %llu\n", sub_interval_ID);
			}
			else
			{
				pos = pos - methy_input_buffer.total_start;
				sub_interval_ID = pos / methy_input_buffer.each_interval_length;
				current_interval = methy_input_buffer.intervals + sub_interval_ID;
			}
			**/
			/**
			fprintf(stderr, "each_interval_length: %llu\n", methy_input_buffer.each_interval_length);
			fprintf(stderr, "sub_interval_ID: %llu\n", sub_interval_ID);
			**/


			///pos = pos - methy_input_buffer.total_start;
			sub_interval_ID = (pos - methy_input_buffer.total_start) / methy_input_buffer.each_interval_length;
			current_interval = methy_input_buffer.intervals + sub_interval_ID;

			///fprintf(stderr, "sub_interval_ID: %d\n", sub_interval_ID);



			pthread_mutex_lock(&methy_input_buffer.Mutex[sub_interval_ID]);
			///如果满了，则该线程本身要wait，并通知消费者线程消费数据
			while (current_interval->current_size
				== current_interval->total_size)
			{



				///按道理这个信号量似乎不用发
				///因为队列不可能一边满一边空
				pthread_cond_signal(&methy_input_buffer.stallCond[sub_interval_ID]);


				pthread_cond_wait(&methy_input_buffer.flushCond[sub_interval_ID],
					&methy_input_buffer.Mutex[sub_interval_ID]);
			}

			/**************************read1***********************************/
			tmp_swap.buffer = current_interval->R[0].reads_3_bit[current_interval->current_size];
			tmp_swap.size = current_interval->R[0].r_size_3_bit[current_interval->current_size];

			current_interval->R[0].reads_3_bit[current_interval->current_size] = tmp_read.buffer;
			current_interval->R[0].r_size_3_bit[current_interval->current_size] = tmp_read.size;

			tmp_read = tmp_swap;


			current_interval->R[0].r_length[current_interval->current_size] = flag;
			current_interval->R[0].sites[current_interval->current_size] = pos;
			current_interval->R[0].r_real_length[current_interval->current_size] = seq_length;
			/**************************read1***********************************/

			
			current_interval->current_size++;


			pthread_cond_signal(&methy_input_buffer.stallCond[sub_interval_ID]);
			pthread_mutex_unlock(&methy_input_buffer.Mutex[sub_interval_ID]);

		}


		methy_input_buffer.end = 1;
		int i;
		for (i = 0; i < methy_input_buffer.number_of_intervals; i++)
		{
			pthread_cond_signal(&methy_input_buffer.stallCond[i]);
		}


		pthread_mutex_lock(&methy_input_buffer.all_completed_Mutex);
		if (methy_input_buffer.all_completed == THREAD_COUNT)
		{
			pthread_mutex_lock(&methy_input_buffer.main_thread_Mutex);
			pthread_cond_signal(&methy_input_buffer.main_thread_flushCond);
			pthread_mutex_unlock(&methy_input_buffer.main_thread_Mutex);
		}
		methy_input_buffer.all_completed++;
		pthread_mutex_unlock(&methy_input_buffer.all_completed_Mutex);


		pthread_mutex_lock(&methy_input_buffer.input_thread_Mutex);
		pthread_cond_wait(&methy_input_buffer.input_thread_flushCond, &methy_input_buffer.input_thread_Mutex);
		pthread_mutex_unlock(&methy_input_buffer.input_thread_Mutex);

		round_ID++;

	}
}


