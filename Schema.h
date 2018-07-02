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

struct result_queue{
    TAILQ_ENTRY(result_queue) result_tailq;
    tmp_result single_result;
};

TAILQ_HEAD(Result_Queue,result_queue);

//struct Result_Queue queue_header;  ///queue_header是队列的头部

extern long long			mappingCnt[MAX_Thread];
extern long long			mappedSeqCnt[MAX_Thread];
extern long long			completedSeqCnt[MAX_Thread];

void Prepare_alignment(char *genFileName,_rg_name_l *chhy_ih_refGenName,int chhy_refChromeCont);
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


void get_mapping_informations(long long* number_of_read, long long* number_of_unique_mapped_read,
	long long* number_of_ambiguous_mapped_read, long long* number_of_unmapped_read);



///下面是我到函数



unsigned int get_total_reads_number();




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





#endif
