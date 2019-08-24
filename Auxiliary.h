/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/


#ifndef __COMMON__
#define __COMMON__

#include "bwt.h"


#define READS_QUENUE_MAX_LENGTH		10000
#define SEQ_MAX_LENGTH		1000			// Seq Max Length
#define NAME_LENGTH	1000			// Filename Max Length
#define MAX_Thread 1024


extern unsigned char	WINDOW_SIZE;		// WINDOW SIZE for indexing/searching
extern unsigned char    default_ws;
extern unsigned char	MERGE_STEP;
extern unsigned short	SEQ_LENGTH;		// Sequence(read) length

extern char				*versionN;


extern unsigned int		THREAD_COUNT;
extern int				is_index;
extern int				is_search;
extern int				is_pairedEnd;
extern int				cropSize;
extern char 			*Read_File1;
extern char				*Read_File2;
extern char				*seqUnmapped;
extern char				*Mapped_File;
extern char 			*Mapped_FilePath;
extern char                *folder_path;
extern char				*Un_Mapped_File;
extern unsigned char	thread_e;
extern double thread_e_f;
extern unsigned char	maxHits;
extern int				minDistance_pair;
extern int				maxDistance_pair;
extern char				fileName[2][NAME_LENGTH];
extern int				fileCnt;
extern int each_length;
extern double bs_score_threshold;
extern double bs_edit_distance_threshold;
extern int is_local;
extern int bs_available_seed_length;
extern int pbat;
extern char rc_table[128];
extern int sub_block_inner_size;
extern int sub_block_number;
extern int ouput_buffer_size;
extern int	over_all_seed_length;
extern int methylation_size;
extern int output_methy;
extern bitmapper_bs_iter genome_cuts;
extern bitmapper_bs_iter cut_length;
extern double methylation_buffer_times;
extern int is_methy;
extern int bam_output;
extern double maxVariantFrac;
extern int minVariantDepth;
extern int CpG;
extern int CHG;
extern int CHH;
extern int unmapped_out;
extern int ambiguous_out;
extern int mapstats;


extern int GapOpenPenalty;
extern int GapExtensionPenalty;
extern int MistMatchPenaltyMax;
extern int MistMatchPenaltyMin;
extern int N_Penalty;
extern int Q_base;

extern char *Mapstats_File;
extern char *Mapstats_File_Path;



double	Get_T(void);
void	reverseComplement (char *seq, char *rcSeq , int length);
///void    reverse_pattern(char* pattern, char* rev_pattern, int length);
void 	reverse (char *seq, char *rcSeq , int length);
void 	stripPath(char *full, char **path, char **fileName);
inline void reverse_pattern(char* pattern, char* rev_pattern, int length)
{
	int i = 0;

	for (i = 0; i < length; i++)
	{
		rev_pattern[i] = pattern[length - i - 1];
	}
	rev_pattern[i] = '\0';
}
#endif
