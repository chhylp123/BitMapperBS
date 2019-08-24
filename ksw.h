#ifndef __AC_KSW_H
#define __AC_KSW_H

#include <stdint.h>
#include "bwt.h"

#define KSW_XBYTE  0x10000
#define KSW_XSTOP  0x20000
#define KSW_XSUBO  0x40000
#define KSW_XSTART 0x80000

#define MINUS_INF -0x40000000
//#define MINUS_INF_SSE -32767
#define MINUS_INF_SSE -30000


struct _kswq_t;
typedef struct _kswq_t kswq_t;

/********************
 *** SW extension ***
 ********************/

typedef struct {
	int32_t h, e;
} eh_t;

typedef struct {
	int score; // best score
	int te, qe; // target end and query end
	int score2, te2; // second best score and ending position on the target
	int tb, qb; // target start and query start
	int NM;

	///a column of h and z
	eh_t *eh;
	int eh_size;

	int8_t *qp; // query profile
	int8_t *qp_diff; // query profile
	///qp_size = size(qp) + size(qp_diff)
	int qp_size;

	uint8_t *z; 
	int z_size;


	uint32_t *cigar;
	int cigar_size;

} kswr_t;

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * Aligning two sequences
	 *
	 * @param qlen    length of the query sequence (typically <tlen)
	 * @param query   query sequence with 0 <= query[i] < m
	 * @param tlen    length of the target sequence
	 * @param target  target sequence
	 * @param m       number of residue types
	 * @param mat     m*m scoring matrix in one-dimention array
	 * @param gapo    gap open penalty; a gap of length l cost "-(gapo+l*gape)"
	 * @param gape    gap extension penalty
	 * @param xtra    extra information (see below)
	 * @param qry     query profile (see below)
	 *
	 * @return        alignment information in a struct; unset values to -1
	 *
	 * When xtra==0, ksw_align() uses a signed two-byte integer to store a
	 * score and only finds the best score and the end positions. The 2nd best
	 * score or the start positions are not attempted. The default behavior can
	 * be tuned by setting KSW_X* flags:
	 *
	 *   KSW_XBYTE:  use an unsigned byte to store a score. If overflow occurs,
	 *               kswr_t::score will be set to 255
	 *
	 *   KSW_XSUBO:  track the 2nd best score and the ending position on the
	 *               target if the 2nd best is higher than (xtra&0xffff)
	 *
	 *   KSW_XSTOP:  stop if the maximum score is above (xtra&0xffff)
	 *
	 *   KSW_XSTART: find the start positions
	 *
	 * When *qry==NULL, ksw_align() will compute and allocate the query profile
	 * and when the function returns, *qry will point to the profile, which can
	 * be deallocated simply by free(). If one query is aligned against multiple
	 * target sequences, *qry should be set to NULL during the first call and
	 * freed after the last call. Note that qry can equal 0. In this case, the
	 * query profile will be deallocated in ksw_align().
	 */
	kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry);

	int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int h0, int *_qle, int *_tle);
	int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *_n_cigar, uint32_t **_cigar);
	int ksw_semi_global(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);
	int ksw_semi_global_sse(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);
	int ksw_semi_global_sse_debug(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);
int ksw_semi_global_sse_debug(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);
int ksw_semi_global_sse_debug_debug(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);


int ksw_semi_global_sse_new(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_, int flag);

int ksw_semi_global_path(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);

int fast_recalculate_bs_Cigar(
	char *pattern, int p_length,
	char *text, int t_length,
	unsigned short errthold,
	int end_site,
	int error,
	int* return_start_site,
	bitmapper_bs_iter* return_end_site,
	unsigned int* return_err,
	int* score,
	char* cigar,
	int is_forward_strand,
	int8_t* mat,
	int8_t* mat_diff,
	int gapo, 
	int gape,
	int MistMatchPenaltyMax,
	int MistMatchPenaltyMin,
	int N_Penalty,
	char* qulity,
	int need_r_quality,
	int quality_base
	);

int ksw_semi_global_quality(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat,
const int8_t *mat_diff, int gapo, int gape, int w, char* qulity, int quality_base, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);



int ksw_semi_global_quality_back(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat,
const int8_t *mat_diff, int gapo, int gape, int w, char* qulity, int quality_base, kswr_t* qry, int *n_cigar_, uint32_t **cigar_);


inline int MismatchPenaltyByQuality(int MistMatchPenaltyMax, int MistMatchPenaltyMin, int quality_base, int Q)
{
	double Phred;
	Phred = Q - quality_base;
	if (Phred > 40)
	{
		Phred = 40;
	}
	Phred = Phred/40;
	int return_Q = Phred*(MistMatchPenaltyMax - MistMatchPenaltyMin);
	return_Q = return_Q + MistMatchPenaltyMin;

	return return_Q;
}

inline void init_qry_total(kswr_t* qry)
{
	qry->eh_size = 0;
	qry->eh = NULL;

	qry->qp_size = 0;
	qry->qp = NULL;
	qry->qp_diff = NULL;

	qry->z_size = 0;
	qry->z = NULL;

	qry->cigar_size = 0;
	qry->cigar = NULL;
}

#ifdef __cplusplus
}
#endif

#endif
