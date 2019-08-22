/* The MIT License

   Copyright (c) 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <stdio.h>
#include <string.h>
#include "ksw.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif


unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


const kswr_t g_defr = { 0, -1, -1, -1, -1, -1, -1 };








/********************
 * Global alignment *
 ********************/


static inline uint32_t *push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = (uint32_t *)realloc(cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

typedef union __m128i_16 {
    __m128i m;
    int16_t v[8];
} __m128i_16_t;

typedef struct {
	__m128i* alloc;
	__m128i* mem;
} aligned_sse;

inline void malloc_16bit_aligned(aligned_sse* ptr, long long size)
{
	ptr->alloc = (__m128i*)malloc(sizeof(__m128i)*size + 2);
	ptr->mem = (__m128i*)(((uintptr_t)ptr->alloc + 1) & (~ (uintptr_t)1));
}

inline void free_16bit_aligned(aligned_sse* ptr)
{
	free(ptr->alloc);
}

///qlen = tlen + 2*w
int ksw_semi_global_sse(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	int beg, end;
	int i, j, k = MINUS_INF_SSE, gapoe = gapo + gape, score;
	///可以把band_length扩展到能被4整除的数
	int band_length = 2*w+1;
	///单个sse element大小为16 bit, 所以sse vector中共有8个元素
    const int32_t segWidth = 8; /* number of values in vector unit */
    ///sse vector的个数
    const int32_t segLen = (band_length + segWidth - 1) / segWidth;
	aligned_sse vProfile, vProfile_tmp, pvHStore, pvHLoad, pvE;
	malloc_16bit_aligned(&vProfile, m*segLen);
	malloc_16bit_aligned(&vProfile_tmp, m);
	malloc_16bit_aligned(&pvHStore, m*segLen);
	malloc_16bit_aligned(&pvHLoad, m*segLen);
	malloc_16bit_aligned(&pvE, m*segLen);

	__m128i vNegLimit = _mm_set1_epi16(MINUS_INF_SSE);
	const __m128i vGapO = _mm_set1_epi16((short)gapo);
    const __m128i vGapE = _mm_set1_epi16((short)gape);

	int segNum = 0;
	int index = 0;

	///n是字符集大小
	///band_length可以向上扩展到4的倍数
    for (k=0; k<m; ++k) 
	{
		const int8_t *p = &mat[k * m];
        ///sse vector个数
        for (i=0; i<segLen; ++i) 
		{
			///是个union
            __m128i_16_t t;
			j = i;
			for (segNum=0; segNum<segWidth; ++segNum) 
			{
				///注意这里是qlen, 而不是band_length
				t.v[segNum] = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
                j += segLen;
			}
			_mm_store_si128(&vProfile.mem[index], t.m);
            ++index;
		}
	}

	/* initialize H and E */
	__m128i_16_t h;
	__m128i_16_t e;
	for (segNum=0; segNum<segWidth; ++segNum) {
		h.v[segNum] = 0;
		e.v[segNum] = h.v[segNum] - gapoe;
	}
	index = 0;
	for (i=0; i<segLen; ++i) {		
		_mm_store_si128(&pvHStore.mem[index], h.m);
		_mm_store_si128(&pvE.mem[index], e.m);
		++index;
	}

	__m128i tmp_E;
	int8_t tmp_score;

	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		__m128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m128i vF = vNegLimit;

		__m128i vH = pvHStore.mem[0];

		/* Correct part of the vProfile */
        const __m128i* vP = vProfile.mem + seq_nt4_table[(int)target[i]] * segLen;

		/* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad.mem;
        pvHLoad.mem = pvHStore.mem;
        pvHStore.mem = pv;

		/**
		beg = i;
		end = i + band_length; // only loop through [beg,end) of the query sequence
		**/
		
		/* inner loop to process the query sequence */
        for (j=0; j<segLen; ++j) 
		{
			///vH[i, j] = vH[i-1, j-1] + vP[i, j]
            vH = _mm_add_epi16(vH, _mm_load_si128(vP + j));
			///vE[i, j]
            vE = _mm_load_si128(pvE.mem + j);
			/* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
			/* Save vH values. */
			///vH应该是一点都不要动的
            _mm_store_si128(pvHStore.mem + j, vH);


			/* Update vE value. */
            vH = _mm_sub_epi16(vH, vGapO);
            vE = _mm_sub_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            ///calculate vE[i+1, j] from vE[i, j] & vH[i, j]
			///vE的存储是要移位的啊
			///vE[j]是要存到vE[j-1]处的
            _mm_store_si128(pvE.mem + j, vE);

			/* Update vF value. */
            ///calculate vF[i, j+1] from vF[i, j]
            vF = _mm_sub_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);

			/* Load the next vH. */
            vH = _mm_load_si128(pvHLoad.mem + j);
		}

		///pvE和vP最好都要动一下
		///pvE好办，但是vP难办啊, 因为vP有m个
		for (k=0; k<m; ++k) 
		{
			_mm_store_si128(vProfile_tmp.mem+k, _mm_srli_si128(vProfile.mem[k * segLen], 2));
			j = segWidth*segLen+i;
			const int8_t *p = &mat[k * m];
			tmp_score = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
			vProfile_tmp.mem[k] = _mm_insert_epi16(vProfile_tmp.mem[k], tmp_score, 7);
		}
		
		tmp_E = _mm_srli_si128(pvE.mem[0], 2);
		tmp_E = _mm_insert_epi16(tmp_E, MINUS_INF_SSE, 7);

		__m128i* tmp_vP;
		for (j=1; j<segLen; ++j) 
		{
			_mm_store_si128(pvE.mem + j - 1, pvE.mem[j]);

			for (k=0; k<m; ++k) 
			{
				tmp_vP = vProfile.mem + k * segLen;
				_mm_store_si128(tmp_vP + j - 1, tmp_vP[j]);
			}
		}

		_mm_store_si128(pvE.mem + segLen-1, tmp_E);
		for (k=0; k<m; ++k) 
		{
			tmp_vP = vProfile.mem + k * segLen;
			_mm_store_si128(tmp_vP+segLen-1, vProfile_tmp.mem[k]);
		}

		///下面correct步不涉及到vP和vE, 所以应该没关系
		/* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        ///外循环是单个sse vector的长度
        ///不更新E会导致E偏小
		for (k=0; k<segWidth; ++k) 
		{
			///此时的vF是最后一个vF sse vector
            ///现在要用这个vF sse vector来矫正第一个vF sse vector
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, MINUS_INF_SSE, 0);
			///刚开始的时候是vF[0, j]
            ///内循环是sse vector的个数
            for (j=0; j<segLen; ++j) 
			{
				///vH[0, j]和vF[0, j]
                ///谁大谁就是真的
                vH = _mm_load_si128(pvHStore.mem + j);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore.mem + j, vH);

				vH = _mm_sub_epi16(vH, vGapO);
                vF = _mm_sub_epi16(vF, vGapE);

				if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;

			}
		}

	end:
		i;

	}

	int16_t *t = (int16_t*)pvHStore.mem;
	score = t[0];

	for (i = 0; i < band_length; i++)
	{
		if (t[i] > score)
		{
			score = t[i];
		}
	}
	





	free_16bit_aligned(&vProfile);
	free_16bit_aligned(&vProfile_tmp);
	free_16bit_aligned(&pvHStore);
	free_16bit_aligned(&pvHLoad);
	free_16bit_aligned(&pvE);


	qry->score = score;
	return score;
}

void output_stripe(__m128i* vector, int32_t segWidth, int32_t segLen, int32_t startIndice)
{
	int i, j;

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 0));
		startIndice++;
	}

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 1));
		startIndice++;
	}

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 2));
		startIndice++;
	}

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 3));
		startIndice++;
	}

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 4));
		startIndice++;
	}

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 5));
		startIndice++;
	}

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 6));
		startIndice++;
	}

	for (j = 0; j < segLen; j++)
	{
		fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector[j], 7));
		startIndice++;
	}
		
	
	
}


void output_stripe_single_word(__m128i vector, int32_t segLen, int32_t startIndice, int32_t word_j)
{

	startIndice+=word_j;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 0));
	startIndice+=segLen;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 1));
	startIndice+=segLen;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 2));
	startIndice+=segLen;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 3));
	startIndice+=segLen;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 4));
	startIndice+=segLen;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 5));
	startIndice+=segLen;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 6));
	startIndice+=segLen;
	fprintf(stderr, "j: %d, %d\n", startIndice, (int16_t)_mm_extract_epi16(vector, 7));
}


///qlen = tlen + 2*w
int ksw_semi_global_sse_debug_debug(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	int beg, end;
	int i, j, k = MINUS_INF_SSE, gapoe = gapo + gape, score;
	///可以把band_length扩展到能被4整除的数
	int band_length = 2*w+1;
	///单个sse element大小为16 bit, 所以sse vector中共有8个元素
    const int32_t segWidth = 8; /* number of values in vector unit */
    ///sse vector的个数
    const int32_t segLen = (band_length + segWidth - 1) / segWidth;
	aligned_sse vProfile, vProfile_tmp, pvHStore, pvHLoad, pvE, pvE_tmp;
	malloc_16bit_aligned(&vProfile, m*segLen);
	malloc_16bit_aligned(&vProfile_tmp, m);
	malloc_16bit_aligned(&pvHStore, m*segLen);
	malloc_16bit_aligned(&pvHLoad, m*segLen);
	malloc_16bit_aligned(&pvE, m*segLen);
	malloc_16bit_aligned(&pvE_tmp, m*segLen);
	

	__m128i vNegLimit = _mm_set1_epi16(MINUS_INF_SSE);
	const __m128i vGapOE = _mm_set1_epi16((short)gapoe);
    const __m128i vGapE = _mm_set1_epi16((short)gape);

	int segNum = 0;
	int index = 0;

	///n是字符集大小
	///band_length可以向上扩展到4的倍数
    for (k=0; k<m; ++k) 
	{
		const int8_t *p = &mat[k * m];
        ///sse vector个数
        for (i=0; i<segLen; ++i) 
		{
			///是个union
            __m128i_16_t t;
			j = i;
			for (segNum=0; segNum<segWidth; ++segNum) 
			{
				///注意这里是qlen, 而不是band_length
				t.v[segNum] = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
                j += segLen;
			}
			_mm_store_si128(&vProfile.mem[index], t.m);
            ++index;
		}
	}
	/* initialize H and E */
	__m128i_16_t h;
	__m128i_16_t e;
	for (segNum=0; segNum<segWidth; ++segNum) {
		h.v[segNum] = 0;
		e.v[segNum] = h.v[segNum] - gapoe;
		///e.v[segNum] = h.v[segNum] - gapo;
	}
	index = 0;
	for (i=0; i<segLen; ++i) {		
		_mm_store_si128(&pvHStore.mem[index], h.m);
		_mm_store_si128(&pvE.mem[index], e.m);
		++index;
	}

	__m128i tmp_E;
	int8_t tmp_score;

	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		__m128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m128i vF = vNegLimit;

		__m128i vH = pvHStore.mem[0];
		__m128i vHm;

		/* Correct part of the vProfile */
        //const __m128i* vP = vProfile.mem + seq_nt4_table[(int)target[i]] * segLen;
		__m128i* vP = vProfile.mem + seq_nt4_table[(int)target[i]] * segLen;

		/* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad.mem;
        pvHLoad.mem = pvHStore.mem;
        pvHStore.mem = pv;

		///if (i == 2 || i == 3)
		{
			fprintf(stderr, "*******i(h): %d*******\n", i);
			output_stripe(pvHLoad.mem, segWidth, segLen, i);
			fprintf(stderr, "*******i(p): %d*******\n", i);
			output_stripe(vP, segWidth, segLen, i);
			fprintf(stderr, "*******i(e): %d*******\n", i);
			output_stripe(pvE.mem, segWidth, segLen, i);
			fprintf(stderr, "segLen: %d, segWidth: %d\n", segLen, segWidth);
		}
		
		
		
		memcpy(pvE_tmp.mem, pvE.mem, segLen*sizeof(__m128i));

		/* inner loop to process the query sequence */
        for (j=0; j<segLen; ++j) 
		{
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Cells are computed in the following order:
			//   M(i,j)   = H(i-1,j-1) + S(i,j)
			//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
			// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
			// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
			// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
			// In practice, this should happen very rarely given a reasonable scoring system.

			///M(i,j) = H(i-1,j-1) + S(i,j)
			vHm = _mm_add_epi16(vH, _mm_load_si128(vP + j));
			///E(i, j)
			vE = _mm_load_si128(pvE.mem + j);
			///H(i,j) = max{M(i,j), E(i,j), F(i,j)}
			vH = _mm_max_epi16(vHm, vE);
			vH = _mm_max_epi16(vH, vF);
			_mm_store_si128(pvHStore.mem + j, vH);

			///M(i,j)-gapo-gape
			vHm = _mm_sub_epi16(vHm, vGapOE);
			///E(i,j)-gape
			vE = _mm_sub_epi16(vE, vGapE);
			///E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
			vE = _mm_max_epi16(vE, vHm);
			_mm_store_si128(pvE.mem + j, vE);

			// if (i == 2)
			// {
			// 	fprintf(stderr, "*******-i(f): %d, j: %d*******\n", i, j);
			// 	output_stripe_single_word(vF, segLen, i, j);
			// }


			///F(i,j) - gape
			vF = _mm_sub_epi16(vF, vGapE);
			///F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
            vF = _mm_max_epi16(vF, vHm);

			// if (i == 2)
			// {
			// 	fprintf(stderr, "*******+i(f): %d, j: %d*******\n", i, j+1);
			// 	output_stripe_single_word(vF, segLen, i, j+1);
			// }

			/* Load the next vH. */
            vH = _mm_load_si128(pvHLoad.mem + j);
		}

		///下面correct步不涉及到vP和vE, 所以应该没关系
		/* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
		//   M(i,j)   = H(i-1,j-1) + S(i,j)
		//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
		//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
		//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
        ///外循环是单个sse vector的长度
        ///不更新E会导致E偏小
		/**
		fprintf(stderr, "*******i(h)+: %d*******\n", i);
		output_stripe(pvHStore.mem, segWidth, segLen, i);
		**/
		for (k=0; k<segWidth; ++k) 
		{
			///此时的vF是最后一个vF sse vector
            ///现在要用这个vF sse vector来矫正第一个vF sse vector
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, MINUS_INF_SSE, 0);

			

			///刚开始的时候是vF[0, j]
            ///内循环是sse vector的个数
            for (j=0; j<segLen; ++j) 
			{
				/**
				if (i == 2)
				{
					fprintf(stderr, "*******new sbsb vF, i: %d, j: %d*******\n", i, j);
					output_stripe_single_word(vF, segLen, i, j);
				}
				**/
				
				
				/// H(i,j) = max{M(i,j), E(i,j), F(i,j)}
				///此时H(i,j)至少是max{M(i,j), E(i,j)}
				///所以只要和F(i,j)}比大小
				// vH = _mm_load_si128(pvHStore.mem + j);
				// if (i == 2)
				// {
				// 	fprintf(stderr, "*******new sbsb vH+, i: %d, j: %d*******\n", i, j);
				// 	output_stripe_single_word(vH, segLen, i, j);
				// }

                // vH = _mm_max_epi16(vH,vF);

				// if (i == 2)
				// {
				// 	fprintf(stderr, "*******new sbsb vH-, i: %d, j: %d*******\n", i, j);
				// 	output_stripe_single_word(vH, segLen, i, j);
				// }

                // _mm_store_si128(pvHStore.mem + j, vH);
				// ///M(i,j) = H(i-1,j-1) + S(i,j)
				// vHm = _mm_load_si128(pvHLoad.mem + j);
				// vH = _mm_add_epi16(vHm, _mm_load_si128(vP + j));


				
				///M(i,j) = H(i-1,j-1) + S(i,j)
				vHm = _mm_load_si128(pvHLoad.mem + j);
				vH = _mm_add_epi16(vHm, _mm_load_si128(vP + j));
				vE = _mm_load_si128(pvE_tmp.mem + j);
				/// H(i,j) = max{M(i,j), E(i,j), F(i,j)}
				vHm = _mm_max_epi16(vH,vF);
				vHm = _mm_max_epi16(vHm,vE);
				_mm_store_si128(pvHStore.mem + j, vHm);
				/**
				if (i == 2)
				{
					fprintf(stderr, "*******new sbsb vH+, i: %d, j: %d*******\n", i, j);
					output_stripe_single_word(vH, segLen, i, j);
				}
				**/
				


				///M(i,j) = H(i-1,j-1) + S(i,j)
				// vHm = _mm_load_si128(pvHLoad.mem + j);
				// vH = _mm_add_epi16(vHm, _mm_load_si128(vP + j));
				/// F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				vH = _mm_sub_epi16(vH, vGapOE);
				vF = _mm_sub_epi16(vF, vGapE);
				
				///_mm_cmpgt_epi16(vF, vH): 如果vF > vH, 则返回结果全是1
				//_mm_movemask_epi8等于就是提取元素了呗
				///只要vH里所有元素都>=vF, 就跳出循环
				if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) 
				{
					/**
					vF = _mm_max_epi16(vF, vH);
					goto end;
					**/
					///fprintf(stderr, "haha\n");
				}
			
				vF = _mm_max_epi16(vF, vH);
			}
		}

	end:

		

		/**
		fprintf(stderr, "*******i(h)-: %d*******\n", i);
		output_stripe(pvHStore.mem, segWidth, segLen, i);
		**/
		///pvE和vP最好都要动一下
		///pvE好办，但是vP难办啊, 因为vP有m个
		for (k=0; k<m; ++k) 
		{
			_mm_store_si128(vProfile_tmp.mem+k, _mm_srli_si128(vProfile.mem[k * segLen], 2));
			j = segWidth*segLen+i;
			const int8_t *p = &mat[k * m];
			tmp_score = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
			vProfile_tmp.mem[k] = _mm_insert_epi16(vProfile_tmp.mem[k], tmp_score, 7);
		}
		
		tmp_E = _mm_srli_si128(pvE.mem[0], 2);
		tmp_E = _mm_insert_epi16(tmp_E, MINUS_INF_SSE, 7);

		__m128i* tmp_vP;
		for (j=1; j<segLen; ++j) 
		{
			_mm_store_si128(pvE.mem + j - 1, pvE.mem[j]);

			for (k=0; k<m; ++k) 
			{
				tmp_vP = vProfile.mem + k * segLen;
				_mm_store_si128(tmp_vP + j - 1, tmp_vP[j]);
			}
		}

		_mm_store_si128(pvE.mem + segLen-1, tmp_E);
		for (k=0; k<m; ++k) 
		{
			tmp_vP = vProfile.mem + k * segLen;
			_mm_store_si128(tmp_vP+segLen-1, vProfile_tmp.mem[k]);
		}



		
	}

	
	fprintf(stderr, "*******i(h): %d*******\n", i);
	output_stripe(pvHStore.mem, segWidth, segLen, i);


	int16_t *t = (int16_t*)pvHStore.mem;
	score = t[0];

	for (i = 0; i < band_length; i++)
	{
		if (t[i] > score)
		{
			score = t[i];
		}
	}
	





	free_16bit_aligned(&vProfile);
	free_16bit_aligned(&vProfile_tmp);
	free_16bit_aligned(&pvHStore);
	free_16bit_aligned(&pvHLoad);
	free_16bit_aligned(&pvE);
	free_16bit_aligned(&pvE_tmp);


	qry->score = score;
	return score;
}

///qlen = tlen + 2*w
int ksw_semi_global_sse_debug(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	int beg, end;
	int i, j, k = MINUS_INF_SSE, gapoe = gapo + gape, score;
	///可以把band_length扩展到能被4整除的数
	int band_length = 2*w+1;
	///单个sse element大小为16 bit, 所以sse vector中共有8个元素
    const int32_t segWidth = 8; /* number of values in vector unit */
    ///sse vector的个数
    const int32_t segLen = (band_length + segWidth - 1) / segWidth;
	aligned_sse vProfile, vProfile_tmp, pvHStore, pvHLoad, pvE, pvE_tmp;
	malloc_16bit_aligned(&vProfile, m*segLen);
	malloc_16bit_aligned(&vProfile_tmp, m);
	malloc_16bit_aligned(&pvHStore, m*segLen);
	malloc_16bit_aligned(&pvHLoad, m*segLen);
	malloc_16bit_aligned(&pvE, m*segLen);
	malloc_16bit_aligned(&pvE_tmp, m*segLen);
	

	__m128i vNegLimit = _mm_set1_epi16(MINUS_INF_SSE);
	const __m128i vGapOE = _mm_set1_epi16((short)gapoe);
    const __m128i vGapE = _mm_set1_epi16((short)gape);

	int segNum = 0;
	int index = 0;

	///n是字符集大小
	///band_length可以向上扩展到4的倍数
    for (k=0; k<m; ++k) 
	{
		const int8_t *p = &mat[k * m];
        ///sse vector个数
        for (i=0; i<segLen; ++i) 
		{
			///是个union
            __m128i_16_t t;
			j = i;
			for (segNum=0; segNum<segWidth; ++segNum) 
			{
				///注意这里是qlen, 而不是band_length
				t.v[segNum] = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
                j += segLen;
			}
			_mm_store_si128(&vProfile.mem[index], t.m);
            ++index;
		}
	}
	/* initialize H and E */
	__m128i_16_t h;
	__m128i_16_t e;
	for (segNum=0; segNum<segWidth; ++segNum) {
		h.v[segNum] = 0;
		e.v[segNum] = h.v[segNum] - gapoe;
		///e.v[segNum] = h.v[segNum] - gapo;
	}
	index = 0;
	for (i=0; i<segLen; ++i) {		
		_mm_store_si128(&pvHStore.mem[index], h.m);
		_mm_store_si128(&pvE.mem[index], e.m);
		++index;
	}

	__m128i tmp_E;
	int8_t tmp_score;

	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		__m128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m128i vF = vNegLimit;

		__m128i vH = pvHStore.mem[0];
		__m128i vHm;

		/* Correct part of the vProfile */
        //const __m128i* vP = vProfile.mem + seq_nt4_table[(int)target[i]] * segLen;
		__m128i* vP = vProfile.mem + seq_nt4_table[(int)target[i]] * segLen;

		/* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad.mem;
        pvHLoad.mem = pvHStore.mem;
        pvHStore.mem = pv;

		// if (i == 149)
		// {
		// 	fprintf(stderr, "*******i(h): %d*******\n", i);
		// 	output_stripe(pvHLoad.mem, segWidth, segLen, i);
		// 	fprintf(stderr, "*******i(p): %d*******\n", i);
		// 	output_stripe(vP, segWidth, segLen, i);
		// 	fprintf(stderr, "*******i(e): %d*******\n", i);
		// 	output_stripe(pvE.mem, segWidth, segLen, i);
		// 	fprintf(stderr, "segLen: %d, segWidth: %d\n", segLen, segWidth);
		// }
		
		
		
		memcpy(pvE_tmp.mem, pvE.mem, segLen*sizeof(__m128i));

		/* inner loop to process the query sequence */
        for (j=0; j<segLen; ++j) 
		{
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Cells are computed in the following order:
			//   M(i,j)   = H(i-1,j-1) + S(i,j)
			//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
			// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
			// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
			// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
			// In practice, this should happen very rarely given a reasonable scoring system.

			///M(i,j) = H(i-1,j-1) + S(i,j)
			vHm = _mm_add_epi16(vH, _mm_load_si128(vP + j));
			///E(i, j)
			vE = _mm_load_si128(pvE.mem + j);
			///H(i,j) = max{M(i,j), E(i,j), F(i,j)}
			vH = _mm_max_epi16(vHm, vE);
			vH = _mm_max_epi16(vH, vF);
			_mm_store_si128(pvHStore.mem + j, vH);

			///M(i,j)-gapo-gape
			vHm = _mm_sub_epi16(vHm, vGapOE);
			///E(i,j)-gape
			vE = _mm_sub_epi16(vE, vGapE);
			///E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
			vE = _mm_max_epi16(vE, vHm);
			_mm_store_si128(pvE.mem + j, vE);

			// if (i == 2)
			// {
			// 	fprintf(stderr, "*******-i(f): %d, j: %d*******\n", i, j);
			// 	output_stripe_single_word(vF, segLen, i, j);
			// }


			///F(i,j) - gape
			vF = _mm_sub_epi16(vF, vGapE);
			///F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
            vF = _mm_max_epi16(vF, vHm);

			// if (i == 2)
			// {
			// 	fprintf(stderr, "*******+i(f): %d, j: %d*******\n", i, j+1);
			// 	output_stripe_single_word(vF, segLen, i, j+1);
			// }

			/* Load the next vH. */
            vH = _mm_load_si128(pvHLoad.mem + j);
		}

		///下面correct步不涉及到vP和vE, 所以应该没关系
		/* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
		//   M(i,j)   = H(i-1,j-1) + S(i,j)
		//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
		//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
		//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
        ///外循环是单个sse vector的长度
        ///不更新E会导致E偏小
		/**
		fprintf(stderr, "*******i(h)+: %d*******\n", i);
		output_stripe(pvHStore.mem, segWidth, segLen, i);
		**/
		for (k=0; k<segWidth; ++k) 
		{
			///此时的vF是最后一个vF sse vector
            ///现在要用这个vF sse vector来矫正第一个vF sse vector
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, MINUS_INF_SSE, 0);

			

			///刚开始的时候是vF[0, j]
            ///内循环是sse vector的个数
            for (j=0; j<segLen; ++j) 
			{

				// if (i == 2)
				// {
				// 	fprintf(stderr, "*******new vF, i: %d, j: %d*******\n", i, j);
				// 	output_stripe_single_word(vF, segLen, i, j);
				// }
				
				/**
				/// H(i,j) = max{M(i,j), E(i,j), F(i,j)}
				///此时H(i,j)至少是max{M(i,j), E(i,j)}
				///所以只要和F(i,j)}比大小
				vH = _mm_load_si128(pvHStore.mem + j);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore.mem + j, vH);
				**/

				///M(i,j) = H(i-1,j-1) + S(i,j)
				vHm = _mm_load_si128(pvHLoad.mem + j);
				vH = _mm_add_epi16(vHm, _mm_load_si128(vP + j));

				vE = _mm_load_si128(pvE_tmp.mem + j);
				/// H(i,j) = max{M(i,j), E(i,j), F(i,j)}
				vHm = _mm_max_epi16(vH,vF);
				vHm = _mm_max_epi16(vHm,vE);
				_mm_store_si128(pvHStore.mem + j, vHm);


				/// F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				vH = _mm_sub_epi16(vH, vGapOE);
				vF = _mm_sub_epi16(vF, vGapE);
				
				///_mm_cmpgt_epi16(vF, vH): 如果vF > vH, 则返回结果全是1
				//_mm_movemask_epi8等于就是提取元素了呗
				///只要vH里所有元素都>=vF, 就跳出循环
				if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) 
				{
					goto end;
					///fprintf(stderr, "haha\n");
				}
			
				vF = _mm_max_epi16(vF, vH);
			}
		}

	end:
		/**
		fprintf(stderr, "*******i(h)-: %d*******\n", i);
		output_stripe(pvHStore.mem, segWidth, segLen, i);
		**/
		///pvE和vP最好都要动一下
		///pvE好办，但是vP难办啊, 因为vP有m个
		for (k=0; k<m; ++k) 
		{
			_mm_store_si128(vProfile_tmp.mem+k, _mm_srli_si128(vProfile.mem[k * segLen], 2));
			j = segWidth*segLen+i;
			const int8_t *p = &mat[k * m];
			tmp_score = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
			vProfile_tmp.mem[k] = _mm_insert_epi16(vProfile_tmp.mem[k], tmp_score, 7);
		}
		
		tmp_E = _mm_srli_si128(pvE.mem[0], 2);
		tmp_E = _mm_insert_epi16(tmp_E, MINUS_INF_SSE, 7);

		__m128i* tmp_vP;
		for (j=1; j<segLen; ++j) 
		{
			_mm_store_si128(pvE.mem + j - 1, pvE.mem[j]);

			for (k=0; k<m; ++k) 
			{
				tmp_vP = vProfile.mem + k * segLen;
				_mm_store_si128(tmp_vP + j - 1, tmp_vP[j]);
			}
		}

		_mm_store_si128(pvE.mem + segLen-1, tmp_E);
		for (k=0; k<m; ++k) 
		{
			tmp_vP = vProfile.mem + k * segLen;
			_mm_store_si128(tmp_vP+segLen-1, vProfile_tmp.mem[k]);
		}



		
	}

	
	// fprintf(stderr, "*******i(h): %d*******\n", i);
	// output_stripe(pvHStore.mem, segWidth, segLen, i);


	int16_t *t = (int16_t*)pvHStore.mem;
	score = t[0];

	for (i = 0; i < band_length; i++)
	{
		if (t[i] > score)
		{
			score = t[i];
		}
	}
	





	free_16bit_aligned(&vProfile);
	free_16bit_aligned(&vProfile_tmp);
	free_16bit_aligned(&pvHStore);
	free_16bit_aligned(&pvHLoad);
	free_16bit_aligned(&pvE);
	free_16bit_aligned(&pvE_tmp);


	qry->score = score;
	return score;
}







int ksw_semi_global_sse_new(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_, int flag)
{
	int beg, end;
	int i, j, k = MINUS_INF_SSE, gapoe = gapo + gape, score;
	///可以把band_length扩展到能被4整除的数
	int band_length = 2*w+1;
	///单个sse element大小为16 bit, 所以sse vector中共有8个元素
    const int32_t segWidth = 8; /* number of values in vector unit */
    ///sse vector的个数
    const int32_t segLen = (band_length + segWidth - 1) / segWidth;
	aligned_sse vProfile, vProfile_tmp, pvHStore, pvHLoad, pvE;
	malloc_16bit_aligned(&vProfile, m*segLen);
	malloc_16bit_aligned(&vProfile_tmp, m);
	malloc_16bit_aligned(&pvHStore, m*segLen);
	malloc_16bit_aligned(&pvHLoad, m*segLen);
	malloc_16bit_aligned(&pvE, m*segLen);
	

	__m128i vNegLimit = _mm_set1_epi16(MINUS_INF_SSE);
	const __m128i vGapOE = _mm_set1_epi16((short)gapoe);
    const __m128i vGapE = _mm_set1_epi16((short)gape);

	int segNum = 0;
	int index = 0;

	///n是字符集大小
	///band_length可以向上扩展到4的倍数
    for (k=0; k<m; ++k) 
	{
		const int8_t *p = &mat[k * m];
        ///sse vector个数
        for (i=0; i<segLen; ++i) 
		{
			///是个union
            __m128i_16_t t;
			j = i;
			for (segNum=0; segNum<segWidth; ++segNum) 
			{
				///注意这里是qlen, 而不是band_length
				t.v[segNum] = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
                j += segLen;
			}
			_mm_store_si128(&vProfile.mem[index], t.m);
            ++index;
		}
	}
	/* initialize H and E */
	__m128i_16_t h;
	__m128i_16_t e;
	for (segNum=0; segNum<segWidth; ++segNum) {
		h.v[segNum] = 0;
		e.v[segNum] = h.v[segNum] - gapoe;
	}
	index = 0;
	for (i=0; i<segLen; ++i) {		
		_mm_store_si128(&pvHStore.mem[index], h.m);
		_mm_store_si128(&pvE.mem[index], e.m);
		++index;
	}

	__m128i tmp_E;
	int8_t tmp_score;

	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		__m128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m128i vF = vNegLimit;

		__m128i vH = pvHStore.mem[0];

		/* Correct part of the vProfile */
        //const __m128i* vP = vProfile.mem + seq_nt4_table[(int)target[i]] * segLen;
		__m128i* vP = vProfile.mem + seq_nt4_table[(int)target[i]] * segLen;

		/* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad.mem;
        pvHLoad.mem = pvHStore.mem;
        pvHStore.mem = pv;

		if(flag)
		{
			
			fprintf(stderr, "*******i(p): %d*******\n", i);
			output_stripe(vP, segWidth, segLen, i);
			

			fprintf(stderr, "0 *******i(h): %d*******\n", i);
			output_stripe(pvHLoad.mem, segWidth, segLen, i);

			
			fprintf(stderr, "*******i(e): %d*******\n", i);
			output_stripe(pvE.mem, segWidth, segLen, i);
			
		}
		
		
		/* inner loop to process the query sequence */
        for (j=0; j<segLen; ++j) 
		{
			///H'(i,j) = H(i-1,j-1) + S(i,j)
			vH = _mm_add_epi16(vH, _mm_load_si128(vP + j));
			///E(i, j)
			vE = _mm_load_si128(pvE.mem + j);
			///H'(i,j) = max{H'(i,j), E(i,j), F(i,j)}
			vH = _mm_max_epi16(vH, vE);
			vH = _mm_max_epi16(vH, vF);
			_mm_store_si128(pvHStore.mem + j, vH);

			///H'(i,j)-gapoe
			vH = _mm_sub_epi16(vH, vGapOE);
			///E(i,j)-gape
			vE = _mm_sub_epi16(vE, vGapE);
			///E(i+1,j) = max{H'(i,j)-gapoe, E(i,j)- gape}
			vE = _mm_max_epi16(vE, vH);
			_mm_store_si128(pvE.mem + j, vE);

			///F(i,j) - gape
			vF = _mm_sub_epi16(vF, vGapE);
			///F(i,j+1) = max{H'(i,j)-gapoe, F(i,j)-gape}
            vF = _mm_max_epi16(vF, vH);

			/* Load the next vH. */
            vH = _mm_load_si128(pvHLoad.mem + j);
		}


		if(flag)
		{
			/**
			fprintf(stderr, "*******i(p): %d*******\n", i);
			output_stripe(vP, segWidth, segLen, i);
			**/

			fprintf(stderr, "1 *******i(h): %d*******\n", i);
			output_stripe(pvHStore.mem, segWidth, segLen, i);

			/**
			fprintf(stderr, "*******i(e): %d*******\n", i);
			output_stripe(pvE.mem, segWidth, segLen, i);
			**/
		}


		/* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
		/* 
		1. E[i+1, j] = max{H[i, j] - vGapOE, E[i, j] - vGapE}
		2. F[i, j+1] = max{H[i, j] - vGapOE, F[i, j] - vGapE}
		3. H[i, j] = max{H[i-1, j-1] + Q[i, j], E[i, j], F[i, j]}
		
		(1) for #2, if H[i, j] - vGapOE >= F[i, j] - vGapE, new F[i, j] will not influence the next H

		(2) F[i, j+1] is always equal to  F[i, j] - vGapE. This is because if F[i, j+1] == H[i, j] - vGapOE,
		the H[i, j] - vGapOE > F[i, j] - vGapE. In this case, the calculation would be terminated.
		
		(3) why the terminated condition is max{vH, vF} - vGapOE >= F[i, j] - vGapE, instead of  vH - vGapOE >= F[i, j] - vGapE?
		the reason is that if vF > vH, then the terminated condition is vF - vGapOE >= vF - vGapE. This cannot be true, 
		since vGapOE > vGapE.
		*/
		for (k=0; k<segWidth; ++k) 
		{
			///此时的vF是最后一个vF sse vector
            ///现在要用这个vF sse vector来矫正第一个vF sse vector
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, MINUS_INF_SSE, 0);

			///the first vF is vF[0]
            for (j=0; j<segLen; ++j) 
			{
				vH = _mm_load_si128(pvHStore.mem + j);
				vH = _mm_max_epi16(vH,vF);
				_mm_store_si128(pvHStore.mem + j, vH);

				vE = _mm_sub_epi16(vH, vGapOE);
				vE = _mm_max_epi16(vE, pvE.mem[j]);
				_mm_store_si128(pvE.mem + j, vE);

				/// F(i,j+1) = max{M(i,j)-gapoe, F(i,j) - gape} 
				vH = _mm_sub_epi16(vH, vGapOE);
				vF = _mm_sub_epi16(vF, vGapE);

				if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) 
				{
					goto end;
				}
			}
		}

	end:


		if(flag)
		{
			/**
			fprintf(stderr, "*******i(p): %d*******\n", i);
			output_stripe(vP, segWidth, segLen, i);
			**/

			fprintf(stderr, "2 *******i(h): %d*******\n", i);
			output_stripe(pvHStore.mem, segWidth, segLen, i);

			/**
			fprintf(stderr, "*******i(e): %d*******\n", i);
			output_stripe(pvE.mem, segWidth, segLen, i);
			**/
		}

		///shift pvE & vProfile
		for (k=0; k<m; ++k) 
		{
			_mm_store_si128(vProfile_tmp.mem+k, _mm_srli_si128(vProfile.mem[k * segLen], 2));
			j = segWidth*segLen+i;
			const int8_t *p = &mat[k * m];
			tmp_score = j >= qlen ? 0 : p[seq_nt4_table[(int)query[j]]];
			vProfile_tmp.mem[k] = _mm_insert_epi16(vProfile_tmp.mem[k], tmp_score, 7);
		}
		
		tmp_E = _mm_srli_si128(pvE.mem[0], 2);
		tmp_E = _mm_insert_epi16(tmp_E, MINUS_INF_SSE, 7);

		__m128i* tmp_vP;
		for (j=1; j<segLen; ++j) 
		{
			_mm_store_si128(pvE.mem + j - 1, pvE.mem[j]);

			for (k=0; k<m; ++k) 
			{
				tmp_vP = vProfile.mem + k * segLen;
				_mm_store_si128(tmp_vP + j - 1, tmp_vP[j]);
			}
		}

		_mm_store_si128(pvE.mem + segLen-1, tmp_E);
		for (k=0; k<m; ++k) 
		{
			tmp_vP = vProfile.mem + k * segLen;
			_mm_store_si128(tmp_vP+segLen-1, vProfile_tmp.mem[k]);
		}



		
	}

	
	// fprintf(stderr, "*******i(h): %d*******\n", i);
	// output_stripe(pvHStore.mem, segWidth, segLen, i);


	int16_t *t = (int16_t*)pvHStore.mem;
	score = t[0];

	for (i = 0; i < band_length; i++)
	{
		if (t[i] > score)
		{
			score = t[i];
		}
	}
	





	free_16bit_aligned(&vProfile);
	free_16bit_aligned(&vProfile_tmp);
	free_16bit_aligned(&pvHStore);
	free_16bit_aligned(&pvHLoad);
	free_16bit_aligned(&pvE);


	qry->score = score;
	return score;
}











































///qlen = tlen + 2*w
int ksw_semi_global_back(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	///这个就是一列的值
	eh_t *eh;
	///列，也就是query的profile
	int8_t *qp; // query profile
	int i, j, k = MINUS_INF, gapoe = gapo + gape, score;
	int band_length = 2*w+1;
	int32_t beg, end;


	///这些都是求cigar用的
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	int n_col = band_length; // maximum #columns of the backtrack matrix
	z = n_cigar_ && cigar_? (uint8_t *)malloc((long)n_col * tlen) : 0;

	///m是字符集大小
	qp = (int8_t *)malloc(qlen * m);
	///第一列的值
	eh = (eh_t *)calloc(qlen + 1, 8);

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		//for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
		for (j = 0; j < qlen; ++j) qp[i++] = p[seq_nt4_table[(int)query[j]]];
	}


	// fill the first row
	// eh[0].h = 0; eh[0].e = MINUS_INF;
	// for (j = 1; j <= qlen && j <= w; ++j)
	// 	eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	// 我觉得最后还是要看下到底是<还是<=
	for (j = 0; j < band_length; ++j)
	{
		eh[j].h = 0;
		///感觉应该是这个东西
		eh[j].e = eh[j].h - gapoe;
		///eh[j].e = eh[j].h - gapo;
	}
	for (; j <= qlen; ++j) 
	{
		eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	}

	/// 0 <= i <= (tlen - 1)
	/// 0 <= j <= qlen
	///所以i是从矩阵的第1列而不是第0列开始
	///但是j是从第0行开始的
	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		int32_t f = MINUS_INF, h1, t;
		int8_t *q = &qp[seq_nt4_table[(int)target[i]] * qlen];
		beg = i;
		end = i + band_length; // only loop through [beg,end) of the query sequence
		///h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		h1 = MINUS_INF;

		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				m += q[j];
				///match/mismatch: d = 0
				///e: d = 1
				///f: d = 2
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				///这个是区分e是gap open 还是gap extension
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;

				t = m - gapoe;
				f -= gape;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				m += q[j];
				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				e  = e > t? e : t;
				p->e = e;
				t = m - gapoe;
				f -= gape;
				f  = f > t? f : t;
			}
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}

	///score = eh[qlen].h;
	int max_i = tlen + w;
	score = eh[max_i].h;
	for (i = end; i > beg; i--)
	{
		if (eh[i].h > score)
		{
			score = eh[i].h;
			max_i = i;
		}
	}

	qry->score = score;
	qry->qe = max_i - 1;

	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;

		// (i,k) points to the last cell
		i = tlen - 1; 
		//k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1;
		k = max_i - 1;
		///k = max_i;

		while (i >= 0 && k >= 0) {
			///which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			which = z[(long)i * n_col + (k - i)] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		///if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;

		qry->qb = k + 1;
		qry->NM = 0;

		char hah[1000] = {0};
		int op, opl;
		int query_s = qry->qb;
		int target_s = 0;
		for (i = 0; i < *n_cigar_; ++i) // print CIGAR
		{
			op = (*cigar_)[i]&0xf;
			opl = (*cigar_)[i]>>4;

			sprintf(hah + strlen(hah), "%d%c", opl, "MDISH"[op]);

			if (op == 0)
			{
				for (k = 0; k < opl; k++)
				{
					if(query[query_s] != target[target_s] && 
					!(query[query_s] == 'C' && target[target_s] == 'T'))
					{
						qry->NM++;
					}
					query_s++;
					target_s++;
				}
			}
			else if (op == 1)
			{
				query_s += opl;
				qry->NM += opl;
			}
			else
			{
				target_s += opl;
				qry->NM += opl;
			}
		}

		/**
		if (query_s != qry->qe + 1)
		{
			fprintf(stderr, "query_s: %d, qry->qe: %d, cigar: %s\n", query_s, qry->qe, hah);
		}
		**/
	}



	


	

	free(eh); 
	free(qp); 
	free(z);
	return score;
}







///avoid malloc
///if we just have one candidate location, directly use this method instead of bit-vector
///and let the size of alphat to be fixed (option m)
int ksw_semi_global(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	///这个就是一列的值
	eh_t *eh;
	///列，也就是query的profile
	int8_t *qp; // query profile
	int i, j, k = MINUS_INF, gapoe = gapo + gape, score;
	int band_length = 2*w+1;
	int32_t beg, end;


	///这些都是求cigar用的
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	int n_col = band_length; // maximum #columns of the backtrack matrix
	z = n_cigar_ && cigar_? (uint8_t *)malloc((long)n_col * tlen) : 0;

	///m是字符集大小
	qp = (int8_t *)malloc(qlen * m);
	///第一列的值
	eh = (eh_t *)calloc(qlen + 1, 8);

	
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		//for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
		for (j = 0; j < qlen; ++j) qp[i++] = p[seq_nt4_table[(int)query[j]]];
	}
	/**
	int query_c;
	for (j = 0; j < qlen; ++j)
	{
		query_c = seq_nt4_table[(int)query[j]];
		i = j;
		for (k = 0; k < m; ++k)
		{
			///qp[k*qlen + j] = mat[k * m + query_c];
			qp[i] = mat[query_c];
			i += qlen;
			query_c += m;
		}

	}
	**/


	// fill the first row
	// eh[0].h = 0; eh[0].e = MINUS_INF;
	// for (j = 1; j <= qlen && j <= w; ++j)
	// 	eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	// 我觉得最后还是要看下到底是<还是<=
	for (j = 0; j < band_length; ++j)
	{
		eh[j].h = 0;
		///感觉应该是这个东西
		eh[j].e = eh[j].h - gapoe;
		///eh[j].e = eh[j].h - gapo;
	}
	for (; j <= qlen; ++j) 
	{
		eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	}

	/// 0 <= i <= (tlen - 1)
	/// 0 <= j <= qlen
	///所以i是从矩阵的第1列而不是第0列开始
	///但是j是从第0行开始的
	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		int32_t f = MINUS_INF, h1, t;
		int8_t *q = &qp[seq_nt4_table[(int)target[i]] * qlen];
		beg = i;
		end = i + band_length; // only loop through [beg,end) of the query sequence
		///h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		h1 = MINUS_INF;

		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				m += q[j];
				///match/mismatch: d = 0
				///e: d = 1
				///f: d = 2
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				///这个是区分e是gap open 还是gap extension
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;

				t = m - gapoe;
				f -= gape;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				m += q[j];
				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				e  = e > t? e : t;
				p->e = e;
				t = m - gapoe;
				f -= gape;
				f  = f > t? f : t;
			}
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}

	///score = eh[qlen].h;
	int max_i = tlen + w;
	score = eh[max_i].h;
	for (i = end; i > beg; i--)
	{
		if (eh[i].h > score)
		{
			score = eh[i].h;
			max_i = i;
		}
	}

	qry->score = score;
	qry->qe = max_i - 1;

	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;

		// (i,k) points to the last cell
		i = tlen - 1; 
		//k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1;
		k = max_i - 1;
		///k = max_i;

		while (i >= 0 && k >= 0) {
			///which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			which = z[(long)i * n_col + (k - i)] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		///if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;

		qry->qb = k + 1;
		/**
		qry->NM = 0;

		char hah[1000] = {0};
		int op, opl;
		int query_s = qry->qb;
		int target_s = 0;
		for (i = 0; i < *n_cigar_; ++i) // print CIGAR
		{
			op = (*cigar_)[i]&0xf;
			opl = (*cigar_)[i]>>4;

			sprintf(hah + strlen(hah), "%d%c", opl, "MDISH"[op]);

			if (op == 0)
			{
				for (k = 0; k < opl; k++)
				{
					if(query[query_s] != target[target_s] && 
					!(query[query_s] == 'C' && target[target_s] == 'T'))
					{
						qry->NM++;
					}
					query_s++;
					target_s++;
				}
			}
			else if (op == 1)
			{
				query_s += opl;
				qry->NM += opl;
			}
			else
			{
				target_s += opl;
				qry->NM += opl;
			}
		}

		
		// if (query_s != qry->qe + 1)
		// {
		// 	fprintf(stderr, "query_s: %d, qry->qe: %d, cigar: %s\n", query_s, qry->qe, hah);
		// }
		**/
		
	}



	


	

	free(eh); 
	free(qp); 
	free(z);
	return score;
}


inline int PenaltyDiffByQuality(int diff, int quality_base, int Q)
{
	if(diff == 0)
	{
		return 0;
	}
	double Phred;
	Phred = Q - quality_base;
	if (Phred > 40)
	{
		Phred = 40;
	}
	Phred = Phred/40;
	int return_Q = Phred*(diff);
	return return_Q;
}

///avoid malloc
///if we just have one candidate location, directly use this method instead of bit-vector
///and let the size of alphat to be fixed (option m)
int ksw_semi_global_quality_back(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat,
const int8_t *mat_diff, int gapo, int gape, int w, char* qulity, int quality_base, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	///这个就是一列的值
	eh_t *eh;
	///列，也就是query的profile
	int8_t *qp; // query profile
	int8_t *qp_diff; // query profile
	int i, j, k = MINUS_INF, gapoe = gapo + gape, score;
	int band_length = 2*w+1;
	int32_t beg, end;
	double Phred;


	///这些都是求cigar用的
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	int n_col = band_length; // maximum #columns of the backtrack matrix
	z = n_cigar_ && cigar_? (uint8_t *)malloc((long)n_col * tlen) : 0;

	///m是字符集大小
	///qp = (int8_t *)malloc(qlen * m);
	qp = (int8_t *)malloc(qlen * m * 2);
	qp_diff = qp + qlen * m;
	///第一列的值
	eh = (eh_t *)calloc(qlen + 1, 8);

	
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		const int8_t *p_diff = &mat_diff[k * m];
		for (j = 0; j < qlen; ++j) 
		{
			qp_diff[i] = p_diff[seq_nt4_table[(int)query[j]]];
			qp[i++] = p[seq_nt4_table[(int)query[j]]];
		}
	}


	// fill the first row
	// eh[0].h = 0; eh[0].e = MINUS_INF;
	// for (j = 1; j <= qlen && j <= w; ++j)
	// 	eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	// 我觉得最后还是要看下到底是<还是<=
	for (j = 0; j < band_length; ++j)
	{
		eh[j].h = 0;
		///感觉应该是这个东西
		eh[j].e = eh[j].h - gapoe;
		///eh[j].e = eh[j].h - gapo;
	}
	for (; j <= qlen; ++j) 
	{
		eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	}

	/// 0 <= i <= (tlen - 1)
	/// 0 <= j <= qlen
	///所以i是从矩阵的第1列而不是第0列开始
	///但是j是从第0行开始的
	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		int32_t f = MINUS_INF, h1, t;
		int8_t *q = &qp[seq_nt4_table[(int)target[i]] * qlen];
		int8_t *q_diff = &qp_diff[seq_nt4_table[(int)target[i]] * qlen];


		Phred = qulity[i] - quality_base;
		if (Phred > 40)
		{
			Phred = 40;
		}
		Phred = Phred/40;
		


		beg = i;
		end = i + band_length; // only loop through [beg,end) of the query sequence
		///h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		h1 = MINUS_INF;

		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				///m += q[j];
				///m = m + q[j] - q_diff[j];
				m = m + q[j] - (int)(q_diff[j] * Phred);


				///match/mismatch: d = 0
				///e: d = 1
				///f: d = 2
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				///这个是区分e是gap open 还是gap extension
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;

				t = m - gapoe;
				f -= gape;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				///m += q[j];
				///m = m + q[j] - q_diff[j];
				m = m + q[j] - (int)(q_diff[j] * Phred);



				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				e  = e > t? e : t;
				p->e = e;
				t = m - gapoe;
				f -= gape;
				f  = f > t? f : t;
			}
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}

	///score = eh[qlen].h;
	int max_i = tlen + w;
	score = eh[max_i].h;
	for (i = end; i > beg; i--)
	{
		if (eh[i].h > score)
		{
			score = eh[i].h;
			max_i = i;
		}
	}

	qry->score = score;
	qry->qe = max_i - 1;

	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;

		// (i,k) points to the last cell
		i = tlen - 1; 
		//k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1;
		k = max_i - 1;
		///k = max_i;

		while (i >= 0 && k >= 0) {
			///which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			which = z[(long)i * n_col + (k - i)] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		///if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;

		qry->qb = k + 1;		
	}

	free(eh); 
	free(qp); 
	free(z);
	return score;
}

inline void init_data(kswr_t* qry, int *n_cigar_, uint32_t **cigar_, int n_col, int q_len, int t_len, int m)
{
	///for eh
	if(qry->eh_size < (q_len + 1))
	{
		qry->eh_size = (q_len + 1);
		qry->eh = (eh_t *)realloc(qry->eh, sizeof(eh_t)*qry->eh_size);
	}
	memset(qry->eh, 0, sizeof(eh_t)*qry->eh_size);

	///for qp and qp_diff
	if(qry->qp_size < q_len * m * 2)
	{
		qry->qp_size = q_len * m * 2;
		qry->qp = (int8_t *)realloc(qry->qp, qry->qp_size);
	}
	qry->qp_diff = qry->qp + q_len * m;

	///for z
	if(n_cigar_ && cigar_ && qry->z_size < n_col * t_len)
	{
		qry->z_size = n_col * t_len;
		qry->z = (uint8_t *)realloc(qry->z, (long)qry->z_size);
	}
}
///avoid malloc
///if we just have one candidate location, directly use this method instead of bit-vector
///and let the size of alphat to be fixed (option m)
int ksw_semi_global_quality(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat,
const int8_t *mat_diff, int gapo, int gape, int w, char* qulity, int quality_base, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	///这个就是一列的值
	///eh_t *eh;
	///列，也就是query的profile
	///int8_t *qp; // query profile
	///int8_t *qp_diff; // query profile
	int i, j, k = MINUS_INF, gapoe = gapo + gape, score;
	int band_length = 2*w+1;
	int32_t beg, end;
	double Phred;


	///这些都是求cigar用的
	///uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	int n_col = band_length; // maximum #columns of the backtrack matrix
	///z = n_cigar_ && cigar_? (uint8_t *)malloc((long)n_col * tlen) : 0;

	///m是字符集大小
	///qp = (int8_t *)malloc(qlen * m);
	///qp = (int8_t *)malloc(qlen * m * 2);
	///qp_diff = qp + qlen * m;
	///第一列的值
	///eh = (eh_t *)calloc(qlen + 1, 8);

	init_data(qry, n_cigar_, cigar_, n_col, qlen, tlen, m);

	
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		const int8_t *p_diff = &mat_diff[k * m];
		for (j = 0; j < qlen; ++j) 
		{
			qry->qp_diff[i] = p_diff[seq_nt4_table[(int)query[j]]];
			qry->qp[i++] = p[seq_nt4_table[(int)query[j]]];
		}
	}


	// fill the first row
	// eh[0].h = 0; eh[0].e = MINUS_INF;
	// for (j = 1; j <= qlen && j <= w; ++j)
	// 	eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	// 我觉得最后还是要看下到底是<还是<=
	for (j = 0; j < band_length; ++j)
	{
		qry->eh[j].h = 0;
		///感觉应该是这个东西
		qry->eh[j].e = qry->eh[j].h - gapoe;
		///eh[j].e = eh[j].h - gapo;
	}
	for (; j <= qlen; ++j) 
	{
		qry->eh[j].h = qry->eh[j].e = MINUS_INF; // everything is -inf outside the band
	}

	/// 0 <= i <= (tlen - 1)
	/// 0 <= j <= qlen
	///所以i是从矩阵的第1列而不是第0列开始
	///但是j是从第0行开始的
	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		int32_t f = MINUS_INF, h1, t;
		int8_t *q = &qry->qp[seq_nt4_table[(int)target[i]] * qlen];
		int8_t *q_diff = &qry->qp_diff[seq_nt4_table[(int)target[i]] * qlen];


		Phred = qulity[i] - quality_base;
		if (Phred > 40)
		{
			Phred = 40;
		}
		Phred = Phred/40;
		


		beg = i;
		end = i + band_length; // only loop through [beg,end) of the query sequence
		///h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		h1 = MINUS_INF;

		if (n_cigar_ && cigar_) {
			uint8_t *zi = &qry->z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &qry->eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				///m += q[j];
				///m = m + q[j] - q_diff[j];
				m = m + q[j] - (int)(q_diff[j] * Phred);


				///match/mismatch: d = 0
				///e: d = 1
				///f: d = 2
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				///这个是区分e是gap open 还是gap extension
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;

				t = m - gapoe;
				f -= gape;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &qry->eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				///m += q[j];
				///m = m + q[j] - q_diff[j];
				m = m + q[j] - (int)(q_diff[j] * Phred);



				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				e  = e > t? e : t;
				p->e = e;
				t = m - gapoe;
				f -= gape;
				f  = f > t? f : t;
			}
		}
		qry->eh[end].h = h1; qry->eh[end].e = MINUS_INF;
	}

	///score = eh[qlen].h;
	int max_i = tlen + w;
	score = qry->eh[max_i].h;
	for (i = end; i > beg; i--)
	{
		if (qry->eh[i].h > score)
		{
			score = qry->eh[i].h;
			max_i = i;
		}
	}

	qry->score = score;
	qry->qe = max_i - 1;

	if (n_cigar_ && cigar_) { // backtrack
		///int n_cigar = 0, m_cigar = 0, which = 0;
		///uint32_t *cigar = 0, tmp;
		int n_cigar = 0, which = 0;
		uint32_t tmp;

		// (i,k) points to the last cell
		i = tlen - 1; 
		//k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1;
		k = max_i - 1;
		///k = max_i;

		while (i >= 0 && k >= 0) {
			///which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			which = qry->z[(long)i * n_col + (k - i)] >> (which<<1) & 3;
			if (which == 0)      qry->cigar = push_cigar(&n_cigar, &qry->cigar_size, qry->cigar, 0, 1), --i, --k;
			else if (which == 1) qry->cigar = push_cigar(&n_cigar, &qry->cigar_size, qry->cigar, 2, 1), --i;
			else                 qry->cigar = push_cigar(&n_cigar, &qry->cigar_size, qry->cigar, 1, 1), --k;
		}
		if (i >= 0) qry->cigar = push_cigar(&n_cigar, &qry->cigar_size, qry->cigar, 2, i + 1);
		///if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = qry->cigar[i], qry->cigar[i] = qry->cigar[n_cigar-1-i], qry->cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = qry->cigar;

		qry->qb = k + 1;		
	}

	///free(eh); 
	///free(qp); 
	///free(z);
	return score;
}




///qlen = tlen + 2*w
int ksw_semi_global_path(int qlen, const char *query, int tlen, const char *target, int m, const int8_t *mat, 
int gapo, int gape, int w, kswr_t* qry, int *n_cigar_, uint32_t **cigar_)
{
	///这个就是一列的值
	eh_t *eh;
	///列，也就是query的profile
	int8_t *qp; // query profile
	int i, j, k = MINUS_INF, gapoe = gapo + gape, score;
	int band_length = 2*w+1;
	int32_t beg, end;


	///这些都是求cigar用的
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	int n_col = band_length; // maximum #columns of the backtrack matrix
	z = n_cigar_ && cigar_? (uint8_t *)malloc((long)n_col * tlen) : 0;

	///m是字符集大小
	qp = (int8_t *)malloc(qlen * m);
	///第一列的值
	eh = (eh_t *)calloc(qlen + 1, 8);

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		//for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
		for (j = 0; j < qlen; ++j) qp[i++] = p[seq_nt4_table[(int)query[j]]];
	}


	// fill the first row
	// eh[0].h = 0; eh[0].e = MINUS_INF;
	// for (j = 1; j <= qlen && j <= w; ++j)
	// 	eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	// 我觉得最后还是要看下到底是<还是<=
	for (j = 0; j < band_length; ++j)
	{
		eh[j].h = 0;
		///感觉应该是这个东西
		eh[j].e = eh[j].h - gapoe;
		///eh[j].e = eh[j].h - gapo;
	}
	for (; j <= qlen; ++j) 
	{
		eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	}

	/// 0 <= i <= (tlen - 1)
	/// 0 <= j <= qlen
	///所以i是从矩阵的第1列而不是第0列开始
	///但是j是从第0行开始的
	for (i = 0; LIKELY(i < tlen); ++i)  // target sequence is in the outer loop
	{
		int32_t f = MINUS_INF, h1, t;
		int8_t *q = &qp[seq_nt4_table[(int)target[i]] * qlen];
		beg = i;
		end = i + band_length; // only loop through [beg,end) of the query sequence
		///h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		h1 = MINUS_INF;

		fprintf(stderr, "i: %d\n", i);

		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				m += q[j];
				///match/mismatch: d = 0
				///e: d = 1
				///f: d = 2
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				///这个是区分e是gap open 还是gap extension
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;

				t = m - gapoe;
				f -= gape;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell

				fprintf(stderr, "j: %d, d: %d\n", j, d&3);
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				m += q[j];
				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - gapoe;
				e -= gape;
				e  = e > t? e : t;
				p->e = e;
				t = m - gapoe;
				f -= gape;
				f  = f > t? f : t;
			}
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}

	///score = eh[qlen].h;
	int max_i = tlen + w;
	score = eh[max_i].h;
	for (i = end; i > beg; i--)
	{
		if (eh[i].h > score)
		{
			score = eh[i].h;
			max_i = i;
		}
	}




	


	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;

		// (i,k) points to the last cell
		i = tlen - 1; 
		//k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1;
		k = max_i - 1;
		///k = max_i;


		while (i >= 0 && k >= 0) {
			///which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			which = z[(long)i * n_col + (k - i)] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;


			fprintf(stderr, "i: %d, k: %d\n", i, k);

			
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		///if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;





		/**
		while (i >= 0 && k >= 0) {
			///which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			///which = [0, 1, 2]
			///which<<1 = [0, 2, 4]
			///z里面的每一个元素, 第0位和第1位是h, 都2位和第3位e, 第4位和第6位是f
			which = z[(long)i * n_col + (k - beg)] >> (which<<1) & 3;
			///说明是匹配
			if (which == 0)
			{
				cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1);
				--i;
				--k;
				--beg;
			}      
			else if (which == 1) ///向左
			{
				cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1);
				--i;
				--beg;
			}
			else   ///向上
			{
				cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1);
				--k;
			}                 
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		///if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;
		**/
	}



	qry->score = score;
	qry->qe = max_i - 1;
	qry->qb = k + 1;

	free(eh); 
	free(qp); 
	free(z);
	return score;
}












///p_length > t_length; 
///p_length = t_length + 2 * error_threshold;
inline int try_cigar_without_path(char *pattern, int p_length,
	char *text, int t_length, int end_site, int err, int* return_start_site, int* score,
	int MistMatchPenaltyMax, int MistMatchPenaltyMin, int N_Penalty, char* qulity, int quality_base)
{
	int i = 0;
	int tmp_err = 0;
	int start_site = end_site - t_length + 1;
	(*score) = 0;


	if (start_site >= 0)
	{

		for (i = 0; i < t_length; i++)
		{


			if (text[i] != pattern[i + start_site])
			{

				if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
				{
					tmp_err++;

					if (tmp_err > err)
					{
						return 0;
					}

					if(text[i] == 'N' || pattern[i + start_site] == 'N')
					{
						(*score) -= N_Penalty;
					}
					else
					{
						(*score) -= MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, quality_base, qulity[i]);
					}
			
				}



			}
		}


		if (tmp_err == err)
		{
			(*return_start_site) = start_site;
			return 6;
		}
	}


	return 0;
}



///p_length > t_length; 
///p_length = t_length + 2 * error_threshold;
///end_site and error are old values
///return_start_site, return_end_site, return_err and 
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
	)
{

	int start_site = end_site - t_length + 1;
	int i = 0;
	int tmp_err = 0;

	if ((*return_err) == 0)
	{
		(*score) = 0;
		(*return_start_site) = start_site;
		(*return_end_site) = end_site;
		(*return_err) = error;
		sprintf(cigar, "%dM", t_length);

		return 6;
	}

	if(need_r_quality)
	{
		int end_i = t_length / 2;
		char tmp_c;
		for (i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[t_length - i - 1];
			qulity[t_length - i - 1] = tmp_c;
		}
	}

	if(try_cigar_without_path(pattern, p_length, text, t_length, end_site, error, return_start_site, score, MistMatchPenaltyMax,
	MistMatchPenaltyMin, N_Penalty, qulity, quality_base))
	{
		(*return_end_site) = end_site;
		(*return_err) = error;
		sprintf(cigar, "%dM", t_length);

		if(need_r_quality)
		{
			int end_i = t_length / 2;
			char tmp_c;
			for (i = 0; i < end_i; i++)
			{
				tmp_c = qulity[i];
				qulity[i] = qulity[t_length - i - 1];
				qulity[t_length - i - 1] = tmp_c;
			}
		}

		return 6;
	}

	uint32_t *compress_cigar;
	int n_cigar;
	kswr_t qry;

	///avoid malloc
	///if we just have one candidate location, directly use this method instead of bit-vector
	///and let the size of alphat to be fixed (option m)
	/**
	ksw_semi_global_quality(p_length, pattern, t_length, text, 5, mat, mat_diff, gapo, gape, errthold, qulity, 
	quality_base, qry_fix, &n_cigar, &compress_cigar);
	**/
	///it seems that malloc based method is faster than malloc-free method ...
	///not easy to understand :( 
	ksw_semi_global_quality_back(p_length, pattern, t_length, text, 5, mat, mat_diff, gapo, gape, errthold, qulity, 
	quality_base, &qry, &n_cigar, &compress_cigar);
	//ksw_semi_global(p_length, pattern, t_length, text, 5, mat, gapo, gape, errthold, &qry, &n_cigar, &compress_cigar);
	///to calculate cigar from compress_cigar
	///calculate new edit distance
	qry.NM = 0;
	cigar[0] = '\0';
	int op, opl;
	int insertion_length = 0;

	///int debug_reamove_I = 0;

	/************************************convert 'I' at the beginning of cigar to 'M'****************************************/
	insertion_length = 0;
	///convert 'I' at the begin of cigar to 'M'
	for (i = 0; i < n_cigar; ++i)
	{
		op = compress_cigar[i]&0xf;
		opl = compress_cigar[i]>>4;

		///if this operation is not insertion
		if(op!=2)
		{
			break;
		}
		else
		{
			insertion_length = insertion_length + opl;
		}
	}

	///if i == 0, there is no insertion at the beginning of cigar
	if (i != 0)
	{
		//debug_reamove_I = 1;

		op = compress_cigar[i]&0xf;
		opl = compress_cigar[i]>>4;

		///if next operation is still a 'M', merge
		if (op == 0)
		{
			opl = opl + insertion_length;
		}
		else
		{
			op = 0;
			opl = insertion_length;
			i--;
		}

		compress_cigar[i] = opl;
		compress_cigar[i] = compress_cigar[i] << 4;
		compress_cigar[i] = compress_cigar[i] | op;
		qry.qb = qry.qb - insertion_length;
	}

	int cigar_b = i;
	/************************************convert 'I' at the beginning of cigar to 'M'****************************************/



	/************************************convert 'I' at the end of cigar to 'M'****************************************/
	insertion_length = 0;
	for (i = n_cigar - 1; i >= cigar_b; --i)
	{
		op = compress_cigar[i]&0xf;
		opl = compress_cigar[i]>>4;

		///if this operation is not insertion
		if(op!=2)
		{
			break;
		}
		else
		{
			insertion_length = insertion_length + opl;
		}
	}

	///if i == n_cigar - 1, there is no insertion at the end of cigar
	if (i != n_cigar - 1)
	{
		///debug_reamove_I = 1;

		op = compress_cigar[i]&0xf;
		opl = compress_cigar[i]>>4;

		///if next operation is still a 'M', merge
		if (op == 0)
		{
			opl = opl + insertion_length;
		}
		else
		{
			op = 0;
			opl = insertion_length;
			i++;
		}

		compress_cigar[i] = opl;
		compress_cigar[i] = compress_cigar[i] << 4;
		compress_cigar[i] = compress_cigar[i] | op;
		qry.qe = qry.qe + insertion_length;
	}

	int cigar_e = i;
	/************************************convert 'I' at the end of cigar to 'M'****************************************/
	 


	
	///forward strand
	///if (site < _msf_refGenLength)
	if(is_forward_strand)
	{
		int query_s = qry.qb;
		int target_s = 0;
		int k;

		for (i = cigar_b; i <= cigar_e; ++i) // print CIGAR
		{
			op = compress_cigar[i]&0xf;
			opl = compress_cigar[i]>>4;

			sprintf(cigar + strlen(cigar), "%d%c", opl, "MDISH"[op]);

			if (op == 0)
			{
				for (k = 0; k < opl; k++)
				{
					if((pattern[query_s] != text[target_s]) && 
					!(pattern[query_s] == 'C' && text[target_s] == 'T'))
					{
						qry.NM++;
					}
					query_s++;
					target_s++;
				}
			}
			else if (op == 1)
			{
				query_s += opl;
				qry.NM += opl;
			}
			else
			{
				target_s += opl;
				qry.NM += opl;
			}
		}
		free(compress_cigar);
	}
	else  ///reverse complement strand
	{
		int query_e = qry.qe;
		int target_e = t_length - 1;
		int k;

		///reverse cigar
		for (i = cigar_e; i >= cigar_b; --i) // print CIGAR
		{
			op = compress_cigar[i]&0xf;
			opl = compress_cigar[i]>>4;

			sprintf(cigar + strlen(cigar), "%d%c", opl, "MDISH"[op]);

			if (op == 0)
			{
				for (k = 0; k < opl; k++)
				{
					if((pattern[query_e] != text[target_e]) && 
					!(pattern[query_e] == 'C' && text[target_e] == 'T'))
					{
						qry.NM++;
					}
					query_e--;
					target_e--;
				}
			}
			else if (op == 1)
			{
				query_e -= opl;
				qry.NM += opl;
			}
			else
			{
				target_e -= opl;
				qry.NM += opl;
			}
		}
		free(compress_cigar);
	}
	
	///cigar has been calculated
	(*score) = qry.score;
	(*return_start_site) = qry.qb;
	(*return_end_site) = qry.qe;
	(*return_err) = qry.NM;
	
	if(need_r_quality)
	{
		int end_i = t_length / 2;
		char tmp_c;
		for (i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[t_length - i - 1];
			qulity[t_length - i - 1] = tmp_c;
		}
	}


















	/**
	if(need_r_quality)
	{
		int end_i = t_length / 2;
		char tmp_c;
		for (i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[t_length - i - 1];
			qulity[t_length - i - 1] = tmp_c;
		}
	}

	if (!is_forward_strand)
	{
		char rc_ATGC[256];

		rc_ATGC['A'] = 'T';
		rc_ATGC['C'] = 'G';
		rc_ATGC['G'] = 'C';
		rc_ATGC['T'] = 'A';

		for (i = 0; i < p_length/2; i++)
		{

			char tmp_k = pattern[i];
			pattern[i] = pattern[p_length - i - 1];
			pattern[p_length - i - 1] = tmp_k;
		}

		for (i = 0; i < t_length/2; i++)
		{
			char tmp_k = text[i];
			text[i] = text[t_length - i - 1];
			text[t_length - i - 1] = tmp_k;
		}

		int end_i = t_length / 2;
		char tmp_c;
		for (i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[t_length - i - 1];
			qulity[t_length - i - 1] = tmp_c;
		}
	}
	

	
	

	int debug_cigar_length = strlen(cigar);
	int pre_i = 0;
	char c = 0;
	int p_i, t_i;
	int tmp_error = 0;
	int debug_score;
	int error_sites_text[1000];
	int error_sites_pattern[1000];
	
	debug_score = 0;
	if(is_forward_strand)
	{
		p_i = (*return_start_site);
	}
	else
	{
		p_i = p_length - (*return_end_site) - 1;
	}
	
	
	t_i = 0;
	for (i = 0; i < debug_cigar_length; i++)
	{
		if(cigar[i] > '9' || cigar[i] < '0')
		{
			

			c = cigar[i];

			if(pre_i == 0 && c == 'I')
			{
				fprintf(stderr, "error beginning: cigar: %s\n", cigar);
			}


			cigar[i] = '\0';
			int cLen = atoi(cigar + pre_i);
			cigar[i] = c;
			if(c == 'M')
			{
				int k;
				for (k = 0; k < cLen; k++)
				{
					if (pattern[p_i] != text[t_i]
						&& 
						!(pattern[p_i] == 'C' && text[t_i] == 'T'))
					{
						error_sites_text[tmp_error] = t_i;
						error_sites_pattern[tmp_error] = p_i;

						tmp_error++;

						if(pattern[p_i] == 'N' || text[t_i] == 'N')
						{
							debug_score -= N_Penalty;
						}
						else
						{
							debug_score -= 
								MismatchPenaltyByQuality(MistMatchPenaltyMax, MistMatchPenaltyMin, quality_base, qulity[t_i]);
						}
					}
					
					
					
					p_i++;
					t_i++;
				}
			}
			else if(c == 'I')
			{
				t_i += cLen;
				tmp_error += cLen;

				debug_score = debug_score - gapo - gape*cLen;
			}
			else if(c == 'D')
			{
				p_i += cLen;
				tmp_error += cLen;

				debug_score = debug_score - gapo - gape*cLen;
			}

			pre_i = i + 1;
		}
	}

	if(c == 'I')
	{
		fprintf(stderr, "error end: cigar: %s\n", cigar);
	}

	if (debug_score != (*score) && debug_reamove_I == 0)
	{
		fprintf(stderr, "debug_score: %d, score: %d, is_forward_strand: %d\n", debug_score, (*score), is_forward_strand);
		fprintf(stderr, "cigar: %s, debug_reamove_I: %d\n", cigar, debug_reamove_I);
		fprintf(stderr, "tmp_error: %d\n", tmp_error);
		fprintf(stderr, "p_i: %d\n", p_i);
		fprintf(stderr, "t_i: %d\n", t_i);
		for (i = 0; i < tmp_error; i++)
		{
			fprintf(stderr, "text_i: %d, pattern_i: %d\n", error_sites_text[i], error_sites_pattern[i]);
			fprintf(stderr, "qulity: %d\n", qulity[error_sites_text[i]]);
			
		}

		fprintf(stderr, "pattern:\n");
		for (i = 0; i < p_length; i++)
		{
			fprintf(stderr, "%c", pattern[i]);
		}
		fprintf(stderr, "\n");


		fprintf(stderr, "text:\n");
		for (i = 0; i < t_length; i++)
		{
			fprintf(stderr, "%c", text[i]);
		}
		fprintf(stderr, "\n");
		
	}
	

	if (tmp_error != (*return_err))
	{
		fprintf(stderr, "\nis_forward_strand: %d\n", is_forward_strand);
		fprintf(stderr, "tmp_error: %d\n", tmp_error);
		fprintf(stderr, "(*return_err): %d\n", (*return_err));
		fprintf(stderr, "cigar: %s\n", cigar);

		fprintf(stderr, "(*return_start_site): %d\n",(*return_start_site));
		fprintf(stderr, "(*return_end_site): %d\n",(*return_end_site));

		fprintf(stderr, "pattern: \n");
		for (i = 0; i < p_length; i++)
		{
			fprintf(stderr, "%c", pattern[i]);
		}

		fprintf(stderr, "\n");

		fprintf(stderr, "text: \n");
		for (i = 0; i < t_length; i++)
		{
			fprintf(stderr, "%c", text[i]);
		}

		fprintf(stderr, "\n");

	}


	if (!is_forward_strand)
	{
		char rc_ATGC[256];

		rc_ATGC['A'] = 'T';
		rc_ATGC['C'] = 'G';
		rc_ATGC['G'] = 'C';
		rc_ATGC['T'] = 'A';

		for (i = 0; i < p_length/2; i++)
		{
			char tmp_k = pattern[i];
			pattern[i] = pattern[p_length - i - 1];
			pattern[p_length - i - 1] = tmp_k;
		}

		for (i = 0; i < t_length/2; i++)
		{
			char tmp_k = text[i];
			text[i] = text[t_length - i - 1];
			text[t_length - i - 1] = tmp_k;
		}

		int end_i = t_length / 2;
		char tmp_c;
		for (i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[t_length - i - 1];
			qulity[t_length - i - 1] = tmp_c;
		}
	}

	if(need_r_quality)
	{
		int end_i = t_length / 2;
		char tmp_c;
		for (i = 0; i < end_i; i++)
		{
			tmp_c = qulity[i];
			qulity[i] = qulity[t_length - i - 1];
			qulity[t_length - i - 1] = tmp_c;
		}
	}
	**/
	


	return 1;
}

