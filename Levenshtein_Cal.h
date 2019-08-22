/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/

#ifndef __LEVENSHTEIN_
#define __LEVENSHTEIN_

#include "Auxiliary.h"
#include <iostream>
#include "emmintrin.h"
#include "nmmintrin.h"
#include "smmintrin.h"
#include <immintrin.h>
#include "ksw.h"


///int each_length;
typedef uint64_t Word;
typedef uint16_t Small_Word;
typedef uint32_t Word_32;

///Word D0_arry_64[2000];
///Word HP_arry_64[2000];

///int Route_Size_Whole[2000];
///char Route_Char_Whole[2000];

///char err_match[2000];



inline void dec_bit(Word input)
{
	unsigned int bit_2[64];
	int i = 0;
	for (i = 0; i < 64; i++)
	{

		Word tmp = (Word)((Word)input&(Word)1);

		input = input >> 1;

		if (tmp == (Word)0)
		{
			bit_2[i] = 0;
		}
		else
		{
			bit_2[i] = 1;
		}
	}

	for (i = 63; i >= 0; i--)
	{

		fprintf(stderr, "%u", bit_2[i]);
	}

	fprintf(stderr, "\n");

}

int Calcu_Cigar_MD(int Peq_i, char *pattern, int p_length, char *text, int t_length,
	unsigned short errthold, unsigned short band_down, unsigned short band_below, unsigned short band_length, int err, int match_site, char* cigar, char* MD_Z, int thread_id);

int Brief_2_Banded_BPM_Non_SSE(char *pattern1,char *pattern2,int p_length,char *text,int t_length,
                              int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);



int Brief_Reserve_Banded_BPM_4_high_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);
int Brief_Reserve_Banded_BPM_3_high_SSE(char *pattern1,char *pattern2,char *pattern3,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);

int Start_location_Brief_Reserve_Banded_BPM_4_high_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);
int Start_location_Brief_Reserve_Banded_BPM_3_high_SSE(char *pattern1,char *pattern2,char *pattern3,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);
int Start_location_Reserve_Banded_BPM_high(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err, int thread_id);



int Start_location_Calcu_Cigar_MD(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int pre_err, char* cigar,char* MD_Z,int thread_id,int return_site);


int Start_location_Brief_Reserve_Banded_BPM_8_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,char *pattern8,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);

int Brief_Reserve_Banded_BPM_8_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,char *pattern8,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);

int Brief_Reserve_Banded_BPM_7_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);

int Start_location_Brief_Reserve_Banded_BPM_7_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);

int Brief_Reserve_Banded_BPM_6_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);

int Start_location_Brief_Reserve_Banded_BPM_6_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);


int Brief_Reserve_Banded_BPM_5_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);

int Start_location_Brief_Reserve_Banded_BPM_5_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id);


int Reserve_Banded_BPM(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err, int thread_id);

int Start_location_Reserve_Banded_BPM(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err, int thread_id);




/**************************************************�����ȶ�***************************************************************************************************/

///�������ʵ�ʾͼ���pattern�����п�����N
inline int BS_Reserve_Banded_BPM_back
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, unsigned int* return_err)
{
	(*return_err) = (unsigned int)-1;

	///������Ǹ���ҪԤ����������
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///���ǰ�pattern��ǰ2k + 1���ַ�Ԥ����
	///pattern[0]��ӦPeq[0]
	///pattern[2k]��ӦPeq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	while (i<t_length_1)
	{

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);
		///���б�Խ��߷���ƥ����D0��1
		///�����ƥ����D0��0
		///�����˼����������������Խ����ϵ�б�Խ��߷���������,��ִ���ڲ�����
		///
		if (!(D0&err_mask))
		{
			++err;

			///��ʹȫ���ݼ���Ҳ�ͼ�2k
			if ((err - last_high)>errthold)
				return -1;
		}

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		//����ʵ�����ǰ�pattern[0]�Ƶ���
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///���ǰ��µ�pattern[2k]�ӽ���, ��ò���Ǽӵ�Peq[2k]����
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];
	}


	///fprintf(stderr, "sucess(1)\n");


	///���ѭ���ó�����Ϊ�˷�ֹ�ڴ�й¶
	///��ʵҲ����ѭ��������һ������
	///��ȫ���԰�pattern����һλ
	///��������Ҳ�ã����Լ��ټ��㿪��
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);
	if (!(D0&err_mask))
	{
		++err;
		if ((err - last_high)>errthold)
			return -1;
	}


	////fprintf(stderr, "sucess(2)\n");

	/// last_high = 2k
	/// site = (SEQ_LENGTH + 2k) - 2k -1
	/// site = SEQ_LENGTH - 1
	///��ʱ���siteò���������������Խ��ߵ�λ��
	///int site = p_length - last_high - 1;
	int site = t_length - 1;
	int return_site = -1;
	if ((err <= errthold) && (err<=*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	int i_last = i;
	i = 0;

	while (i<errthold)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err <= *return_err))
		{
			*return_err = err;
			return_site = site + i;
		}
	}


	unsigned int ungap_err;
	ungap_err = err;




	while (i<last_high)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err<=*return_err))
		{
			*return_err = err;
			return_site = site + i;
		}


	}


	if ((ungap_err <= errthold) && (ungap_err == *return_err))
	{
		return_site = site + errthold;
	}


	return return_site;

}










/**************************************************�����ȶ�***************************************************************************************************/

///�������ʵ�ʾͼ���pattern�����п�����N
inline int BS_Reserve_Banded_BPM
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, unsigned int* return_err)
{
	(*return_err) = (unsigned int)-1;

	///������Ǹ���ҪԤ����������
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///���ǰ�pattern��ǰ2k + 1���ַ�Ԥ����
	///pattern[0]��ӦPeq[0]
	///pattern[2k]��ӦPeq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	

	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	while (i<t_length_1)
	{
		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);
		///���б�Խ��߷���ƥ����D0��1
		///�����ƥ����D0��0
		///�����˼����������������Խ����ϵ�б�Խ��߷���������,��ִ���ڲ�����
		///
		if (!(D0&err_mask))
		{
			++err;

			///��ʹȫ���ݼ���Ҳ�ͼ�2k
			if ((err - last_high)>errthold)
				return -1;
		}

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		//����ʵ�����ǰ�pattern[0]�Ƶ���
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///���ǰ��µ�pattern[2k]�ӽ���, ��ò���Ǽӵ�Peq[2k]����
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];
	}




	///fprintf(stderr, "sucess(1)\n");


	///���ѭ���ó�����Ϊ�˷�ֹ�ڴ�й¶
	///��ʵҲ����ѭ��������һ������
	///��ȫ���԰�pattern����һλ
	///��������Ҳ�ã����Լ��ټ��㿪��
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);
	if (!(D0&err_mask))
	{
		++err;
		if ((err - last_high)>errthold)
			return -1;
	}





	////fprintf(stderr, "sucess(2)\n");

	/// last_high = 2k
	/// site = (SEQ_LENGTH + 2k) - 2k -1
	/// site = SEQ_LENGTH - 1
	///��ʱ���siteò���������������Խ��ߵ�λ��
	///int site = p_length - last_high - 1;
	int site = t_length - 1;
	int return_site = -1;
	if ((err <= errthold) && (err<=*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	int i_last = i;
	i = 0;




	while (i<errthold)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err <= *return_err))
		{
			*return_err = err;
			return_site = site + i;
		}
	}


	unsigned int ungap_err;
	ungap_err = err;


	while (i<last_high)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err<=*return_err))
		{
			*return_err = err;
			return_site = site + i;
		}




	}


	if ((ungap_err <= errthold) && (ungap_err == *return_err))
	{
		return_site = site + errthold;
	}

	return return_site;

}








/***************************************************2���ȶ�**************************************************/
inline int BS_Reserve_Banded_BPM_2_SSE_back(char *pattern1, char *pattern2, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold)
{
	return_sites_error[0] = (unsigned int)-1;
	return_sites_error[1] = (unsigned int)-1;
    return_sites[0] = -1;
    return_sites[1] = -1;



	Word Peq[256][2];
	int band_length = (errthold << 1) + 1;

	__m128i Peq_SSE[256];

	int i;

	Word tmp_Peq_1 = 1;




	Peq['A'][0] = 0;
	Peq['A'][1] = 0;

	Peq['C'][0] = 0;
	Peq['C'][1] = 0;

	Peq['G'][0] = 0;
	Peq['G'][1] = 0;

	Peq['T'][0] = 0;
	Peq['T'][1] = 0;



	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'][0] = Peq['T'][0] | Peq['C'][0];
	Peq['T'][1] = Peq['T'][1] | Peq['C'][1];


	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

	Peq_SSE['A'] = _mm_set_epi64x(Peq['A'][1], Peq['A'][0]);
	Peq_SSE['C'] = _mm_set_epi64x(Peq['C'][1], Peq['C'][0]);
	Peq_SSE['G'] = _mm_set_epi64x(Peq['G'][1], Peq['G'][0]);
	Peq_SSE['T'] = _mm_set_epi64x(Peq['T'][1], Peq['T'][0]);




	Word Mask = ((Word)1 << (errthold << 1));

	__m128i Mask1 = _mm_set_epi64x(0, Mask);
	__m128i Mask2 = _mm_set_epi64x(Mask, 0);

	__m128i VP = _mm_setzero_si128();
	__m128i VN = _mm_setzero_si128();
	__m128i X = _mm_setzero_si128();
	__m128i D0 = _mm_setzero_si128();
	__m128i HN = _mm_setzero_si128();
	__m128i HP = _mm_setzero_si128();
	__m128i tmp_process;
    __m128i tmp_process1;

	__m128i Err_2 = _mm_set_epi64x(0, 0);
	__m128i err_mask = _mm_set_epi64x(1, 1);
	///for_not li quan shi 1
	__m128i for_not=_mm_set1_epi32(-1);
	__m128i err_arry;
	__m128i cmp_result;




	/**********************************�����Ϊ��ֹ���****************************/
    Word a_mask=((Word)-1)>>(64-band_length);
    ///zui gao wei shi 0, fang zhi yi chu
    __m128i add_mask=_mm_set_epi64x(a_mask,a_mask);
    /******************************************************************************************/




	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
    int err2;

	__m128i pre_end=_mm_set_epi64x(last_high+errthold,last_high+errthold);

    i=0;
	while (i<t_length_1)
	{


        ///X = Peq[text[i]] | VN;
		X = _mm_or_si128(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
        ///X&VP
		tmp_process1 = _mm_and_si128(X, VP);
		/*****************************************�����������Ǵ����ֹ���*************************/
		tmp_process1 = _mm_and_si128(tmp_process1, add_mask);
		VP = _mm_and_si128(VP, add_mask);
		/*************************************************************************************/
		///(VP + (X&VP))
		tmp_process = _mm_add_epi64(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm_xor_si128(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm_or_si128(tmp_process, X);
        /*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

        ///HN = VP&D0;
		HN = _mm_and_si128(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm_or_si128(D0, VP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		HP = _mm_or_si128(tmp_process, VN);


        ///X = D0 >> 1;
		X = _mm_srli_epi64(D0, 1);
        ///VN = X&HP;
		VN = _mm_and_si128(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm_or_si128(X, HP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		VP = _mm_or_si128(HN, tmp_process);

        ///D0&err_mask
		err_arry = _mm_and_si128(D0, err_mask);
		Err_2 = _mm_add_epi64(Err_2, err_mask);
		Err_2 = _mm_sub_epi64(Err_2, err_arry);

        ///shi ji shang zhe ge zhi hen xiao d
		cmp_result = _mm_cmpgt_epi64(Err_2, pre_end);

        ///jian zhi
		if (_mm_extract_epi64(cmp_result, 0) && _mm_extract_epi64(cmp_result, 1))
			return 1;

		Peq_SSE['A'] = _mm_srli_epi64(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm_srli_epi64(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm_srli_epi64(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm_srli_epi64(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm_or_si128(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm_or_si128(Mask2, Peq_SSE[pattern2[i_bd]]);


		Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);

	}







        ///X = Peq[text[i]] | VN;
    X = _mm_or_si128(Peq_SSE[text[i]], VN);

    /*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
    ///X&VP
    tmp_process1 = _mm_and_si128(X, VP);
    /*****************************************�����������Ǵ����ֹ���*************************/
    tmp_process1 = _mm_and_si128(tmp_process1, add_mask);
    VP = _mm_and_si128(VP, add_mask);
    /*************************************************************************************/
    ///(VP + (X&VP))
    tmp_process = _mm_add_epi64(tmp_process1, VP);
    ///((VP + (X&VP)) ^ VP)
    tmp_process = _mm_xor_si128(tmp_process, VP);
    ///((VP + (X&VP)) ^ VP) | X
    D0 = _mm_or_si128(tmp_process, X);
    /*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

    ///HN = VP&D0;
    HN = _mm_and_si128(D0, VP);

    ///HP = VN | ~(VP | D0);
    tmp_process = _mm_or_si128(D0, VP);
    tmp_process = _mm_andnot_si128(tmp_process, for_not);
    HP = _mm_or_si128(tmp_process, VN);


    ///X = D0 >> 1;
    X = _mm_srli_epi64(D0, 1);
    ///VN = X&HP;
    VN = _mm_and_si128(X, HP);
    ///VP = HN | ~(X | HP);
    tmp_process = _mm_or_si128(X, HP);
    tmp_process = _mm_andnot_si128(tmp_process, for_not);
    VP = _mm_or_si128(HN, tmp_process);

    ///D0&err_mask
    err_arry = _mm_and_si128(D0, err_mask);
    Err_2 = _mm_add_epi64(Err_2, err_mask);
    Err_2 = _mm_sub_epi64(Err_2, err_arry);

    ///shi ji shang zhe ge zhi hen xiao d
    cmp_result = _mm_cmpgt_epi64(Err_2, pre_end);

    ///jian zhi
    if (_mm_extract_epi64(cmp_result, 0) && _mm_extract_epi64(cmp_result, 1))
        return 1;









    int site = t_length - 1;
	err1=_mm_extract_epi64(Err_2,0);
    err2=_mm_extract_epi64(Err_2,1);

    if((err1<=errthold)&&(err1<=return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<=return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }

    i=0;



	while (i<errthold)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm_srli_epi64(VP, i);
		tmp_process = _mm_and_si128(tmp_process, err_mask);
		Err_2 = _mm_add_epi64(Err_2, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm_srli_epi64(VN, i);
		tmp_process1 = _mm_and_si128(tmp_process1, err_mask);
		Err_2 = _mm_sub_epi64(Err_2, tmp_process1);
		++i;

		err1 = _mm_extract_epi64(Err_2, 0);
		err2 = _mm_extract_epi64(Err_2, 1);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
	}



	unsigned int ungap_err1;
	unsigned int ungap_err2;

	ungap_err1 = err1;
	ungap_err2 = err2;



    while(i<last_high)
    {
        ///err = err + ((VP >> i)&(Word)1);
        tmp_process=_mm_srli_epi64(VP,i);
        tmp_process=_mm_and_si128(tmp_process,err_mask);
        Err_2=_mm_add_epi64(Err_2,tmp_process);

        ///err = err - ((VN >> i)&(Word)1);
        tmp_process1=_mm_srli_epi64(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,err_mask);
        Err_2=_mm_sub_epi64(Err_2,tmp_process1);
        ++i;

        err1=_mm_extract_epi64(Err_2,0);
        err2=_mm_extract_epi64(Err_2,1);


        if((err1<=errthold)&&(err1<=return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<=return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
    }



	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}


	return 1;

}








/***************************************************2���ȶ�**************************************************/
inline int BS_Reserve_Banded_BPM_2_SSE(char *pattern1, char *pattern2, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold)
{
	return_sites_error[0] = (unsigned int)-1;
    return_sites_error[1] = (unsigned int)-1;
    return_sites[0] = -1;
    return_sites[1] = -1;



	Word Peq[256][2];
	int band_length = (errthold << 1) + 1;

	__m128i Peq_SSE[256];

	int i;

	Word tmp_Peq_1 = 1;




	Peq['A'][0] = 0;
	Peq['A'][1] = 0;

	Peq['C'][0] = 0;
	Peq['C'][1] = 0;

	Peq['G'][0] = 0;
	Peq['G'][1] = 0;

	Peq['T'][0] = 0;
	Peq['T'][1] = 0;



	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'][0] = Peq['T'][0] | Peq['C'][0];
	Peq['T'][1] = Peq['T'][1] | Peq['C'][1];


	for (i = 0; i < 256; i++)
	{
		Peq_SSE[i] = _mm_setzero_si128();
	}

	Peq_SSE['A'] = _mm_set_epi64x(Peq['A'][1], Peq['A'][0]);
	Peq_SSE['C'] = _mm_set_epi64x(Peq['C'][1], Peq['C'][0]);
	Peq_SSE['G'] = _mm_set_epi64x(Peq['G'][1], Peq['G'][0]);
	Peq_SSE['T'] = _mm_set_epi64x(Peq['T'][1], Peq['T'][0]);




	Word Mask = ((Word)1 << (errthold << 1));

	__m128i Mask1 = _mm_set_epi64x(0, Mask);
	__m128i Mask2 = _mm_set_epi64x(Mask, 0);

	__m128i VP = _mm_setzero_si128();
	__m128i VN = _mm_setzero_si128();
	__m128i X = _mm_setzero_si128();
	__m128i D0 = _mm_setzero_si128();
	__m128i HN = _mm_setzero_si128();
	__m128i HP = _mm_setzero_si128();
	__m128i tmp_process;
    __m128i tmp_process1;

	__m128i Err_2 = _mm_set_epi64x(0, 0);
	__m128i err_mask = _mm_set_epi64x(1, 1);
	///for_not li quan shi 1
	__m128i for_not=_mm_set1_epi32(-1);
	__m128i err_arry;
	__m128i cmp_result;







	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
    int err2;

	__m128i pre_end=_mm_set_epi64x(last_high+errthold,last_high+errthold);

    i=0;
	while (i<t_length_1)
	{


        ///X = Peq[text[i]] | VN;
		X = _mm_or_si128(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
        ///X&VP
		tmp_process1 = _mm_and_si128(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm_add_epi64(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm_xor_si128(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm_or_si128(tmp_process, X);
        /*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

        ///HN = VP&D0;
		HN = _mm_and_si128(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm_or_si128(D0, VP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		HP = _mm_or_si128(tmp_process, VN);


        ///X = D0 >> 1;
		X = _mm_srli_epi64(D0, 1);
        ///VN = X&HP;
		VN = _mm_and_si128(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm_or_si128(X, HP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		VP = _mm_or_si128(HN, tmp_process);

        ///D0&err_mask
		err_arry = _mm_and_si128(D0, err_mask);
		Err_2 = _mm_add_epi64(Err_2, err_mask);
		Err_2 = _mm_sub_epi64(Err_2, err_arry);

        ///shi ji shang zhe ge zhi hen xiao d
		cmp_result = _mm_cmpgt_epi64(Err_2, pre_end);

        ///jian zhi
		if (_mm_extract_epi64(cmp_result, 0) && _mm_extract_epi64(cmp_result, 1))
			return 1;

		Peq_SSE['A'] = _mm_srli_epi64(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm_srli_epi64(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm_srli_epi64(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm_srli_epi64(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm_or_si128(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm_or_si128(Mask2, Peq_SSE[pattern2[i_bd]]);


		Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);

	}







        ///X = Peq[text[i]] | VN;
    X = _mm_or_si128(Peq_SSE[text[i]], VN);

    /*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
    ///X&VP
    tmp_process1 = _mm_and_si128(X, VP);
    ///(VP + (X&VP))
    tmp_process = _mm_add_epi64(tmp_process1, VP);
    ///((VP + (X&VP)) ^ VP)
    tmp_process = _mm_xor_si128(tmp_process, VP);
    ///((VP + (X&VP)) ^ VP) | X
    D0 = _mm_or_si128(tmp_process, X);
    /*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

    ///HN = VP&D0;
    HN = _mm_and_si128(D0, VP);

    ///HP = VN | ~(VP | D0);
    tmp_process = _mm_or_si128(D0, VP);
    tmp_process = _mm_andnot_si128(tmp_process, for_not);
    HP = _mm_or_si128(tmp_process, VN);


    ///X = D0 >> 1;
    X = _mm_srli_epi64(D0, 1);
    ///VN = X&HP;
    VN = _mm_and_si128(X, HP);
    ///VP = HN | ~(X | HP);
    tmp_process = _mm_or_si128(X, HP);
    tmp_process = _mm_andnot_si128(tmp_process, for_not);
    VP = _mm_or_si128(HN, tmp_process);

    ///D0&err_mask
    err_arry = _mm_and_si128(D0, err_mask);
    Err_2 = _mm_add_epi64(Err_2, err_mask);
    Err_2 = _mm_sub_epi64(Err_2, err_arry);

    ///shi ji shang zhe ge zhi hen xiao d
    cmp_result = _mm_cmpgt_epi64(Err_2, pre_end);

    ///jian zhi
    if (_mm_extract_epi64(cmp_result, 0) && _mm_extract_epi64(cmp_result, 1))
        return 1;









    int site = t_length - 1;
	err1=_mm_extract_epi64(Err_2,0);
    err2=_mm_extract_epi64(Err_2,1);

    if((err1<=errthold)&&(err1<=return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<=return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }

    i=0;



	while (i<errthold)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm_srli_epi64(VP, i);
		tmp_process = _mm_and_si128(tmp_process, err_mask);
		Err_2 = _mm_add_epi64(Err_2, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm_srli_epi64(VN, i);
		tmp_process1 = _mm_and_si128(tmp_process1, err_mask);
		Err_2 = _mm_sub_epi64(Err_2, tmp_process1);
		++i;

		err1 = _mm_extract_epi64(Err_2, 0);
		err2 = _mm_extract_epi64(Err_2, 1);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
	}


	unsigned int ungap_err1;
	unsigned int ungap_err2;

	ungap_err1 = err1;
	ungap_err2 = err2;




    while(i<last_high)
    {
        ///err = err + ((VP >> i)&(Word)1);
        tmp_process=_mm_srli_epi64(VP,i);
        tmp_process=_mm_and_si128(tmp_process,err_mask);
        Err_2=_mm_add_epi64(Err_2,tmp_process);

        ///err = err - ((VN >> i)&(Word)1);
        tmp_process1=_mm_srli_epi64(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,err_mask);
        Err_2=_mm_sub_epi64(Err_2,tmp_process1);
        ++i;

        err1=_mm_extract_epi64(Err_2,0);
        err2=_mm_extract_epi64(Err_2,1);


        if((err1<=errthold)&&(err1<=return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<=return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
    }



	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}




	return 1;

}



/*************************************************************************************************************************************************************************/







#if defined __AVX2__


/***************************************************4���ȶ�**************************************************/
inline int BS_Reserve_Banded_BPM_4_SSE_back(char *pattern1, char *pattern2, char *pattern3, char *pattern4, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m256i* Peq_SSE)
{
	return_sites_error[0] = (unsigned int)-1;
	return_sites_error[1] = (unsigned int)-1;
	return_sites_error[2] = (unsigned int)-1;
	return_sites_error[3] = (unsigned int)-1;

	return_sites[0] = -1;
	return_sites[1] = -1;
	return_sites[2] = -1;
	return_sites[3] = -1;



	Word Peq[256][2];
	int band_length = (errthold << 1) + 1;

	int i;

	Word tmp_Peq_1 = 1;




	Peq['A'][0] = 0;
	Peq['A'][1] = 0;
	Peq['A'][2] = 0;
	Peq['A'][3] = 0;

	Peq['C'][0] = 0;
	Peq['C'][1] = 0;
	Peq['C'][2] = 0;
	Peq['C'][3] = 0;

	Peq['G'][0] = 0;
	Peq['G'][1] = 0;
	Peq['G'][2] = 0;
	Peq['G'][3] = 0;

	Peq['T'][0] = 0;
	Peq['T'][1] = 0;
	Peq['T'][2] = 0;
	Peq['T'][3] = 0;




	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;
		Peq[pattern3[i]][2] = Peq[pattern3[i]][2] | tmp_Peq_1;
		Peq[pattern4[i]][3] = Peq[pattern4[i]][3] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'][0] = Peq['T'][0] | Peq['C'][0];
	Peq['T'][1] = Peq['T'][1] | Peq['C'][1];
	Peq['T'][2] = Peq['T'][2] | Peq['C'][2];
	Peq['T'][3] = Peq['T'][3] | Peq['C'][3];



	Peq_SSE['A'] = _mm256_set_epi64x(Peq['A'][3], Peq['A'][2], Peq['A'][1], Peq['A'][0]);
	Peq_SSE['C'] = _mm256_set_epi64x(Peq['C'][3], Peq['C'][2], Peq['C'][1], Peq['C'][0]);
	Peq_SSE['G'] = _mm256_set_epi64x(Peq['G'][3], Peq['G'][2], Peq['G'][1], Peq['G'][0]);
	Peq_SSE['T'] = _mm256_set_epi64x(Peq['T'][3], Peq['T'][2], Peq['T'][1], Peq['T'][0]);



	Word Mask = ((Word)1 << (errthold << 1));

	__m256i Mask1 = _mm256_set_epi64x(0, 0, 0, Mask);
	__m256i Mask2 = _mm256_set_epi64x(0, 0, Mask, 0);
	__m256i Mask3 = _mm256_set_epi64x(0, Mask, 0, 0);
	__m256i Mask4 = _mm256_set_epi64x(Mask, 0, 0, 0);

	__m256i VP = _mm256_setzero_si256();
	__m256i VN = _mm256_setzero_si256();
	__m256i X = _mm256_setzero_si256();
	__m256i D0 = _mm256_setzero_si256();
	__m256i HN = _mm256_setzero_si256();
	__m256i HP = _mm256_setzero_si256();
	__m256i tmp_process;
	__m256i tmp_process1;

	__m256i Err_4 = _mm256_setzero_si256();
	__m256i err_mask = _mm256_set_epi64x(1, 1, 1, 1);
	///for_not li quan shi 1
	__m256i for_not = _mm256_set1_epi32(-1);
	__m256i err_arry;
	__m256i cmp_result;







	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
	int err2;
	int err3;
	int err4;

	__m256i pre_end = _mm256_set_epi64x(last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold);

	i = 0;
	while (i<t_length_1)
	{


		///X = Peq[text[i]] | VN;
		X = _mm256_or_si256(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
		///X&VP
		tmp_process1 = _mm256_and_si256(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm256_add_epi64(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm256_xor_si256(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm256_or_si256(tmp_process, X);
		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

		///HN = VP&D0;
		HN = _mm256_and_si256(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm256_or_si256(D0, VP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		HP = _mm256_or_si256(tmp_process, VN);


		///X = D0 >> 1;
		X = _mm256_srli_epi64(D0, 1);
		///VN = X&HP;
		VN = _mm256_and_si256(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm256_or_si256(X, HP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		VP = _mm256_or_si256(HN, tmp_process);

		///D0&err_mask
		err_arry = _mm256_and_si256(D0, err_mask);
		Err_4 = _mm256_add_epi64(Err_4, err_mask);
		Err_4 = _mm256_sub_epi64(Err_4, err_arry);

		///shi ji shang zhe ge zhi hen xiao d
		cmp_result = _mm256_cmpgt_epi64(Err_4, pre_end);

		///jian zhi
		if (_mm256_extract_epi64(cmp_result, 0) && _mm256_extract_epi64(cmp_result, 1)&
			_mm256_extract_epi64(cmp_result, 2) && _mm256_extract_epi64(cmp_result, 3))
			return 1;


		Peq_SSE['A'] = _mm256_srli_epi64(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm256_srli_epi64(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm256_srli_epi64(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm256_srli_epi64(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm256_or_si256(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm256_or_si256(Mask2, Peq_SSE[pattern2[i_bd]]);
		Peq_SSE[pattern3[i_bd]] = _mm256_or_si256(Mask3, Peq_SSE[pattern3[i_bd]]);
		Peq_SSE[pattern4[i_bd]] = _mm256_or_si256(Mask4, Peq_SSE[pattern4[i_bd]]);


		Peq_SSE['T'] = _mm256_or_si256(Peq_SSE['T'], Peq_SSE['C']);

	}





	///X = Peq[text[i]] | VN;
	X = _mm256_or_si256(Peq_SSE[text[i]], VN);



	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
	///X&VP
	tmp_process1 = _mm256_and_si256(X, VP);
	///(VP + (X&VP))
	tmp_process = _mm256_add_epi64(tmp_process1, VP);
	///((VP + (X&VP)) ^ VP)
	tmp_process = _mm256_xor_si256(tmp_process, VP);
	///((VP + (X&VP)) ^ VP) | X
	D0 = _mm256_or_si256(tmp_process, X);
	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

	///HN = VP&D0;
	HN = _mm256_and_si256(D0, VP);

	///HP = VN | ~(VP | D0);
	tmp_process = _mm256_or_si256(D0, VP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	HP = _mm256_or_si256(tmp_process, VN);


	///X = D0 >> 1;
	X = _mm256_srli_epi64(D0, 1);
	///VN = X&HP;
	VN = _mm256_and_si256(X, HP);
	///VP = HN | ~(X | HP);
	tmp_process = _mm256_or_si256(X, HP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	VP = _mm256_or_si256(HN, tmp_process);

	///D0&err_mask
	err_arry = _mm256_and_si256(D0, err_mask);
	Err_4 = _mm256_add_epi64(Err_4, err_mask);
	Err_4 = _mm256_sub_epi64(Err_4, err_arry);

	///shi ji shang zhe ge zhi hen xiao d
	cmp_result = _mm256_cmpgt_epi64(Err_4, pre_end);

	///jian zhi
	if (_mm256_extract_epi64(cmp_result, 0) && _mm256_extract_epi64(cmp_result, 1)&
		_mm256_extract_epi64(cmp_result, 2) && _mm256_extract_epi64(cmp_result, 3))
		return 1;











	int site = t_length - 1;
	err1 = _mm256_extract_epi64(Err_4, 0);
	err2 = _mm256_extract_epi64(Err_4, 1);
	err3 = _mm256_extract_epi64(Err_4, 2);
	err4 = _mm256_extract_epi64(Err_4, 3);

	if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}
	if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
	{
		return_sites[2] = site;
		return_sites_error[2] = err3;
	}
	if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
	{
		return_sites[3] = site;
		return_sites_error[3] = err4;
	}


	i = 0;



	while (i<errthold)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm256_srli_epi64(VP, i);
		tmp_process = _mm256_and_si256(tmp_process, err_mask);
		Err_4 = _mm256_add_epi64(Err_4, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm256_srli_epi64(VN, i);
		tmp_process1 = _mm256_and_si256(tmp_process1, err_mask);
		Err_4 = _mm256_sub_epi64(Err_4, tmp_process1);

		++i;

		err1 = _mm256_extract_epi64(Err_4, 0);
		err2 = _mm256_extract_epi64(Err_4, 1);
		err3 = _mm256_extract_epi64(Err_4, 2);
		err4 = _mm256_extract_epi64(Err_4, 3);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}




	unsigned int ungap_err1;
	unsigned int ungap_err2;
	unsigned int ungap_err3;
	unsigned int ungap_err4;
	ungap_err1 = err1;
	ungap_err2 = err2;
	ungap_err3 = err3;
	ungap_err4 = err4;




	while (i<last_high)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm256_srli_epi64(VP, i);
		tmp_process = _mm256_and_si256(tmp_process, err_mask);
		Err_4 = _mm256_add_epi64(Err_4, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm256_srli_epi64(VN, i);
		tmp_process1 = _mm256_and_si256(tmp_process1, err_mask);
		Err_4 = _mm256_sub_epi64(Err_4, tmp_process1);

		++i;

		err1 = _mm256_extract_epi64(Err_4, 0);
		err2 = _mm256_extract_epi64(Err_4, 1);
		err3 = _mm256_extract_epi64(Err_4, 2);
		err4 = _mm256_extract_epi64(Err_4, 3);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}




	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}


	if ((ungap_err3 <= errthold) && ungap_err3 == return_sites_error[2])
	{
		return_sites[2] = site + errthold;
	}


	if ((ungap_err4 <= errthold) && ungap_err4 == return_sites_error[3])
	{
		return_sites[3] = site + errthold;
	}


	return 1;

}






/***************************************************4���ȶ�**************************************************/
inline int BS_Reserve_Banded_BPM_4_SSE(char *pattern1, char *pattern2, char *pattern3, char *pattern4, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m256i* Peq_SSE)
{
	/**
	return_sites_error[0] = (unsigned int)-1;
	return_sites_error[1] = (unsigned int)-1;
	return_sites_error[2] = (unsigned int)-1;
	return_sites_error[3] = (unsigned int)-1;

	return_sites[0] = -1;
	return_sites[1] = -1;
	return_sites[2] = -1;
	return_sites[3] = -1;
	**/
	memset(return_sites, -1, sizeof(int)* 4);
	memset(return_sites_error, -1, sizeof(unsigned int)* 4);

	Word Peq[256][4];
	int band_length = (errthold << 1) + 1;

	int i;

	Word tmp_Peq_1 = 1;


	memset(Peq['A'], 0, sizeof(Word)* 4);
	memset(Peq['C'], 0, sizeof(Word)* 4);
	memset(Peq['G'], 0, sizeof(Word)* 4);
	memset(Peq['T'], 0, sizeof(Word)* 4);
	/**
	Peq['A'][0] = 0;
	Peq['A'][1] = 0;
	Peq['A'][2] = 0;
	Peq['A'][3] = 0;

	Peq['C'][0] = 0;
	Peq['C'][1] = 0;
	Peq['C'][2] = 0;
	Peq['C'][3] = 0;

	Peq['G'][0] = 0;
	Peq['G'][1] = 0;
	Peq['G'][2] = 0;
	Peq['G'][3] = 0;

	Peq['T'][0] = 0;
	Peq['T'][1] = 0;
	Peq['T'][2] = 0;
	Peq['T'][3] = 0;
	**/



	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;
		Peq[pattern3[i]][2] = Peq[pattern3[i]][2] | tmp_Peq_1;
		Peq[pattern4[i]][3] = Peq[pattern4[i]][3] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}
	/**
	Peq['T'][0] = Peq['T'][0] | Peq['C'][0];
	Peq['T'][1] = Peq['T'][1] | Peq['C'][1];
	Peq['T'][2] = Peq['T'][2] | Peq['C'][2];
	Peq['T'][3] = Peq['T'][3] | Peq['C'][3];
	**/


	Peq_SSE['A'] = _mm256_set_epi64x(Peq['A'][3], Peq['A'][2], Peq['A'][1], Peq['A'][0]);
	Peq_SSE['C'] = _mm256_set_epi64x(Peq['C'][3], Peq['C'][2], Peq['C'][1], Peq['C'][0]);
	Peq_SSE['G'] = _mm256_set_epi64x(Peq['G'][3], Peq['G'][2], Peq['G'][1], Peq['G'][0]);
	Peq_SSE['T'] = _mm256_set_epi64x(Peq['T'][3], Peq['T'][2], Peq['T'][1], Peq['T'][0]);



	Peq_SSE['T'] = _mm256_or_si256(Peq_SSE['T'], Peq_SSE['C']);

	Word Mask = ((Word)1 << (errthold << 1));

	__m256i Mask1 = _mm256_set_epi64x(0, 0, 0, Mask);
	__m256i Mask2 = _mm256_set_epi64x(0, 0, Mask, 0);
	__m256i Mask3 = _mm256_set_epi64x(0, Mask, 0, 0);
	__m256i Mask4 = _mm256_set_epi64x(Mask, 0, 0, 0);

	__m256i VP = _mm256_setzero_si256();
	__m256i VN = _mm256_setzero_si256();
	__m256i X = _mm256_setzero_si256();
	__m256i D0 = _mm256_setzero_si256();
	__m256i HN = _mm256_setzero_si256();
	__m256i HP = _mm256_setzero_si256();
	__m256i tmp_process;
	__m256i tmp_process1;

	__m256i Err_4 = _mm256_setzero_si256();
	__m256i err_mask = _mm256_set_epi64x(1, 1, 1, 1);
	///for_not li quan shi 1
	__m256i for_not = _mm256_set1_epi32(-1);
	__m256i err_arry;
	__m256i cmp_result;







	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
	int err2;
	int err3;
	int err4;

	__m256i pre_end = _mm256_set_epi64x(last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold);

	i = 0;
	while (i<t_length_1)
	{


		///X = Peq[text[i]] | VN;
		X = _mm256_or_si256(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
		///X&VP
		tmp_process1 = _mm256_and_si256(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm256_add_epi64(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm256_xor_si256(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm256_or_si256(tmp_process, X);
		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

		///HN = VP&D0;
		HN = _mm256_and_si256(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm256_or_si256(D0, VP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		HP = _mm256_or_si256(tmp_process, VN);


		///X = D0 >> 1;
		X = _mm256_srli_epi64(D0, 1);
		///VN = X&HP;
		VN = _mm256_and_si256(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm256_or_si256(X, HP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		VP = _mm256_or_si256(HN, tmp_process);

		///D0&err_mask
		err_arry = _mm256_and_si256(D0, err_mask);
		Err_4 = _mm256_add_epi64(Err_4, err_mask);
		Err_4 = _mm256_sub_epi64(Err_4, err_arry);

		/**
		///shi ji shang zhe ge zhi hen xiao d
		cmp_result = _mm256_cmpgt_epi64(Err_4, pre_end);

		///jian zhi
		if (_mm256_extract_epi64(cmp_result, 0) && _mm256_extract_epi64(cmp_result, 1)&
		_mm256_extract_epi64(cmp_result, 2) && _mm256_extract_epi64(cmp_result, 3))
		return 1;
		**/

		Peq_SSE['A'] = _mm256_srli_epi64(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm256_srli_epi64(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm256_srli_epi64(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm256_srli_epi64(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm256_or_si256(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm256_or_si256(Mask2, Peq_SSE[pattern2[i_bd]]);
		Peq_SSE[pattern3[i_bd]] = _mm256_or_si256(Mask3, Peq_SSE[pattern3[i_bd]]);
		Peq_SSE[pattern4[i_bd]] = _mm256_or_si256(Mask4, Peq_SSE[pattern4[i_bd]]);


		Peq_SSE['T'] = _mm256_or_si256(Peq_SSE['T'], Peq_SSE['C']);

	}





	///X = Peq[text[i]] | VN;
	X = _mm256_or_si256(Peq_SSE[text[i]], VN);



	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
	///X&VP
	tmp_process1 = _mm256_and_si256(X, VP);
	///(VP + (X&VP))
	tmp_process = _mm256_add_epi64(tmp_process1, VP);
	///((VP + (X&VP)) ^ VP)
	tmp_process = _mm256_xor_si256(tmp_process, VP);
	///((VP + (X&VP)) ^ VP) | X
	D0 = _mm256_or_si256(tmp_process, X);
	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

	///HN = VP&D0;
	HN = _mm256_and_si256(D0, VP);

	///HP = VN | ~(VP | D0);
	tmp_process = _mm256_or_si256(D0, VP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	HP = _mm256_or_si256(tmp_process, VN);


	///X = D0 >> 1;
	X = _mm256_srli_epi64(D0, 1);
	///VN = X&HP;
	VN = _mm256_and_si256(X, HP);
	///VP = HN | ~(X | HP);
	tmp_process = _mm256_or_si256(X, HP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	VP = _mm256_or_si256(HN, tmp_process);

	///D0&err_mask
	err_arry = _mm256_and_si256(D0, err_mask);
	Err_4 = _mm256_add_epi64(Err_4, err_mask);
	Err_4 = _mm256_sub_epi64(Err_4, err_arry);

	///shi ji shang zhe ge zhi hen xiao d
	cmp_result = _mm256_cmpgt_epi64(Err_4, pre_end);

	///jian zhi
	if (_mm256_extract_epi64(cmp_result, 0) && _mm256_extract_epi64(cmp_result, 1)&
		_mm256_extract_epi64(cmp_result, 2) && _mm256_extract_epi64(cmp_result, 3))
		return 1;











	int site = t_length - 1;
	err1 = _mm256_extract_epi64(Err_4, 0);
	err2 = _mm256_extract_epi64(Err_4, 1);
	err3 = _mm256_extract_epi64(Err_4, 2);
	err4 = _mm256_extract_epi64(Err_4, 3);

	if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}
	if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
	{
		return_sites[2] = site;
		return_sites_error[2] = err3;
	}
	if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
	{
		return_sites[3] = site;
		return_sites_error[3] = err4;
	}


	i = 0;





	while (i<errthold)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm256_srli_epi64(VP, i);
		tmp_process = _mm256_and_si256(tmp_process, err_mask);
		Err_4 = _mm256_add_epi64(Err_4, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm256_srli_epi64(VN, i);
		tmp_process1 = _mm256_and_si256(tmp_process1, err_mask);
		Err_4 = _mm256_sub_epi64(Err_4, tmp_process1);

		++i;

		err1 = _mm256_extract_epi64(Err_4, 0);
		err2 = _mm256_extract_epi64(Err_4, 1);
		err3 = _mm256_extract_epi64(Err_4, 2);
		err4 = _mm256_extract_epi64(Err_4, 3);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}

	unsigned int ungap_err1;
	unsigned int ungap_err2;
	unsigned int ungap_err3;
	unsigned int ungap_err4;
	ungap_err1 = err1;
	ungap_err2 = err2;
	ungap_err3 = err3;
	ungap_err4 = err4;


	while (i<last_high)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm256_srli_epi64(VP, i);
		tmp_process = _mm256_and_si256(tmp_process, err_mask);
		Err_4 = _mm256_add_epi64(Err_4, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm256_srli_epi64(VN, i);
		tmp_process1 = _mm256_and_si256(tmp_process1, err_mask);
		Err_4 = _mm256_sub_epi64(Err_4, tmp_process1);

		++i;

		err1 = _mm256_extract_epi64(Err_4, 0);
		err2 = _mm256_extract_epi64(Err_4, 1);
		err3 = _mm256_extract_epi64(Err_4, 2);
		err4 = _mm256_extract_epi64(Err_4, 3);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}




	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}


	if ((ungap_err3 <= errthold) && ungap_err3 == return_sites_error[2])
	{
		return_sites[2] = site + errthold;
	}


	if ((ungap_err4 <= errthold) && ungap_err4 == return_sites_error[3])
	{
		return_sites[3] = site + errthold;
	}

	return 1;

}





/***************************************************8���ȶ�**************************************************/
///���û��Ϊungap��·���Ż������Ҫ�����Բο�BS_Reserve_Banded_BPM_4_SSE
inline int BS_Reserve_Banded_BPM_8_SSE(
	char *pattern1, char *pattern2, char *pattern3, char *pattern4,
	char *pattern5, char *pattern6, char *pattern7, char *pattern8,
	int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m256i* Peq_SSE)
{

	memset(return_sites, -1, sizeof(int)* 8);
	memset(return_sites_error, -1, sizeof(unsigned int)* 8);

	Word_32 Peq[256][8];
	int band_length = (errthold << 1) + 1;

	int i;

	Word_32 tmp_Peq_1 = 1;

	memset(Peq['A'], 0, sizeof(Word_32)* 8);
	memset(Peq['C'], 0, sizeof(Word_32)* 8);
	memset(Peq['G'], 0, sizeof(Word_32)* 8);
	memset(Peq['T'], 0, sizeof(Word_32)* 8);


	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;
		Peq[pattern3[i]][2] = Peq[pattern3[i]][2] | tmp_Peq_1;
		Peq[pattern4[i]][3] = Peq[pattern4[i]][3] | tmp_Peq_1;


		Peq[pattern5[i]][4] = Peq[pattern5[i]][4] | tmp_Peq_1;
		Peq[pattern6[i]][5] = Peq[pattern6[i]][5] | tmp_Peq_1;
		Peq[pattern7[i]][6] = Peq[pattern7[i]][6] | tmp_Peq_1;
		Peq[pattern8[i]][7] = Peq[pattern8[i]][7] | tmp_Peq_1;


		tmp_Peq_1 = tmp_Peq_1 << 1;
	}




	Peq_SSE['A'] =
		_mm256_set_epi32(
		Peq['A'][7], Peq['A'][6], Peq['A'][5], Peq['A'][4],
		Peq['A'][3], Peq['A'][2], Peq['A'][1], Peq['A'][0]);

	Peq_SSE['C'] =
		_mm256_set_epi32(
		Peq['C'][7], Peq['C'][6], Peq['C'][5], Peq['C'][4],
		Peq['C'][3], Peq['C'][2], Peq['C'][1], Peq['C'][0]);



	Peq_SSE['G'] =
		_mm256_set_epi32(
		Peq['G'][7], Peq['G'][6], Peq['G'][5], Peq['G'][4],
		Peq['G'][3], Peq['G'][2], Peq['G'][1], Peq['G'][0]);



	Peq_SSE['T'] =
		_mm256_set_epi32(
		Peq['T'][7], Peq['T'][6], Peq['T'][5], Peq['T'][4],
		Peq['T'][3], Peq['T'][2], Peq['T'][1], Peq['T'][0]);

	Peq_SSE['T'] = _mm256_or_si256(Peq_SSE['T'], Peq_SSE['C']);




	Word_32 Mask = ((Word_32)1 << (errthold << 1));

	__m256i Mask1 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, Mask);
	__m256i Mask2 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, Mask, 0);
	__m256i Mask3 = _mm256_set_epi32(0, 0, 0, 0, 0, Mask, 0, 0);
	__m256i Mask4 = _mm256_set_epi32(0, 0, 0, 0, Mask, 0, 0, 0);
	__m256i Mask5 = _mm256_set_epi32(0, 0, 0, Mask, 0, 0, 0, 0);
	__m256i Mask6 = _mm256_set_epi32(0, 0, Mask, 0, 0, 0, 0, 0);
	__m256i Mask7 = _mm256_set_epi32(0, Mask, 0, 0, 0, 0, 0, 0);
	__m256i Mask8 = _mm256_set_epi32(Mask, 0, 0, 0, 0, 0, 0, 0);


	__m256i VP = _mm256_setzero_si256();
	__m256i VN = _mm256_setzero_si256();
	__m256i X = _mm256_setzero_si256();
	__m256i D0 = _mm256_setzero_si256();
	__m256i HN = _mm256_setzero_si256();
	__m256i HP = _mm256_setzero_si256();
	__m256i tmp_process;
	__m256i tmp_process1;

	__m256i Err_32 = _mm256_setzero_si256();
	__m256i err_mask = _mm256_set_epi32(1, 1, 1, 1, 1, 1, 1, 1);
	///for_not li quan shi 1
	__m256i for_not = _mm256_set1_epi32(-1);
	__m256i err_arry;
	__m256i cmp_result;







	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
	int err2;
	int err3;
	int err4;
	int err5;
	int err6;
	int err7;
	int err8;


	__m256i pre_end = _mm256_set_epi32
		(last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold,
		last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold);

	i = 0;
	while (i<t_length_1)
	{


		///X = Peq[text[i]] | VN;
		X = _mm256_or_si256(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
		///X&VP
		tmp_process1 = _mm256_and_si256(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm256_add_epi32(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm256_xor_si256(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm256_or_si256(tmp_process, X);
		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

		///HN = VP&D0;
		HN = _mm256_and_si256(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm256_or_si256(D0, VP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		HP = _mm256_or_si256(tmp_process, VN);


		///X = D0 >> 1;
		X = _mm256_srli_epi32(D0, 1);
		///VN = X&HP;
		VN = _mm256_and_si256(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm256_or_si256(X, HP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		VP = _mm256_or_si256(HN, tmp_process);

		///D0&err_mask
		err_arry = _mm256_and_si256(D0, err_mask);
		Err_32 = _mm256_add_epi32(Err_32, err_mask);
		Err_32 = _mm256_sub_epi32(Err_32, err_arry);



		Peq_SSE['A'] = _mm256_srli_epi32(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm256_srli_epi32(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm256_srli_epi32(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm256_srli_epi32(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm256_or_si256(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm256_or_si256(Mask2, Peq_SSE[pattern2[i_bd]]);
		Peq_SSE[pattern3[i_bd]] = _mm256_or_si256(Mask3, Peq_SSE[pattern3[i_bd]]);
		Peq_SSE[pattern4[i_bd]] = _mm256_or_si256(Mask4, Peq_SSE[pattern4[i_bd]]);
		Peq_SSE[pattern5[i_bd]] = _mm256_or_si256(Mask5, Peq_SSE[pattern5[i_bd]]);
		Peq_SSE[pattern6[i_bd]] = _mm256_or_si256(Mask6, Peq_SSE[pattern6[i_bd]]);
		Peq_SSE[pattern7[i_bd]] = _mm256_or_si256(Mask7, Peq_SSE[pattern7[i_bd]]);
		Peq_SSE[pattern8[i_bd]] = _mm256_or_si256(Mask8, Peq_SSE[pattern8[i_bd]]);


		Peq_SSE['T'] = _mm256_or_si256(Peq_SSE['T'], Peq_SSE['C']);

	}



	///X = Peq[text[i]] | VN;
	X = _mm256_or_si256(Peq_SSE[text[i]], VN);



	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
	///X&VP
	tmp_process1 = _mm256_and_si256(X, VP);
	///(VP + (X&VP))
	tmp_process = _mm256_add_epi32(tmp_process1, VP);
	///((VP + (X&VP)) ^ VP)
	tmp_process = _mm256_xor_si256(tmp_process, VP);
	///((VP + (X&VP)) ^ VP) | X
	D0 = _mm256_or_si256(tmp_process, X);
	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

	///HN = VP&D0;
	HN = _mm256_and_si256(D0, VP);

	///HP = VN | ~(VP | D0);
	tmp_process = _mm256_or_si256(D0, VP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	HP = _mm256_or_si256(tmp_process, VN);


	///X = D0 >> 1;
	X = _mm256_srli_epi32(D0, 1);
	///VN = X&HP;
	VN = _mm256_and_si256(X, HP);
	///VP = HN | ~(X | HP);
	tmp_process = _mm256_or_si256(X, HP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	VP = _mm256_or_si256(HN, tmp_process);

	///D0&err_mask
	err_arry = _mm256_and_si256(D0, err_mask);
	Err_32 = _mm256_add_epi32(Err_32, err_mask);
	Err_32 = _mm256_sub_epi32(Err_32, err_arry);

	///shi ji shang zhe ge zhi hen xiao d
	cmp_result = _mm256_cmpgt_epi32(Err_32, pre_end);

	///jian zhi
	if (_mm256_extract_epi32(cmp_result, 0) && _mm256_extract_epi32(cmp_result, 1)&
		_mm256_extract_epi32(cmp_result, 2) && _mm256_extract_epi32(cmp_result, 3)&
		_mm256_extract_epi32(cmp_result, 4) && _mm256_extract_epi32(cmp_result, 5)&
		_mm256_extract_epi32(cmp_result, 6) && _mm256_extract_epi32(cmp_result, 7))
		return 1;




	int site = t_length - 1;
	err1 = _mm256_extract_epi32(Err_32, 0);
	err2 = _mm256_extract_epi32(Err_32, 1);
	err3 = _mm256_extract_epi32(Err_32, 2);
	err4 = _mm256_extract_epi32(Err_32, 3);

	err5 = _mm256_extract_epi32(Err_32, 4);
	err6 = _mm256_extract_epi32(Err_32, 5);
	err7 = _mm256_extract_epi32(Err_32, 6);
	err8 = _mm256_extract_epi32(Err_32, 7);


	if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}
	if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
	{
		return_sites[2] = site;
		return_sites_error[2] = err3;
	}
	if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
	{
		return_sites[3] = site;
		return_sites_error[3] = err4;
	}
	if ((err5 <= errthold) && (err5 <= return_sites_error[4]))
	{
		return_sites[4] = site;
		return_sites_error[4] = err5;
	}
	if ((err6 <= errthold) && (err6 <= return_sites_error[5]))
	{
		return_sites[5] = site;
		return_sites_error[5] = err6;
	}
	if ((err7 <= errthold) && (err7 <= return_sites_error[6]))
	{
		return_sites[6] = site;
		return_sites_error[6] = err7;
	}
	if ((err8 <= errthold) && (err8 <= return_sites_error[7]))
	{
		return_sites[7] = site;
		return_sites_error[7] = err8;
	}




	i = 0;

	while (i < errthold)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm256_srli_epi32(VP, i);
		tmp_process = _mm256_and_si256(tmp_process, err_mask);
		Err_32 = _mm256_add_epi32(Err_32, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm256_srli_epi32(VN, i);
		tmp_process1 = _mm256_and_si256(tmp_process1, err_mask);
		Err_32 = _mm256_sub_epi32(Err_32, tmp_process1);

		++i;

		err1 = _mm256_extract_epi32(Err_32, 0);
		err2 = _mm256_extract_epi32(Err_32, 1);
		err3 = _mm256_extract_epi32(Err_32, 2);
		err4 = _mm256_extract_epi32(Err_32, 3);

		err5 = _mm256_extract_epi32(Err_32, 4);
		err6 = _mm256_extract_epi32(Err_32, 5);
		err7 = _mm256_extract_epi32(Err_32, 6);
		err8 = _mm256_extract_epi32(Err_32, 7);



		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
		if ((err5 <= errthold) && (err5 <= return_sites_error[4]))
		{
			return_sites[4] = site + i;
			return_sites_error[4] = err5;
		}
		if ((err6 <= errthold) && (err6 <= return_sites_error[5]))
		{
			return_sites[5] = site + i;
			return_sites_error[5] = err6;
		}
		if ((err7 <= errthold) && (err7 <= return_sites_error[6]))
		{
			return_sites[6] = site + i;
			return_sites_error[6] = err7;
		}
		if ((err8 <= errthold) && (err8 <= return_sites_error[7]))
		{
			return_sites[7] = site + i;
			return_sites_error[7] = err8;
		}

	}

	unsigned int ungap_err1;
	unsigned int ungap_err2;
	unsigned int ungap_err3;
	unsigned int ungap_err4;
	unsigned int ungap_err5;
	unsigned int ungap_err6;
	unsigned int ungap_err7;
	unsigned int ungap_err8;
	ungap_err1 = err1;
	ungap_err2 = err2;
	ungap_err3 = err3;
	ungap_err4 = err4;
	ungap_err5 = err5;
	ungap_err6 = err6;
	ungap_err7 = err7;
	ungap_err8 = err8;


	while (i<last_high)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm256_srli_epi32(VP, i);
		tmp_process = _mm256_and_si256(tmp_process, err_mask);
		Err_32 = _mm256_add_epi32(Err_32, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm256_srli_epi32(VN, i);
		tmp_process1 = _mm256_and_si256(tmp_process1, err_mask);
		Err_32 = _mm256_sub_epi32(Err_32, tmp_process1);

		++i;

		err1 = _mm256_extract_epi32(Err_32, 0);
		err2 = _mm256_extract_epi32(Err_32, 1);
		err3 = _mm256_extract_epi32(Err_32, 2);
		err4 = _mm256_extract_epi32(Err_32, 3);

		err5 = _mm256_extract_epi32(Err_32, 4);
		err6 = _mm256_extract_epi32(Err_32, 5);
		err7 = _mm256_extract_epi32(Err_32, 6);
		err8 = _mm256_extract_epi32(Err_32, 7);



		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
		if ((err5 <= errthold) && (err5 <= return_sites_error[4]))
		{
			return_sites[4] = site + i;
			return_sites_error[4] = err5;
		}
		if ((err6 <= errthold) && (err6 <= return_sites_error[5]))
		{
			return_sites[5] = site + i;
			return_sites_error[5] = err6;
		}
		if ((err7 <= errthold) && (err7 <= return_sites_error[6]))
		{
			return_sites[6] = site + i;
			return_sites_error[6] = err7;
		}
		if ((err8 <= errthold) && (err8 <= return_sites_error[7]))
		{
			return_sites[7] = site + i;
			return_sites_error[7] = err8;
		}

	}


	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}


	if ((ungap_err3 <= errthold) && ungap_err3 == return_sites_error[2])
	{
		return_sites[2] = site + errthold;
	}


	if ((ungap_err4 <= errthold) && ungap_err4 == return_sites_error[3])
	{
		return_sites[3] = site + errthold;
	}

	if ((ungap_err5 <= errthold) && ungap_err5 == return_sites_error[4])
	{
		return_sites[4] = site + errthold;
	}

	if ((ungap_err6 <= errthold) && ungap_err6 == return_sites_error[5])
	{
		return_sites[5] = site + errthold;
	}

	if ((ungap_err7 <= errthold) && ungap_err7 == return_sites_error[6])
	{
		return_sites[6] = site + errthold;
	}

	if ((ungap_err8 <= errthold) && ungap_err8 == return_sites_error[7])
	{
		return_sites[7] = site + errthold;
	}

	return 1;

}






/***************************************************16���ȶ�**************************************************/
///���û��Ϊungap��·���Ż������Ҫ�����Բο�BS_Reserve_Banded_BPM_4_SSE
inline int BS_Reserve_Banded_BPM_16_SSE(
	char *pattern1, char *pattern2, char *pattern3, char *pattern4,
	char *pattern5, char *pattern6, char *pattern7, char *pattern8,
	char *pattern9, char *pattern10, char *pattern11, char *pattern12,
	char *pattern13, char *pattern14, char *pattern15, char *pattern16,
	int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m256i* Peq_SSE)
{

	memset(return_sites, -1, sizeof(int)* 16);
	memset(return_sites_error, -1, sizeof(unsigned int)* 16);

	Small_Word Peq[256][16];
	int band_length = (errthold << 1) + 1;

	int i;

	Small_Word tmp_Peq_1 = 1;

	memset(Peq['A'], 0, sizeof(Small_Word)* 16);
	memset(Peq['C'], 0, sizeof(Small_Word)* 16);
	memset(Peq['G'], 0, sizeof(Small_Word)* 16);
	memset(Peq['T'], 0, sizeof(Small_Word)* 16);


	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;
		Peq[pattern3[i]][2] = Peq[pattern3[i]][2] | tmp_Peq_1;
		Peq[pattern4[i]][3] = Peq[pattern4[i]][3] | tmp_Peq_1;


		Peq[pattern5[i]][4] = Peq[pattern5[i]][4] | tmp_Peq_1;
		Peq[pattern6[i]][5] = Peq[pattern6[i]][5] | tmp_Peq_1;
		Peq[pattern7[i]][6] = Peq[pattern7[i]][6] | tmp_Peq_1;
		Peq[pattern8[i]][7] = Peq[pattern8[i]][7] | tmp_Peq_1;


		Peq[pattern9[i]][8] = Peq[pattern9[i]][8] | tmp_Peq_1;
		Peq[pattern10[i]][9] = Peq[pattern10[i]][9] | tmp_Peq_1;
		Peq[pattern11[i]][10] = Peq[pattern11[i]][10] | tmp_Peq_1;
		Peq[pattern12[i]][11] = Peq[pattern12[i]][11] | tmp_Peq_1;


		Peq[pattern13[i]][12] = Peq[pattern13[i]][12] | tmp_Peq_1;
		Peq[pattern14[i]][13] = Peq[pattern14[i]][13] | tmp_Peq_1;
		Peq[pattern15[i]][14] = Peq[pattern15[i]][14] | tmp_Peq_1;
		Peq[pattern16[i]][15] = Peq[pattern16[i]][15] | tmp_Peq_1;


		tmp_Peq_1 = tmp_Peq_1 << 1;
	}




	Peq_SSE['A'] =
		_mm256_set_epi16(
		Peq['A'][15], Peq['A'][14], Peq['A'][13], Peq['A'][12],
		Peq['A'][11], Peq['A'][10], Peq['A'][9], Peq['A'][8],
		Peq['A'][7], Peq['A'][6], Peq['A'][5], Peq['A'][4],
		Peq['A'][3], Peq['A'][2], Peq['A'][1], Peq['A'][0]);

	Peq_SSE['C'] =
		_mm256_set_epi16(
		Peq['C'][15], Peq['C'][14], Peq['C'][13], Peq['C'][12],
		Peq['C'][11], Peq['C'][10], Peq['C'][9], Peq['C'][8],
		Peq['C'][7], Peq['C'][6], Peq['C'][5], Peq['C'][4],
		Peq['C'][3], Peq['C'][2], Peq['C'][1], Peq['C'][0]);



	Peq_SSE['G'] =
		_mm256_set_epi16(
		Peq['G'][15], Peq['G'][14], Peq['G'][13], Peq['G'][12],
		Peq['G'][11], Peq['G'][10], Peq['G'][9], Peq['G'][8],
		Peq['G'][7], Peq['G'][6], Peq['G'][5], Peq['G'][4],
		Peq['G'][3], Peq['G'][2], Peq['G'][1], Peq['G'][0]);



	Peq_SSE['T'] =
		_mm256_set_epi16(
		Peq['T'][15], Peq['T'][14], Peq['T'][13], Peq['T'][12],
		Peq['T'][11], Peq['T'][10], Peq['T'][9], Peq['T'][8],
		Peq['T'][7], Peq['T'][6], Peq['T'][5], Peq['T'][4],
		Peq['T'][3], Peq['T'][2], Peq['T'][1], Peq['T'][0]);

	Peq_SSE['T'] = _mm256_or_si256(Peq_SSE['T'], Peq_SSE['C']);


	/**
	if (strcmp(pattern1, "CTGATTCTTCCTAACCGTGAGCATGGAATATTCTTCCATTGGTTTGTGTCCTCTTTTATTTCGTGGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCGC") == 0)
	{
	fprintf(stderr, "\n\n\n\n\n\nPeq_SSE['A']: %llu\n", _mm256_extract_epi16(Peq_SSE['A'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['A'], 0));

	fprintf(stderr, "Peq_SSE['C']: %llu\n", _mm256_extract_epi16(Peq_SSE['C'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['C'], 0));

	fprintf(stderr, "Peq_SSE['G']: %llu\n", _mm256_extract_epi16(Peq_SSE['G'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['G'], 0));

	fprintf(stderr, "Peq_SSE['T']: %llu\n", _mm256_extract_epi16(Peq_SSE['T'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['T'], 0));

	fprintf(stderr, "#####################\n");
	}
	**/


	Small_Word Mask = ((Small_Word)1 << (errthold << 1));

	__m256i Mask1 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mask);
	__m256i Mask2 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mask, 0);
	__m256i Mask3 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mask, 0, 0);
	__m256i Mask4 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mask, 0, 0, 0);
	__m256i Mask5 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mask, 0, 0, 0, 0);
	__m256i Mask6 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Mask, 0, 0, 0, 0, 0);
	__m256i Mask7 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, Mask, 0, 0, 0, 0, 0, 0);
	__m256i Mask8 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, Mask, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask9 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, Mask, 0, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask10 = _mm256_set_epi16(0, 0, 0, 0, 0, 0, Mask, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask11 = _mm256_set_epi16(0, 0, 0, 0, 0, Mask, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask12 = _mm256_set_epi16(0, 0, 0, 0, Mask, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask13 = _mm256_set_epi16(0, 0, 0, Mask, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask14 = _mm256_set_epi16(0, 0, Mask, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask15 = _mm256_set_epi16(0, Mask, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	__m256i Mask16 = _mm256_set_epi16(Mask, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	__m256i VP = _mm256_setzero_si256();
	__m256i VN = _mm256_setzero_si256();
	__m256i X = _mm256_setzero_si256();
	__m256i D0 = _mm256_setzero_si256();
	__m256i HN = _mm256_setzero_si256();
	__m256i HP = _mm256_setzero_si256();
	__m256i tmp_process;
	__m256i tmp_process1;

	__m256i Err_16 = _mm256_setzero_si256();
	__m256i err_mask = _mm256_set_epi16(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	///for_not li quan shi 1
	__m256i for_not = _mm256_set1_epi32(-1);
	__m256i err_arry;
	__m256i cmp_result;







	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
	int err2;
	int err3;
	int err4;
	int err5;
	int err6;
	int err7;
	int err8;
	int err9;
	int err10;
	int err11;
	int err12;
	int err13;
	int err14;
	int err15;
	int err16;

	__m256i pre_end = _mm256_set_epi16
		(last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold,
		last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold,
		last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold,
		last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold);

	i = 0;
	while (i<t_length_1)
	{


		///X = Peq[text[i]] | VN;
		X = _mm256_or_si256(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
		///X&VP
		tmp_process1 = _mm256_and_si256(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm256_add_epi16(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm256_xor_si256(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm256_or_si256(tmp_process, X);
		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

		///HN = VP&D0;
		HN = _mm256_and_si256(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm256_or_si256(D0, VP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		HP = _mm256_or_si256(tmp_process, VN);


		///X = D0 >> 1;
		X = _mm256_srli_epi16(D0, 1);
		///VN = X&HP;
		VN = _mm256_and_si256(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm256_or_si256(X, HP);
		tmp_process = _mm256_andnot_si256(tmp_process, for_not);
		VP = _mm256_or_si256(HN, tmp_process);

		///D0&err_mask
		err_arry = _mm256_and_si256(D0, err_mask);
		Err_16 = _mm256_add_epi16(Err_16, err_mask);
		Err_16 = _mm256_sub_epi16(Err_16, err_arry);

		/**
		///shi ji shang zhe ge zhi hen xiao d
		cmp_result = _mm256_cmpgt_epi64(Err_4, pre_end);

		///jian zhi
		if (_mm256_extract_epi64(cmp_result, 0) && _mm256_extract_epi64(cmp_result, 1)&
		_mm256_extract_epi64(cmp_result, 2) && _mm256_extract_epi64(cmp_result, 3))
		return 1;
		**/

		Peq_SSE['A'] = _mm256_srli_epi16(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm256_srli_epi16(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm256_srli_epi16(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm256_srli_epi16(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm256_or_si256(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm256_or_si256(Mask2, Peq_SSE[pattern2[i_bd]]);
		Peq_SSE[pattern3[i_bd]] = _mm256_or_si256(Mask3, Peq_SSE[pattern3[i_bd]]);
		Peq_SSE[pattern4[i_bd]] = _mm256_or_si256(Mask4, Peq_SSE[pattern4[i_bd]]);
		Peq_SSE[pattern5[i_bd]] = _mm256_or_si256(Mask5, Peq_SSE[pattern5[i_bd]]);
		Peq_SSE[pattern6[i_bd]] = _mm256_or_si256(Mask6, Peq_SSE[pattern6[i_bd]]);
		Peq_SSE[pattern7[i_bd]] = _mm256_or_si256(Mask7, Peq_SSE[pattern7[i_bd]]);
		Peq_SSE[pattern8[i_bd]] = _mm256_or_si256(Mask8, Peq_SSE[pattern8[i_bd]]);
		Peq_SSE[pattern9[i_bd]] = _mm256_or_si256(Mask9, Peq_SSE[pattern9[i_bd]]);
		Peq_SSE[pattern10[i_bd]] = _mm256_or_si256(Mask10, Peq_SSE[pattern10[i_bd]]);
		Peq_SSE[pattern11[i_bd]] = _mm256_or_si256(Mask11, Peq_SSE[pattern11[i_bd]]);
		Peq_SSE[pattern12[i_bd]] = _mm256_or_si256(Mask12, Peq_SSE[pattern12[i_bd]]);
		Peq_SSE[pattern13[i_bd]] = _mm256_or_si256(Mask13, Peq_SSE[pattern13[i_bd]]);
		Peq_SSE[pattern14[i_bd]] = _mm256_or_si256(Mask14, Peq_SSE[pattern14[i_bd]]);
		Peq_SSE[pattern15[i_bd]] = _mm256_or_si256(Mask15, Peq_SSE[pattern15[i_bd]]);
		Peq_SSE[pattern16[i_bd]] = _mm256_or_si256(Mask16, Peq_SSE[pattern16[i_bd]]);

		Peq_SSE['T'] = _mm256_or_si256(Peq_SSE['T'], Peq_SSE['C']);

	}

	/**
	if (strcmp(pattern1, "CTGATTCTTCCTAACCGTGAGCATGGAATATTCTTCCATTGGTTTGTGTCCTCTTTTATTTCGTGGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCGC") == 0)
	{
	fprintf(stderr, "Peq_SSE['A']: %llu\n", _mm256_extract_epi16(Peq_SSE['A'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['A'], 0));

	fprintf(stderr, "Peq_SSE['C']: %llu\n", _mm256_extract_epi16(Peq_SSE['C'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['C'], 0));

	fprintf(stderr, "Peq_SSE['G']: %llu\n", _mm256_extract_epi16(Peq_SSE['G'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['G'], 0));

	fprintf(stderr, "Peq_SSE['T']: %llu\n", _mm256_extract_epi16(Peq_SSE['T'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['T'], 0));

	fprintf(stderr, "D0_SSE: %llu\n", _mm256_extract_epi16(D0,0));
	dec_bit(_mm256_extract_epi16(D0, 0));

	fprintf(stderr, "HN_SSE: %llu\n", _mm256_extract_epi16(HN, 0));
	dec_bit(_mm256_extract_epi16(HN,0));


	fprintf(stderr, "HP_SSE: %llu\n", _mm256_extract_epi16(HP,0));
	dec_bit(_mm256_extract_epi16(HP, 0));

	fprintf(stderr, "VN_SSE: %llu\n", _mm256_extract_epi16(VN, 0));
	dec_bit(_mm256_extract_epi16(VN,0));

	fprintf(stderr, "VP_SSE: %llu\n", _mm256_extract_epi16(VP,0));
	dec_bit(_mm256_extract_epi16(VP,0));


	fprintf(stderr, "#####################\n");
	}
	**/


	///X = Peq[text[i]] | VN;
	X = _mm256_or_si256(Peq_SSE[text[i]], VN);



	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
	///X&VP
	tmp_process1 = _mm256_and_si256(X, VP);
	///(VP + (X&VP))
	tmp_process = _mm256_add_epi16(tmp_process1, VP);
	///((VP + (X&VP)) ^ VP)
	tmp_process = _mm256_xor_si256(tmp_process, VP);
	///((VP + (X&VP)) ^ VP) | X
	D0 = _mm256_or_si256(tmp_process, X);
	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

	///HN = VP&D0;
	HN = _mm256_and_si256(D0, VP);

	///HP = VN | ~(VP | D0);
	tmp_process = _mm256_or_si256(D0, VP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	HP = _mm256_or_si256(tmp_process, VN);


	///X = D0 >> 1;
	X = _mm256_srli_epi16(D0, 1);
	///VN = X&HP;
	VN = _mm256_and_si256(X, HP);
	///VP = HN | ~(X | HP);
	tmp_process = _mm256_or_si256(X, HP);
	tmp_process = _mm256_andnot_si256(tmp_process, for_not);
	VP = _mm256_or_si256(HN, tmp_process);

	///D0&err_mask
	err_arry = _mm256_and_si256(D0, err_mask);
	Err_16 = _mm256_add_epi16(Err_16, err_mask);
	Err_16 = _mm256_sub_epi16(Err_16, err_arry);

	///shi ji shang zhe ge zhi hen xiao d
	cmp_result = _mm256_cmpgt_epi16(Err_16, pre_end);

	///jian zhi
	if (_mm256_extract_epi16(cmp_result, 0) && _mm256_extract_epi16(cmp_result, 1)&
		_mm256_extract_epi16(cmp_result, 2) && _mm256_extract_epi16(cmp_result, 3)&
		_mm256_extract_epi16(cmp_result, 4) && _mm256_extract_epi16(cmp_result, 5)&
		_mm256_extract_epi16(cmp_result, 6) && _mm256_extract_epi16(cmp_result, 7)&
		_mm256_extract_epi16(cmp_result, 8) && _mm256_extract_epi16(cmp_result, 9)&
		_mm256_extract_epi16(cmp_result, 10) && _mm256_extract_epi16(cmp_result, 11)&
		_mm256_extract_epi16(cmp_result, 12) && _mm256_extract_epi16(cmp_result, 13)&
		_mm256_extract_epi16(cmp_result, 14) && _mm256_extract_epi16(cmp_result, 15))
		return 1;





	/**
	if (strcmp(pattern1, "CTGATTCTTCCTAACCGTGAGCATGGAATATTCTTCCATTGGTTTGTGTCCTCTTTTATTTCGTGGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCGC") == 0)
	{
	fprintf(stderr, "Peq_SSE['A']: %llu\n", _mm256_extract_epi16(Peq_SSE['A'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['A'], 0));

	fprintf(stderr, "Peq_SSE['C']: %llu\n", _mm256_extract_epi16(Peq_SSE['C'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['C'], 0));

	fprintf(stderr, "Peq_SSE['G']: %llu\n", _mm256_extract_epi16(Peq_SSE['G'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['G'], 0));

	fprintf(stderr, "Peq_SSE['T']: %llu\n", _mm256_extract_epi16(Peq_SSE['T'], 0));
	dec_bit(_mm256_extract_epi16(Peq_SSE['T'], 0));

	fprintf(stderr, "D0_SSE: %llu\n", _mm256_extract_epi16(D0, 0));
	dec_bit(_mm256_extract_epi16(D0, 0));

	fprintf(stderr, "HN_SSE: %llu\n", _mm256_extract_epi16(HN, 0));
	dec_bit(_mm256_extract_epi16(HN, 0));


	fprintf(stderr, "HP_SSE: %llu\n", _mm256_extract_epi16(HP, 0));
	dec_bit(_mm256_extract_epi16(HP, 0));

	fprintf(stderr, "VN_SSE: %llu\n", _mm256_extract_epi16(VN, 0));
	dec_bit(_mm256_extract_epi16(VN, 0));

	fprintf(stderr, "VP_SSE: %llu\n", _mm256_extract_epi16(VP, 0));
	dec_bit(_mm256_extract_epi16(VP, 0));

	fprintf(stderr, "err: %llu\n", _mm256_extract_epi16(Err_16, 0));

	fprintf(stderr, "#####################\n");
	}
	**/



	int site = t_length - 1;
	err1 = _mm256_extract_epi16(Err_16, 0);
	err2 = _mm256_extract_epi16(Err_16, 1);
	err3 = _mm256_extract_epi16(Err_16, 2);
	err4 = _mm256_extract_epi16(Err_16, 3);

	err5 = _mm256_extract_epi16(Err_16, 4);
	err6 = _mm256_extract_epi16(Err_16, 5);
	err7 = _mm256_extract_epi16(Err_16, 6);
	err8 = _mm256_extract_epi16(Err_16, 7);

	err9 = _mm256_extract_epi16(Err_16, 8);
	err10 = _mm256_extract_epi16(Err_16, 9);
	err11 = _mm256_extract_epi16(Err_16, 10);
	err12 = _mm256_extract_epi16(Err_16, 11);

	err13 = _mm256_extract_epi16(Err_16, 12);
	err14 = _mm256_extract_epi16(Err_16, 13);
	err15 = _mm256_extract_epi16(Err_16, 14);
	err16 = _mm256_extract_epi16(Err_16, 15);

	if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}
	if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
	{
		return_sites[2] = site;
		return_sites_error[2] = err3;
	}
	if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
	{
		return_sites[3] = site;
		return_sites_error[3] = err4;
	}
	if ((err5 <= errthold) && (err5 <= return_sites_error[4]))
	{
		return_sites[4] = site;
		return_sites_error[4] = err5;
	}
	if ((err6 <= errthold) && (err6 <= return_sites_error[5]))
	{
		return_sites[5] = site;
		return_sites_error[5] = err6;
	}
	if ((err7 <= errthold) && (err7 <= return_sites_error[6]))
	{
		return_sites[6] = site;
		return_sites_error[6] = err7;
	}
	if ((err8 <= errthold) && (err8 <= return_sites_error[7]))
	{
		return_sites[7] = site;
		return_sites_error[7] = err8;
	}
	if ((err9 <= errthold) && (err9 <= return_sites_error[8]))
	{
		return_sites[8] = site;
		return_sites_error[8] = err9;
	}
	if ((err10 <= errthold) && (err10 <= return_sites_error[9]))
	{
		return_sites[9] = site;
		return_sites_error[9] = err10;
	}
	if ((err11 <= errthold) && (err11 <= return_sites_error[10]))
	{
		return_sites[10] = site;
		return_sites_error[10] = err11;
	}
	if ((err12 <= errthold) && (err12 <= return_sites_error[11]))
	{
		return_sites[11] = site;
		return_sites_error[11] = err12;
	}
	if ((err13 <= errthold) && (err13 <= return_sites_error[12]))
	{
		return_sites[12] = site;
		return_sites_error[12] = err13;
	}
	if ((err14 <= errthold) && (err14 <= return_sites_error[13]))
	{
		return_sites[13] = site;
		return_sites_error[13] = err14;
	}
	if ((err15 <= errthold) && (err15 <= return_sites_error[14]))
	{
		return_sites[14] = site;
		return_sites_error[14] = err15;
	}
	if ((err16 <= errthold) && (err16 <= return_sites_error[15]))
	{
		return_sites[15] = site;
		return_sites_error[15] = err16;
	}



	i = 0;



	while (i<last_high)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm256_srli_epi16(VP, i);
		tmp_process = _mm256_and_si256(tmp_process, err_mask);
		Err_16 = _mm256_add_epi16(Err_16, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm256_srli_epi16(VN, i);
		tmp_process1 = _mm256_and_si256(tmp_process1, err_mask);
		Err_16 = _mm256_sub_epi16(Err_16, tmp_process1);

		++i;

		err1 = _mm256_extract_epi16(Err_16, 0);
		err2 = _mm256_extract_epi16(Err_16, 1);
		err3 = _mm256_extract_epi16(Err_16, 2);
		err4 = _mm256_extract_epi16(Err_16, 3);

		err5 = _mm256_extract_epi16(Err_16, 4);
		err6 = _mm256_extract_epi16(Err_16, 5);
		err7 = _mm256_extract_epi16(Err_16, 6);
		err8 = _mm256_extract_epi16(Err_16, 7);

		err9 = _mm256_extract_epi16(Err_16, 8);
		err10 = _mm256_extract_epi16(Err_16, 9);
		err11 = _mm256_extract_epi16(Err_16, 10);
		err12 = _mm256_extract_epi16(Err_16, 11);

		err13 = _mm256_extract_epi16(Err_16, 12);
		err14 = _mm256_extract_epi16(Err_16, 13);
		err15 = _mm256_extract_epi16(Err_16, 14);
		err16 = _mm256_extract_epi16(Err_16, 15);

		/**
		if ((err1 <= errthold) && (err1<return_sites_error[0]))
		{
		return_sites[0] = site + i;
		return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2<return_sites_error[1]))
		{
		return_sites[1] = site + i;
		return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3<return_sites_error[2]))
		{
		return_sites[2] = site + i;
		return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4<return_sites_error[3]))
		{
		return_sites[3] = site + i;
		return_sites_error[3] = err4;
		}
		**/



		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
		if ((err5 <= errthold) && (err5 <= return_sites_error[4]))
		{
			return_sites[4] = site + i;
			return_sites_error[4] = err5;
		}
		if ((err6 <= errthold) && (err6 <= return_sites_error[5]))
		{
			return_sites[5] = site + i;
			return_sites_error[5] = err6;
		}
		if ((err7 <= errthold) && (err7 <= return_sites_error[6]))
		{
			return_sites[6] = site + i;
			return_sites_error[6] = err7;
		}
		if ((err8 <= errthold) && (err8 <= return_sites_error[7]))
		{
			return_sites[7] = site + i;
			return_sites_error[7] = err8;
		}
		if ((err9 <= errthold) && (err9 <= return_sites_error[8]))
		{
			return_sites[8] = site + i;
			return_sites_error[8] = err9;
		}
		if ((err10 <= errthold) && (err10 <= return_sites_error[9]))
		{
			return_sites[9] = site + i;
			return_sites_error[9] = err10;
		}
		if ((err11 <= errthold) && (err11 <= return_sites_error[10]))
		{
			return_sites[10] = site + i;
			return_sites_error[10] = err11;
		}
		if ((err12 <= errthold) && (err12 <= return_sites_error[11]))
		{
			return_sites[11] = site + i;
			return_sites_error[11] = err12;
		}
		if ((err13 <= errthold) && (err13 <= return_sites_error[12]))
		{
			return_sites[12] = site + i;
			return_sites_error[12] = err13;
		}
		if ((err14 <= errthold) && (err14 <= return_sites_error[13]))
		{
			return_sites[13] = site + i;
			return_sites_error[13] = err14;
		}
		if ((err15 <= errthold) && (err15 <= return_sites_error[14]))
		{
			return_sites[14] = site + i;
			return_sites_error[14] = err15;
		}
		if ((err16 <= errthold) && (err16 <= return_sites_error[15]))
		{
			return_sites[15] = site + i;
			return_sites_error[15] = err16;
		}

		/**
		if (strcmp(pattern1, "CTGATTCTTCCTAACCGTGAGCATGGAATATTCTTCCATTGGTTTGTGTCCTCTTTTATTTCGTGGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCGC") == 0)
		{
		fprintf(stderr, "i: %llu\n", i);
		fprintf(stderr, "err1: %llu\n", err1);

		fprintf(stderr, "return_sites[0]: %llu\n", return_sites[0]);
		fprintf(stderr, "return_sites_error[0]: %llu\n", return_sites_error[0]);

		fprintf(stderr, "#####################\n");
		}
		**/



	}




	return 1;

}



#endif







inline int bs_Calculate_Cigar(
	char *pattern, int p_length,
	char *text, int t_length,
	unsigned short errthold,
	int* return_err,
	char* cigar,
	uint16_t* matrix,
	int end_site,
	int* return_start_site,
	char* path
	)
{

	int start_site = end_site - t_length + 1;
	int i = 0;

	if ((*return_err) == 0)
	{
		///��ʵ���ﻹҪ��λpath, Ҫ��pathȫ����Ϊ0
		(*return_start_site) = start_site;


		/**
		std::cerr << "start_site  : " << start_site << std::endl;
		std::cerr << "end_site    : " << end_site << std::endl;
		std::cerr << "return_err    : " << *return_err << std::endl;
		**/

		return 1;
	}




	if (start_site>=0)
	{

		int tmp_err = 0;

		for (i = 0; i<t_length; i++)
		{
			path[i] = 0;
			if (text[i] != pattern[i + start_site])
			{
				path[i] = 1;
				tmp_err++;
			}
		}

		if (tmp_err == (*return_err))
		{
			/**
			std::cerr << "start_site  : " << start_site << std::endl;
			std::cerr << "end_site    : " << end_site << std::endl;
			std::cerr << "return_err    : " << *return_err << std::endl;
			std::cerr << "Direction: ";


			for (i = 0; i < t_length; i++)
			{
				std::cerr << (int)path[i];

			}





			std::cerr << std::endl;

			**/





			for (i = 0; i < t_length; i++)
			{
				std::cerr << (int)path[i];

			}
			std::cerr << std::endl;











			(*return_start_site) = start_site;

			return 1;
		}



	}























	///������Ǹ���ҪԤ����������
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	///int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///���ǰ�pattern��ǰ2k + 1���ַ�Ԥ����
	///pattern[0]��ӦPeq[0]
	///pattern[2k]��ӦPeq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);
	memset(matrix, 0, sizeof(uint16_t)*band_length);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	Word inner_i;
	Word column_start;

	while (i<t_length_1)
	{

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		//����ʵ�����ǰ�pattern[0]�Ƶ���
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///���ǰ��µ�pattern[2k]�ӽ���, ��ò���Ǽӵ�Peq[2k]����
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];


		column_start = i*band_length;

		for (inner_i = 0; inner_i < band_length; inner_i++)
		{
			matrix[column_start + inner_i] = matrix[column_start - band_length + inner_i] + (~(D0 >> inner_i)&err_mask);
		}

		/**
		if (strcmp(text, "AAAAAAAAAAATTAGTTGGGTATGGTGGTGCGTGTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGACATCGTATGCCGTCTTCTGCTTGAAAAA") == 0)
		{
			std::cerr << "-i: " << i << std::endl;
			std::cerr << "-D0: " << D0 << std::endl;
			std::cerr << "-VP: " << VP << std::endl;
			std::cerr << "-VN: " << VN << std::endl;
			std::cerr << "-HP: " << HP << std::endl;
			std::cerr << "-HN: " << HN << std::endl;
		}

		**/




	}


	///fprintf(stderr, "sucess(1)\n");


	///���ѭ���ó�����Ϊ�˷�ֹ�ڴ�й¶
	///��ʵҲ����ѭ��������һ������
	///��ȫ���԰�pattern����һλ
	///��������Ҳ�ã����Լ��ټ��㿪��
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);






	///ע��������(i+1)
	column_start = (i + 1)*band_length;

	for (inner_i = 0; inner_i < band_length; inner_i++)
	{
		matrix[column_start + inner_i] = matrix[column_start - band_length + inner_i] + (~(D0 >> inner_i)&err_mask);
	}

























	int back_track_site = band_length - (p_length - end_site);







	/**
	if (strcmp(text, "AAAAAAAAAAATTAGTTGGGTATGGTGGTGCGTGTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGACATCGTATGCCGTCTTCTGCTTGAAAAA") == 0)
	{

		int min_error = 999;
		int min_index = -1;

		for (inner_i = 0; inner_i < band_length; inner_i++)
		{
			if (matrix[column_start + inner_i]< min_error)
			{
				min_error = matrix[column_start + inner_i];
				min_index = inner_i;
			}
		}



		std::cerr << "Matrix: " << std::endl;


		for (i = 0; i < t_length + 1; i++)
		{
			std::cerr << "i: " << i << std::endl;

			column_start = i*band_length;

			for (inner_i = 0; inner_i < band_length; inner_i++)
			{
				std::cerr << matrix[column_start + inner_i] << " ";
			}

			std::cerr << std::endl;
		}


		std::cerr << "min_error  : " << min_error << std::endl;
		std::cerr << "min_index  : " << min_index << std::endl;
		std::cerr << "back_track_site  : " << back_track_site << std::endl;

	}

	**/














	///char cigar[1000];

	Word v_value, h_value, delta_value, min_value;
	Word direction, is_mismatch; ///0 is match, 1 is mismatch, 2 is xiang shang, 3 is shuiping

	i = t_length;
	int path_length = 0;
	start_site = end_site;


	//Word path_site[1000];


    int low_bound = band_length-1;

	while (i>0)
	{

		column_start = i * band_length;

		///h_value = matrix[column_start - band_length + back_track_site + 1];
		delta_value = matrix[column_start - band_length + back_track_site];


        if(back_track_site == 0)
        {
            h_value = matrix[column_start - band_length + back_track_site + 1];

            min_value = delta_value;
			direction = 0;

			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

        }
        else if(back_track_site == low_bound)
        {
            v_value = matrix[column_start + back_track_site - 1];

            min_value = delta_value;
		    direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}

        }
        else
        {
            h_value = matrix[column_start - band_length + back_track_site + 1];
            v_value = matrix[column_start + back_track_site - 1];


            min_value = delta_value;
	        direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}


			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

        }





		/**
		if(strcmp(text, "AAAAAAAAAAATTAGTTGGGTATGGTGGTGCGTGTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGACATCGTATGCCGTCTTCTGCTTGAAAAA") == 0)
		{

		std::cerr << "i: " << i << std::endl;

		std::cerr << "start_site : " << start_site << std::endl;
		std::cerr << "back_track_site: " << back_track_site << std::endl;

		std::cerr << "delta_value : " << delta_value << std::endl;
		std::cerr << "h_value : " << h_value << std::endl;
		std::cerr << "v_value : " << v_value << std::endl;
		std::cerr << "c_value : " << matrix[column_start + back_track_site] << std::endl;
		std::cerr << "min_value : " << min_value << std::endl;

		}

		**/




		if (direction == 0)
		{

			if (delta_value != matrix[column_start + back_track_site])
			{
				direction = 1;
			}

			i--;

			start_site--;

		}
		if (direction == 2)///ru guo xiang shang yi dong, bing bu huan lie
		{
			back_track_site--;
			start_site--;
		}
		else if (direction == 3)///ru guo xiang zuo yi dong
		{
			i--;
			back_track_site++;
		}


		path[path_length++] = direction;



		/**
		if(strcmp(text, "AAAAAAAAAAATTAGTTGGGTATGGTGGTGCGTGTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGACATCGTATGCCGTCTTCTGCTTGAAAAA") == 0)
		{
		std::cerr << "direction: " << direction << std::endl;
		std::cerr << "*************************************" << std::endl;

		}
		**/


	}

	if (direction != 3)
	{
		start_site++;
	}
	(*return_start_site) = start_site;



        ///if(strcmp(text, "AAAAAAAAATTATATATGATATATATTATAGGTTATGTGTTATGTTTTATATGTAATATATAATTTCATATATTTTATATGTTATATATAATTTTTTTTTA") == 0)

	///if (*return_err>20)
	/**
	if (strcmp(text, "AAAAAAAAAAATTAGTTGGGTATGGTGGTGCGTGTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGACATCGTATGCCGTCTTCTGCTTGAAAAA") == 0)
	   {

		std::cerr << "start_site  : " << start_site << std::endl;
		std::cerr << "end_site    : " << end_site << std::endl;
		std::cerr << "return_err    : " << *return_err << std::endl;
		std::cerr << "path_length : " << path_length << std::endl;
		std::cerr << "Direction: ";


		int debug_err = 0;;

		for (i = path_length - 1; i >= 0; i--)
		{
			std::cerr << (int)path[i];

			if ((int)path[i] != 0)
			{
				debug_err++;
			}
		}





		std::cerr << std::endl;


		if (debug_err != *(return_err))
		{
			std::cerr << "ERROR: debug_err = " << debug_err << " return_err = " << *(return_err) << std::endl;
		}
         }


         **/

        /**
        int debug_err = 0;;

		for (i = path_length - 1; i >= 0; i--)
		{
			///std::cerr << (int)path[i];

			if ((int)path[i] != 0)
			{
				debug_err++;
			}
		}





		///std::cerr << std::endl;


		if (debug_err != *(return_err))
		{
			std::cerr << "ERROR: debug_err = " << debug_err << " return_err = " << *(return_err) << std::endl;
		}
         **/







	int debug_err = 0;;

	for (i = path_length - 1; i >= 0; i--)
	{
		std::cerr << (int)path[i];

		if ((int)path[i] != 0)
		{
			debug_err++;
		}
	}





	std::cerr << std::endl;


	if (debug_err != *(return_err))
	{
		std::cerr << "ERROR: debug_err = " << debug_err << " return_err = " << *(return_err) << std::endl;
	}






	return 1;














}









inline int fast_bs_Calculate_Cigar_back_new(
	char *pattern, int p_length,
	char *text, int t_length,
	unsigned short errthold,
	unsigned int* return_err,
	char* cigar,
	int end_site,
	int* return_start_site,
	char* path,
	int* return_path_length,
	Word* matrix_bit
	)
{

	int start_site = end_site - t_length + 1;
	int i = 0;
	int tmp_err = 0;

	if ((*return_err) == 0)
	{
		///��ʵ���ﻹҪ��λpath, Ҫ��pathȫ����Ϊ0
		(*return_start_site) = start_site;



		return 1;
	}




	if (start_site >= 0)
	{



		for (i = 0; i<t_length; i++)
		{
			///path[i] = 0;
			path[t_length - i - 1] = 0;
			if (text[i] != pattern[i + start_site])
			{

				if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
				{
					path[t_length - i - 1] = 1;
					tmp_err++;
				}



			}
		}






		if (tmp_err == (*return_err))
		{


			(*return_path_length) = t_length;

			(*return_start_site) = start_site;

			return 6;
		}






	}



	///������Ǹ���ҪԤ����������
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	///int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///���ǰ�pattern��ǰ2k + 1���ַ�Ԥ����
	///pattern[0]��ӦPeq[0]
	///pattern[2k]��ӦPeq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	///Word inner_i;
	Word column_start;

	while (i<t_length_1)
	{


		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		X = Peq[text[i]] | VN;



		D0 = ((VP + (X&VP)) ^ VP) | X;




		HN = VP&D0;



		HP = VN | ~(VP | D0);


		X = D0 >> 1;


		VN = X&HP;


		VP = HN | ~(X | HP);

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		//����ʵ�����ǰ�pattern[0]�Ƶ���
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///���ǰ��µ�pattern[2k]�ӽ���, ��ò���Ǽӵ�Peq[2k]����
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];

		///ÿ8��Ԫ��һ��
		column_start = i << 3;

		matrix_bit[column_start] = D0;
		matrix_bit[column_start + 1] = VP;
		matrix_bit[column_start + 2] = VN;
		matrix_bit[column_start + 3] = HP;
		matrix_bit[column_start + 4] = HN;



	}



	///���ѭ���ó�����Ϊ�˷�ֹ�ڴ�й¶
	///��ʵҲ����ѭ��������һ������
	///��ȫ���԰�pattern����һλ
	///��������Ҳ�ã����Լ��ټ��㿪��
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);


	///ÿ8��Ԫ��һ��
	column_start = (i + 1) << 3;

	matrix_bit[column_start] = D0;
	matrix_bit[column_start + 1] = VP;
	matrix_bit[column_start + 2] = VN;
	matrix_bit[column_start + 3] = HP;
	matrix_bit[column_start + 4] = HN;





	int back_track_site = band_length - (p_length - end_site);

	Word v_value, h_value, delta_value, min_value, current_value;
	Word direction, is_mismatch; ///0 is match, 1 is mismatch, 2 is xiang shang, 3 is shuiping

	i = t_length;
	int path_length = 0;
	start_site = end_site;

	current_value = *return_err;


	int low_bound = band_length - 1;

	while (i>0)
	{
		if (current_value == 0)
		{
			break;
		}

		///ÿ8��Ԫ��һ��
		column_start = i << 3;


		///���б�Խ��߷���ƥ����D0��1
		///���Ե�D0��1, ��0
		///D0��0, ��1
		delta_value = current_value -
			((~(matrix_bit[column_start] >> back_track_site))&err_mask);


		if (back_track_site == 0)
		{
			///HP
			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&err_mask);
			//HN
			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&err_mask);






			min_value = delta_value;
			direction = 0;

			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}
		else if (back_track_site == low_bound)
		{
			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);

			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}

		}
		else
		{

			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&(Word)1);


			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&(Word)1);


			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);


			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}


			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}


		if (direction == 0)
		{

			if (delta_value != current_value)
			{
				direction = 1;
			}

			i--;

			start_site--;

		}
		if (direction == 2)///ru guo xiang shang yi dong, bing bu huan lie
		{
			back_track_site--;
			start_site--;
		}
		else if (direction == 3)///ru guo xiang zuo yi dong
		{
			i--;
			back_track_site++;
		}


		path[path_length++] = direction;


		current_value = min_value;

	}


	if (i > 0)
	{
		memset(path + path_length, 0, i);
		start_site = start_site - i;
		direction = 0;
		path_length = path_length + i;
	}

	/**
	while (i > 0)
	{
	path[path_length++] = 0;
	i--;
	start_site--;
	direction = 0;

	}
	**/

	if (direction != 3)
	{
		start_site++;
	}
	(*return_start_site) = start_site;
	(*return_path_length) = path_length;




	return 1;
}




inline int try_cigar(char *pattern, int p_length,
	char *text, int t_length, int end_site, char* path,
	unsigned int* return_err, 
	int* return_start_site, 
	int* return_path_length)
{
	int i = 0;
	int tmp_err = 0;
	int start_site = end_site - t_length + 1;

	if (start_site >= 0)
	{

		for (i = 0; i < t_length; i++)
		{
			///path[i] = 0;
			path[t_length - i - 1] = 0;
			if (text[i] != pattern[i + start_site])
			{

				if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
				{
					path[t_length - i - 1] = 1;
					tmp_err++;

					if (tmp_err > (*return_err))
					{
						return 0;
					}
				}



			}
		}






		if (tmp_err == (*return_err))
		{


			(*return_path_length) = t_length;

			(*return_start_site) = start_site;

			return 6;
		}
	}


	return 0;
}


inline int trim_Ns(char *pattern, int p_length,
	char *text, int t_length, int errthold, unsigned int* return_err,
	int* return_start_site, bitmapper_bs_iter* return_end_site, int pre_end_site)
{
	int i = 0;
	int start_Ns = 0;
	int end_Ns = 0;

	for (i = 0; i < t_length; i++)
	{
		if (text[i] != 'N')
		{
			break;
		}
		else
		{
			start_Ns++;
		}

	}

	for (i = t_length - 1; i>=0; i--)
	{
		if (text[i] != 'N')
		{
			break;
		}
		else
		{
			end_Ns++;
		}
	}




	int total_Ns = start_Ns + end_Ns;
	int tmp_err = 0;
	int start_site = errthold;

	t_length = t_length - end_Ns;

	for (i = start_Ns; i < t_length; i++)
	{

		if (text[i] != pattern[i + start_site])
		{

			if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
			{

				tmp_err++;


			}

		}
	}

	/**
	fprintf(stderr, "(*return_err):%d\n", (*return_err));
	fprintf(stderr, "tmp_err:%d\n", tmp_err);
	fprintf(stderr, "total_Ns:%d\n", total_Ns);
	**/


	if (total_Ns + tmp_err <= errthold)
	{
		if ((*return_err) * 4 > total_Ns * 1 + tmp_err*4)
		{

			fprintf(stderr, "err:%d\n", (*return_err));

			(*return_start_site) = errthold;
			(*return_end_site) = errthold + t_length + end_Ns - 1;
			(*return_err) = total_Ns + tmp_err;

			fprintf(stderr, "err:%d\n", (*return_err));

			return 1;
		}
	}


	return 0;

}




inline int trim_Ns_new(char *pattern, int p_length,
	char *text, int t_length, int errthold, unsigned int* return_err,
	int* return_start_site, bitmapper_bs_iter* return_end_site, int pre_end_site,
	int current_score, int mp, int np, char* cigar)
{
	int i = 0;
	int start_Ns = 0;
	int end_Ns = 0;

	for (i = 0; i < t_length; i++)
	{
		if (text[i] != 'N')
		{
			break;
		}
		else
		{
			start_Ns++;
		}

	}

	for (i = t_length - 1; i >= 0; i--)
	{
		if (text[i] != 'N')
		{
			break;
		}
		else
		{
			end_Ns++;
		}
	}




	int total_Ns = start_Ns + end_Ns;
	int tmp_err = 0;
	int start_site = errthold;

	t_length = t_length - end_Ns;

	for (i = start_Ns; i < t_length; i++)
	{

		if (text[i] != pattern[i + start_site])
		{

			if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
			{

				tmp_err++;


			}

		}
	}




	if (total_Ns + tmp_err <= errthold)
	{
		if (current_score > total_Ns * np + tmp_err * mp)
		{


			(*return_start_site) = errthold;
			(*return_end_site) = errthold + t_length + end_Ns - 1;
			(*return_err) = total_Ns + tmp_err;

			sprintf(cigar, "%dM", t_length);

			return 1;
		}
	}


	return 0;

}






inline int trim_Ns_new_all(char *pattern, int p_length,
	char *text, int t_length, int errthold, unsigned int* return_err,
	int* return_start_site, bitmapper_bs_iter* return_end_site, int pre_end_site,
	int current_score, int mp, int np, char* cigar)
{
	int i = 0;



	int total_Ns = 0;
	int tmp_err = 0;
	int start_site = errthold;
	int score = 0;

	for (i = 0; i < t_length; i++)
	{
		if (text[i] == 'N')
		{
			total_Ns++;
			score = score + np;

			
			if (score > current_score)
			{
				return 0;
			}
			
		}
		else if (text[i] != pattern[i + start_site])
		{

			if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
			{
				tmp_err++;
				score = score + mp;

				
				if (score > current_score)
				{
					return 0;
				}
				
			}
		}


	}


	if (total_Ns + tmp_err <= errthold)
	{
		if (current_score >= score)
		{


			(*return_start_site) = errthold;
			(*return_end_site) = errthold + t_length - 1;
			(*return_err) = total_Ns + tmp_err;

			sprintf(cigar, "%dM", t_length);

			return 1;
		}
	}


	return 0;

}






inline void dec_bit_new(Word input)
{
	uint8_t binary[64];

	int i = 0;

	for (i = 0; i < 64; i++)
	{
		binary[i] = input & (Word)1;
		input = input >> 1;
	}

	for (i = 63; i >= 0; i--)
	{
		fprintf(stderr, "%u", binary[i]);
	}
	fprintf(stderr, "\n");

}

inline int verify_edit_distance(char* src, char* dest)
{
#define MAX_STRING_LEN 300 
#define min(x,y)  ( x<y?x:y )


	int i, j;
	int d[MAX_STRING_LEN][MAX_STRING_LEN] = { 0 };

	for (i = 0; i <= (int)strlen(src); ++i) {
		d[i][0] = 0;
	}
	for (j = 0; j <= (int)strlen(dest); ++j) {
		d[0][j] = j;
	}
	for (i = 1; i <= (int)strlen(src); i++){
		for (j = 1; j <= (int)strlen(dest); j++){
			if ((src[i - 1] == dest[j - 1]) 
				|| 
				(dest[j - 1] == 'T' && src[i - 1] == 'C'))
			{
				d[i][j] = d[i - 1][j - 1];
			}
			else{
				int edIns = d[i][j - 1] + 1;
				int edDel = d[i - 1][j] + 1;
				int edRep = d[i - 1][j - 1] + 1;

				d[i][j] = min(min(edIns, edDel), edRep);
			}
		}
	}
	for (int m = 0; m <= strlen(src); m++){
		for (int n = 0; n <= strlen(dest); n++)
			fprintf(stderr, "%d,", d[m][n]);
		fprintf(stderr,"\n");
	}

	return d[strlen(src)][strlen(dest)];


}

inline int fast_bs_Calculate_Cigar(
	char *pattern, int p_length,
	char *text, int t_length,
	unsigned short errthold,
	unsigned int* return_err,
	char* cigar,
	int end_site,
	int* return_start_site,
	char* path,
	int* return_path_length,
	Word* matrix_bit,
	bitmapper_bs_iter* adjust_end_site
	)
{

	int start_site = end_site - t_length + 1;
	int i = 0;
	int tmp_err = 0;

	if ((*return_err) == 0)
	{
		(*return_start_site) = start_site;

		return 1;
	}



	if (try_cigar(pattern, p_length, text, t_length, end_site, path,
		return_err, return_start_site, return_path_length))
	{
		return 6;
	}





	/**
	if (trim_Ns(pattern, p_length, text, t_length, errthold, return_err, return_start_site, adjust_end_site, end_site))
	{
		return 6;
	}
	**/
	


	///������Ǹ���ҪԤ����������
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	///int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;


	///band_length = 2k + 1
	///���ǰ�pattern��ǰ2k + 1���ַ�Ԥ����
	///pattern[0]��ӦPeq[0]
	///pattern[2k]��ӦPeq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;




	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	///Word inner_i;
	Word column_start;

	int current_err = 0;



	while (i<t_length_1)
	{


		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		X = Peq[text[i]] | VN;


		D0 = ((VP + (X&VP)) ^ VP) | X;


		HN = VP&D0;



		HP = VN | ~(VP | D0);


		X = D0 >> 1;


		VN = X&HP;


		VP = HN | ~(X | HP);




		if (!(D0&err_mask))
		{
			++current_err;
		}




		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		//����ʵ�����ǰ�pattern[0]�Ƶ���
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///���ǰ��µ�pattern[2k]�ӽ���, ��ò���Ǽӵ�Peq[2k]����
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];

		///ÿ8��Ԫ��һ��
		column_start = i << 3;

		matrix_bit[column_start] = D0;
		matrix_bit[column_start + 1] = VP;
		matrix_bit[column_start + 2] = VN;
		matrix_bit[column_start + 3] = HP;
		matrix_bit[column_start + 4] = HN;
	}





	///���ѭ���ó�����Ϊ�˷�ֹ�ڴ�й¶
	///��ʵҲ����ѭ��������һ������
	///��ȫ���԰�pattern����һλ
	///��������Ҳ�ã����Լ��ټ��㿪��
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);


	if (!(D0&err_mask))
	{
		++current_err;
	}



	///ÿ8��Ԫ��һ��
	column_start = (i + 1) << 3;

	matrix_bit[column_start] = D0;
	matrix_bit[column_start + 1] = VP;
	matrix_bit[column_start + 2] = VN;
	matrix_bit[column_start + 3] = HP;
	matrix_bit[column_start + 4] = HN;


	///verify_edit_distance(pattern, text);
	


	/**************************************��һ��Ѷ������Բ�Ҫ**********************************************/
	/**
	int site = t_length - 1;
	int return_site = -1;

	

	fprintf(stderr, "read: %s\n", text);
	fprintf(stderr, "ref : %s\n", pattern);
	fprintf(stderr, "errthold: %d\n", errthold);
	fprintf(stderr, "t_length: %d\n", t_length);
	fprintf(stderr, "p_length: %d\n", p_length);
	

	fprintf(stderr, "current_err: %d, site: %d\n", current_err, site);

	if (current_err == *return_err)
	{
		return_site = site;

		if (try_cigar(pattern, p_length, text, t_length, return_site, path,
			return_err, return_start_site, return_path_length))
		{
			*adjust_end_site = return_site;
			return 6;
		}

	}
	int i_last = i;
	i = 0;



	while (i<errthold)
	{
		current_err = current_err + ((VP >> i)&(Word)1);
		current_err = current_err - ((VN >> i)&(Word)1);
		++i;

		fprintf(stderr, "current_err: %d, site: %d\n", current_err, site + i);

		if (current_err == *return_err)
		{
			return_site = site + i;

			if (try_cigar(pattern, p_length, text, t_length, return_site, path,
				return_err, return_start_site, return_path_length))
			{
				*adjust_end_site = return_site;
				return 6;
			}

		}
	}


	while (i<last_high)
	{
		current_err = current_err + ((VP >> i)&(Word)1);
		current_err = current_err - ((VN >> i)&(Word)1);
		++i;

		fprintf(stderr, "current_err: %d, site: %d\n", current_err, site + i);

		if (current_err == *return_err)
		{
			return_site = site + i;

			if (try_cigar(pattern, p_length, text, t_length, return_site, path,
				return_err, return_start_site, return_path_length))
			{
				*adjust_end_site = return_site;
				return 6;
			}

		}
	}
	
	fprintf(stderr, "\n");
	**/
	/**************************************��һ��Ѷ������Բ�Ҫ**********************************************/









































	int back_track_site = band_length - (p_length - end_site);

	Word v_value, h_value, delta_value, min_value, current_value;
	Word direction, is_mismatch; ///0 is match, 1 is mismatch, 2 is xiang shang, 3 is shuiping

	i = t_length;
	int path_length = 0;
	start_site = end_site;

	current_value = *return_err;


	int low_bound = band_length - 1;

	while (i>0)
	{
		if (current_value == 0)
		{
			break;
		}

		///ÿ8��Ԫ��һ��
		column_start = i << 3;


		///���б�Խ��߷���ƥ����D0��1
		///���Ե�D0��1, ��0
		///D0��0, ��1
		delta_value = current_value -
			((~(matrix_bit[column_start] >> back_track_site))&err_mask);


		if (back_track_site == 0)
		{
			///HP
			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&err_mask);
			//HN
			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&err_mask);






			min_value = delta_value;
			direction = 0;

			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}
		else if (back_track_site == low_bound)
		{
			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);

			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}

		}
		else
		{

			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&(Word)1);


			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&(Word)1);


			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);


			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}


			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}


		if (direction == 0)
		{

			if (delta_value != current_value)
			{
				direction = 1;
			}

			i--;

			start_site--;

		}
		if (direction == 2)///ru guo xiang shang yi dong, bing bu huan lie
		{
			back_track_site--;
			start_site--;
		}
		else if (direction == 3)///ru guo xiang zuo yi dong
		{
			i--;
			back_track_site++;
		}


		path[path_length++] = direction;


		current_value = min_value;

	}


	if (i > 0)
	{
		memset(path + path_length, 0, i);
		start_site = start_site - i;
		direction = 0;
		path_length = path_length + i;
	}

	/**
	while (i > 0)
	{
	path[path_length++] = 0;
	i--;
	start_site--;
	direction = 0;

	}
	**/

	if (direction != 3)
	{
		start_site++;
	}
	(*return_start_site) = start_site;
	(*return_path_length) = path_length;




	return 1;
}





inline int fast_bs_Calculate_Cigar_back(
	char *pattern, int p_length,
	char *text, int t_length,
	unsigned short errthold,
	unsigned int* return_err,
	char* cigar,
	int end_site,
	int* return_start_site,
	char* path,
        int* return_path_length,
	Word* matrix_bit
	)
{

	int start_site = end_site - t_length + 1;
	int i = 0;
	int tmp_err = 0;

	if ((*return_err) == 0)
	{
		///��ʵ���ﻹҪ��λpath, Ҫ��pathȫ����Ϊ0
		(*return_start_site) = start_site;



		return 1;
	}




	if (start_site >= 0)
	{



		for (i = 0; i<t_length; i++)
		{
			///path[i] = 0;
			path[t_length - i - 1] = 0;
			if (text[i] != pattern[i + start_site])
			{

				if (!(text[i] == 'T' && pattern[i + start_site] == 'C'))
				{
					path[t_length - i - 1] = 1;
					tmp_err++;
				}



			}
		}






		if (tmp_err == (*return_err))
		{


            (*return_path_length) = t_length;

			(*return_start_site) = start_site;

			return 6;
		}






	}



	///������Ǹ���ҪԤ����������
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	///int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///���ǰ�pattern��ǰ2k + 1���ַ�Ԥ����
	///pattern[0]��ӦPeq[0]
	///pattern[2k]��ӦPeq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	///Word inner_i;
	Word column_start;

	while (i<t_length_1)
	{


		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		X = Peq[text[i]] | VN;



		D0 = ((VP + (X&VP)) ^ VP) | X;




		HN = VP&D0;



		HP = VN | ~(VP | D0);


		X = D0 >> 1;


		VN = X&HP;


		VP = HN | ~(X | HP);

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		//����ʵ�����ǰ�pattern[0]�Ƶ���
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///���ǰ��µ�pattern[2k]�ӽ���, ��ò���Ǽӵ�Peq[2k]����
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];

		///ÿ8��Ԫ��һ��
		column_start = i << 3;

		matrix_bit[column_start] = D0;
		matrix_bit[column_start + 1] = VP;
		matrix_bit[column_start + 2] = VN;
		matrix_bit[column_start + 3] = HP;
		matrix_bit[column_start + 4] = HN;

		/**
		if (strcmp("NCATGGGATTATAGGCGTGAGTTATTGTGTTTGGTTTTATTTTTTTAATGGGAATTTTAGTTTTATAAGAGTTGTAGTTTGATTTAGAGTTTGTGTTAGTTAGTGGAATAGATGTATTTTATATTTGTAGTATTAGGG", text) == 0)
		{
			fprintf(stderr, "i: %d\n", i);
			dec_bit(matrix_bit[column_start + 3]);
			dec_bit(matrix_bit[column_start + 4]);
			dec_bit(HP);
			dec_bit(HN);

		}
		**/




	}



	///���ѭ���ó�����Ϊ�˷�ֹ�ڴ�й¶
	///��ʵҲ����ѭ��������һ������
	///��ȫ���԰�pattern����һλ
	///��������Ҳ�ã����Լ��ټ��㿪��
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);


	///ÿ8��Ԫ��һ��
	column_start = (i + 1) << 3;

	matrix_bit[column_start] = D0;
	matrix_bit[column_start + 1] = VP;
	matrix_bit[column_start + 2] = VN;
	matrix_bit[column_start + 3] = HP;
	matrix_bit[column_start + 4] = HN;





	int back_track_site = band_length - (p_length - end_site);

	Word v_value, h_value, delta_value, min_value, current_value;
	Word direction, is_mismatch; ///0 is match, 1 is mismatch, 2 is xiang shang, 3 is shuiping
	Word pre_direction = 0;

	i = t_length;
	int path_length = 0;
	start_site = end_site;
	

	current_value = *return_err;


	int low_bound = band_length - 1;

	pre_direction = 0;
	direction = 0;

	while (i>0)
	{
		if (current_value == 0)
		{
			break;
		}

		///ÿ8��Ԫ��һ��
		column_start = i << 3;


		///���б�Խ��߷���ƥ����D0��1
		///���Ե�D0��1, ��0
		///D0��0, ��1
		delta_value = current_value -
			((~(matrix_bit[column_start] >> back_track_site))&err_mask);


		if (back_track_site == 0)
		{
			///HP
			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&err_mask);
			//HN
			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&err_mask);






			min_value = delta_value;
			direction = 0;

			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}
		else if (back_track_site == low_bound)
		{
			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);

			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}

		}
		else
		{

			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&(Word)1);


			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&(Word)1);

	
			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);


			min_value = delta_value;
			direction = 0;

			///if (v_value < min_value)
			if ((v_value < min_value)
				||
				(v_value != current_value && pre_direction == 2 && v_value == min_value))
			{
				min_value = v_value;
				direction = 2;
			}
			
			
			

			/**
			if (strcmp("NCATGGGATTATAGGCGTGAGTTATTGTGTTTGGTTTTATTTTTTTAATGGGAATTTTAGTTTTATAAGAGTTGTAGTTTGATTTAGAGTTTGTGTTAGTTAGTGGAATAGATGTATTTTATATTTGTAGTATTAGGG", text) == 0)
			{
				fprintf(stderr, "\npath_length: %d\n", path_length);
				
				fprintf(stderr, "current_value: %d\n", current_value);
				fprintf(stderr, "delta_value: %d\n", delta_value);
				fprintf(stderr, "h_value: %d\n", h_value);
				fprintf(stderr, "min_value: %d\n", min_value);
				fprintf(stderr, "pre_direction: %d\n", pre_direction);
				fprintf(stderr, "back_track_site: %d\n", back_track_site);
				fprintf(stderr, "low_bound: %d\n", low_bound);
				fprintf(stderr, "i: %d\n", i);
				dec_bit(matrix_bit[column_start + 3]);
				dec_bit(matrix_bit[column_start + 4]);
				
			}
			**/


			///if (h_value < min_value)
			if ((h_value < min_value)
				||
				(h_value != current_value && pre_direction == 3 && h_value == min_value))
			{
				min_value = h_value;
				direction = 3;
			}


			

			

		}


		if (direction == 0)
		{

			if (delta_value != current_value)
			{
				direction = 1;
			}

			i--;

			start_site--;

		}
		if (direction == 2)///ru guo xiang shang yi dong, bing bu huan lie
		{
			back_track_site--;
			start_site--;
		}
		else if (direction == 3)///ru guo xiang zuo yi dong
		{
			i--;
			back_track_site++;
		}


		path[path_length++] = direction;


		current_value = min_value;

		pre_direction = direction;

	}

	
	if (i > 0)
	{
		memset(path + path_length, 0, i);
		start_site = start_site - i;
		direction = 0;
		path_length = path_length + i;
	}
	
	/**
	while (i > 0)
	{
		path[path_length++] = 0;
		i--;
		start_site--;
		direction = 0;

	}
	**/

	if (direction != 3)
	{
		start_site++;
	}
	(*return_start_site) = start_site;
    (*return_path_length) = path_length;




	return 1;
}









inline int fast_bs_Calculate_Cigar_score(
	char *pattern, int p_length,
	char *text, int t_length,
	unsigned short errthold,
	int* return_err,
	char* cigar,
	int end_site,
	int* return_start_site,
	char* path,
	int* return_path_length,
	Word* matrix_bit,
	int* local_score
	)
{

	int start_site = end_site - t_length + 1;
	int i = 0;
	int tmp_err = 0;

	int ma = 2;
	int mp = 6;
	int gap_open = 5;
	int gap_extension = 3;

	if ((*return_err) == 0)
	{
		///��ʵ���ﻹҪ��λpath, Ҫ��pathȫ����Ϊ0
		(*return_start_site) = start_site;

		return 1;
	}




	if (start_site >= 0)
	{

		///��һ���ַ������ó�������

		for (i = 0; i<t_length; i++)
		{

			path[t_length - i - 1] = 0;
			if (text[i] != pattern[i + start_site])
			{
				path[t_length - i - 1] = 1;
				tmp_err++;
			}
		}






		if (tmp_err == (*return_err))
		{

			(*return_path_length) = t_length;

			(*return_start_site) = start_site;

			return 1;
		}






	}


	///������Ǹ���ҪԤ����������
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	///int i = 0;
	Word tmp_Peq_1 = (Word)1;

	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///���ǰ�pattern��ǰ2k + 1���ַ�Ԥ����
	///pattern[0]��ӦPeq[0]
	///pattern[2k]��ӦPeq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = ((Word)1 << (errthold << 1));

	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)

	///Word inner_i;
	Word column_start;

	while (i<t_length_1)
	{

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;

		HP = VN | ~(VP | D0);

		X = D0 >> 1;

		VN = X&HP;

		VP = HN | ~(X | HP);

		///pattern[0]��Peq[2k], ��pattern[2k]��Peq[0]
		//����ʵ�����ǰ�pattern[0]�Ƶ���
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///���ǰ��µ�pattern[2k]�ӽ���, ��ò���Ǽӵ�Peq[2k]����
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


		Peq['T'] = Peq['T'] | Peq['C'];

		///ÿ8��Ԫ��һ��
		column_start = i << 3;

		matrix_bit[column_start] = D0;
		matrix_bit[column_start + 1] = VP;
		matrix_bit[column_start + 2] = VN;
		matrix_bit[column_start + 3] = HP;
		matrix_bit[column_start + 4] = HN;

	}



	///���ѭ���ó�����Ϊ�˷�ֹ�ڴ�й¶
	///��ʵҲ����ѭ��������һ������
	///��ȫ���԰�pattern����һλ
	///��������Ҳ�ã����Լ��ټ��㿪��
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);


	///ÿ8��Ԫ��һ��
	column_start = (i + 1) << 3;

	matrix_bit[column_start] = D0;
	matrix_bit[column_start + 1] = VP;
	matrix_bit[column_start + 2] = VN;
	matrix_bit[column_start + 3] = HP;
	matrix_bit[column_start + 4] = HN;


	int back_track_site = band_length - (p_length - end_site);

	Word v_value, h_value, delta_value, min_value, current_value;
	Word direction, is_mismatch; ///0 is match, 1 is mismatch, 2 is xiang shang, 3 is shuiping

	i = t_length;
	int path_length = 0;
	start_site = end_site;

	current_value = *return_err;

	int low_bound = band_length - 1;

	while (i>0)
	{

		///ÿ8��Ԫ��һ��
		column_start = i << 3;


		///���б�Խ��߷���ƥ����D0��1
		///���Ե�D0��1, ��0
		///D0��0, ��1
		delta_value = current_value -
			((~(matrix_bit[column_start] >> back_track_site))&err_mask);


		if (back_track_site == 0)
		{
			///HP
			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&err_mask);
			//HN
			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&err_mask);






			min_value = delta_value;
			direction = 0;

			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}
		else if (back_track_site == low_bound)
		{
			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);

			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}

		}
		else
		{

			h_value = current_value - ((matrix_bit[column_start + 3] >> back_track_site)&(Word)1);

			h_value = h_value + ((matrix_bit[column_start + 4] >> back_track_site)&(Word)1);

			v_value = current_value - ((matrix_bit[column_start + 1] >> (back_track_site - 1))&err_mask);
			v_value = v_value + ((matrix_bit[column_start + 2] >> (back_track_site - 1))&err_mask);


			min_value = delta_value;
			direction = 0;

			if (v_value < min_value)
			{
				min_value = v_value;
				direction = 2;
			}


			if (h_value < min_value)
			{
				min_value = h_value;
				direction = 3;
			}

		}


		if (direction == 0)
		{

			if (delta_value != current_value)
			{
				direction = 1;
			}

			i--;

			start_site--;

		}
		if (direction == 2)///ru guo xiang shang yi dong, bing bu huan lie
		{
			back_track_site--;
			start_site--;
		}
		else if (direction == 3)///ru guo xiang zuo yi dong
		{
			i--;
			back_track_site++;
		}


		path[path_length++] = direction;


		current_value = min_value;

	}

	if (direction != 3)
	{
		start_site++;
	}
	(*return_start_site) = start_site;
	(*return_path_length) = path_length;


	return 1;
}










/***************************************************4���ȶԼ�sse**************************************************/
inline int BS_Reserve_Banded_BPM_4_SSE_only(char *pattern1, char *pattern2, char *pattern3, char *pattern4, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m128i* Peq_SSE)

{


	memset(return_sites, -1, sizeof(int)* 4);
	memset(return_sites_error, -1, sizeof(unsigned int)* 4);

	Word_32 Peq[256][4];
	int band_length = (errthold << 1) + 1;


	int i;

	Word_32 tmp_Peq_1 = 1;


	memset(Peq['A'], 0, sizeof(Word_32)* 4);
	memset(Peq['C'], 0, sizeof(Word_32)* 4);
	memset(Peq['G'], 0, sizeof(Word_32)* 4);
	memset(Peq['T'], 0, sizeof(Word_32)* 4);

	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;
		Peq[pattern3[i]][2] = Peq[pattern3[i]][2] | tmp_Peq_1;
		Peq[pattern4[i]][3] = Peq[pattern4[i]][3] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq_SSE['A'] = _mm_set_epi32(Peq['A'][3], Peq['A'][2], Peq['A'][1], Peq['A'][0]);
	Peq_SSE['C'] = _mm_set_epi32(Peq['C'][3], Peq['C'][2], Peq['C'][1], Peq['C'][0]);
	Peq_SSE['G'] = _mm_set_epi32(Peq['G'][3], Peq['G'][2], Peq['G'][1], Peq['G'][0]);
	Peq_SSE['T'] = _mm_set_epi32(Peq['T'][3], Peq['T'][2], Peq['T'][1], Peq['T'][0]);

	Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);



	Word_32 Mask = ((Word_32)1 << (errthold << 1));

	__m128i Mask1 = _mm_set_epi32(0, 0, 0, Mask);
	__m128i Mask2 = _mm_set_epi32(0, 0, Mask, 0);
	__m128i Mask3 = _mm_set_epi32(0, Mask, 0, 0);
	__m128i Mask4 = _mm_set_epi32(Mask, 0, 0, 0);


	__m128i VP = _mm_setzero_si128();
	__m128i VN = _mm_setzero_si128();
	__m128i X = _mm_setzero_si128();
	__m128i D0 = _mm_setzero_si128();
	__m128i HN = _mm_setzero_si128();
	__m128i HP = _mm_setzero_si128();
	__m128i tmp_process;
	__m128i tmp_process1;



	__m128i Err_4 = _mm_setzero_si128();
	__m128i err_mask = _mm_set_epi32(1, 1, 1, 1);
	///for_not li quan shi 1
	__m128i for_not = _mm_set1_epi32(-1);
	__m128i err_arry;
	__m128i cmp_result;




	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
	int err2;
	int err3;
	int err4;



	__m128i pre_end = _mm_set_epi32
		(last_high + errthold, last_high + errthold, last_high + errthold, last_high + errthold);


	i = 0;
	while (i<t_length_1)
	{


		///X = Peq[text[i]] | VN;
		X = _mm_or_si128(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
		///X&VP
		tmp_process1 = _mm_and_si128(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm_add_epi32(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm_xor_si128(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm_or_si128(tmp_process, X);
		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

		///HN = VP&D0;
		HN = _mm_and_si128(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm_or_si128(D0, VP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		HP = _mm_or_si128(tmp_process, VN);


		///X = D0 >> 1;
		X = _mm_srli_epi32(D0, 1);
		///VN = X&HP;
		VN = _mm_and_si128(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm_or_si128(X, HP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		VP = _mm_or_si128(HN, tmp_process);

		///D0&err_mask
		err_arry = _mm_and_si128(D0, err_mask);
		Err_4 = _mm_add_epi32(Err_4, err_mask);
		Err_4 = _mm_sub_epi32(Err_4, err_arry);



		Peq_SSE['A'] = _mm_srli_epi32(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm_srli_epi32(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm_srli_epi32(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm_srli_epi32(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm_or_si128(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm_or_si128(Mask2, Peq_SSE[pattern2[i_bd]]);
		Peq_SSE[pattern3[i_bd]] = _mm_or_si128(Mask3, Peq_SSE[pattern3[i_bd]]);
		Peq_SSE[pattern4[i_bd]] = _mm_or_si128(Mask4, Peq_SSE[pattern4[i_bd]]);


		Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);

	}







	///X = Peq[text[i]] | VN;
	X = _mm_or_si128(Peq_SSE[text[i]], VN);

	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
	///X&VP
	tmp_process1 = _mm_and_si128(X, VP);
	///(VP + (X&VP))
	tmp_process = _mm_add_epi32(tmp_process1, VP);
	///((VP + (X&VP)) ^ VP)
	tmp_process = _mm_xor_si128(tmp_process, VP);
	///((VP + (X&VP)) ^ VP) | X
	D0 = _mm_or_si128(tmp_process, X);
	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

	///HN = VP&D0;
	HN = _mm_and_si128(D0, VP);

	///HP = VN | ~(VP | D0);
	tmp_process = _mm_or_si128(D0, VP);
	tmp_process = _mm_andnot_si128(tmp_process, for_not);
	HP = _mm_or_si128(tmp_process, VN);


	///X = D0 >> 1;
	X = _mm_srli_epi32(D0, 1);
	///VN = X&HP;
	VN = _mm_and_si128(X, HP);
	///VP = HN | ~(X | HP);
	tmp_process = _mm_or_si128(X, HP);
	tmp_process = _mm_andnot_si128(tmp_process, for_not);
	VP = _mm_or_si128(HN, tmp_process);

	///D0&err_mask
	err_arry = _mm_and_si128(D0, err_mask);
	Err_4 = _mm_add_epi32(Err_4, err_mask);
	Err_4 = _mm_sub_epi32(Err_4, err_arry);

	///shi ji shang zhe ge zhi hen xiao d
	cmp_result = _mm_cmpgt_epi32(Err_4, pre_end);

	///jian zhi
	if (_mm_extract_epi32(cmp_result, 0) && _mm_extract_epi32(cmp_result, 1)
		&& _mm_extract_epi32(cmp_result, 2) && _mm_extract_epi32(cmp_result, 3))
		return 1;








	int site = t_length - 1;
	err1 = _mm_extract_epi32(Err_4, 0);
	err2 = _mm_extract_epi32(Err_4, 1);
	err3 = _mm_extract_epi32(Err_4, 2);
	err4 = _mm_extract_epi32(Err_4, 3);









	/**
	if (strcmp(text, "TTTTAGGTTTAAGTAATTTTTTTGTTTTAGTTTTTTAAGTAGTTGGGGTTATAGGTGTATGTTATTATNTTTAGTTAATTTTTT") == 0
	&&
	strcmp(pattern1, "TAGCTCACTGAAGCCTCAGATGATCCTCCCACCTCAGCCTCCCAAGTAGCTGGGGCTACAGGTGCATGCTACCACGACTGGCTAATTTTAATATTTTTA") == 0)
	{
	fprintf(stderr, "\n\n\n\n*****************************\ntext: %s\n", text);
	fprintf(stderr, "text: %s\n", pattern1);


	fprintf(stderr, "err1: %d\n", err1);
	fprintf(stderr, "err2: %d\n", err2);
	fprintf(stderr, "err3: %d\n", err3);
	fprintf(stderr, "err4: %d\n", err4);

	fprintf(stderr, "return_sites[0]: %d\n", return_sites[0]);
	fprintf(stderr, "return_sites[1]: %d\n", return_sites[1]);
	fprintf(stderr, "return_sites[2]: %d\n", return_sites[2]);
	fprintf(stderr, "return_sites[3]: %d\n", return_sites[3]);
	fprintf(stderr, "site: %d\n", site);

	fprintf(stderr, "*****************************\n\n\n\n\n");
	}
	**/






	if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}
	if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
	{
		return_sites[2] = site;
		return_sites_error[2] = err3;
	}
	if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
	{
		return_sites[3] = site;
		return_sites_error[3] = err4;
	}






	i = 0;



	while (i<errthold)
	{
		///err = err + ((VP >> i)&(Word_32)1);
		tmp_process = _mm_srli_epi32(VP, i);
		tmp_process = _mm_and_si128(tmp_process, err_mask);
		Err_4 = _mm_add_epi32(Err_4, tmp_process);

		///err = err - ((VN >> i)&(Word_32)1);
		tmp_process1 = _mm_srli_epi32(VN, i);
		tmp_process1 = _mm_and_si128(tmp_process1, err_mask);
		Err_4 = _mm_sub_epi32(Err_4, tmp_process1);
		++i;

		err1 = _mm_extract_epi32(Err_4, 0);
		err2 = _mm_extract_epi32(Err_4, 1);
		err3 = _mm_extract_epi32(Err_4, 2);
		err4 = _mm_extract_epi32(Err_4, 3);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}


	/**
	if (strcmp(text, "TTTTAGGTTTAAGTAATTTTTTTGTTTTAGTTTTTTAAGTAGTTGGGGTTATAGGTGTATGTTATTATNTTTAGTTAATTTTTT") == 0
	&&
	strcmp(pattern1, "TAGCTCACTGAAGCCTCAGATGATCCTCCCACCTCAGCCTCCCAAGTAGCTGGGGCTACAGGTGCATGCTACCACGACTGGCTAATTTTAATATTTTTA") == 0)
	{
	fprintf(stderr, "\n\n\n\n*****************************\ntext: %s\n", text);
	fprintf(stderr, "text: %s\n", pattern1);


	fprintf(stderr, "err1: %d\n", err1);
	fprintf(stderr, "err2: %d\n", err2);
	fprintf(stderr, "err3: %d\n", err3);
	fprintf(stderr, "err4: %d\n", err4);

	fprintf(stderr, "return_sites[0]: %d\n", return_sites[0]);
	fprintf(stderr, "return_sites[1]: %d\n", return_sites[1]);
	fprintf(stderr, "return_sites[2]: %d\n", return_sites[2]);
	fprintf(stderr, "return_sites[3]: %d\n", return_sites[3]);

	fprintf(stderr, "*****************************\n\n\n\n\n");
	}
	**/


	unsigned int ungap_err1;
	unsigned int ungap_err2;
	unsigned int ungap_err3;
	unsigned int ungap_err4;
	ungap_err1 = err1;
	ungap_err2 = err2;
	ungap_err3 = err3;
	ungap_err4 = err4;




	while (i<last_high)
	{
		///err = err + ((VP >> i)&(Word_32)1);
		tmp_process = _mm_srli_epi32(VP, i);
		tmp_process = _mm_and_si128(tmp_process, err_mask);
		Err_4 = _mm_add_epi32(Err_4, tmp_process);

		///err = err - ((VN >> i)&(Word_32)1);
		tmp_process1 = _mm_srli_epi32(VN, i);
		tmp_process1 = _mm_and_si128(tmp_process1, err_mask);
		Err_4 = _mm_sub_epi32(Err_4, tmp_process1);
		++i;

		err1 = _mm_extract_epi32(Err_4, 0);
		err2 = _mm_extract_epi32(Err_4, 1);
		err3 = _mm_extract_epi32(Err_4, 2);
		err4 = _mm_extract_epi32(Err_4, 3);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
		if ((err3 <= errthold) && (err3 <= return_sites_error[2]))
		{
			return_sites[2] = site + i;
			return_sites_error[2] = err3;
		}
		if ((err4 <= errthold) && (err4 <= return_sites_error[3]))
		{
			return_sites[3] = site + i;
			return_sites_error[3] = err4;
		}
	}


	/**
	if (strcmp(text, "TTTTAGGTTTAAGTAATTTTTTTGTTTTAGTTTTTTAAGTAGTTGGGGTTATAGGTGTATGTTATTATNTTTAGTTAATTTTTT") == 0
	&&
	strcmp(pattern1, "TAGCTCACTGAAGCCTCAGATGATCCTCCCACCTCAGCCTCCCAAGTAGCTGGGGCTACAGGTGCATGCTACCACGACTGGCTAATTTTAATATTTTTA") == 0)
	{
	fprintf(stderr, "\n\n\n\n*****************************\ntext: %s\n", text);
	fprintf(stderr, "text: %s\n", pattern1);


	fprintf(stderr, "err1: %d\n", err1);
	fprintf(stderr, "err2: %d\n", err2);
	fprintf(stderr, "err3: %d\n", err3);
	fprintf(stderr, "err4: %d\n", err4);

	fprintf(stderr, "return_sites[0]: %d\n", return_sites[0]);
	fprintf(stderr, "return_sites[1]: %d\n", return_sites[1]);
	fprintf(stderr, "return_sites[2]: %d\n", return_sites[2]);
	fprintf(stderr, "return_sites[3]: %d\n", return_sites[3]);

	fprintf(stderr, "*****************************\n\n\n\n\n");
	}
	**/



	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}


	if ((ungap_err3 <= errthold) && ungap_err3 == return_sites_error[2])
	{
		return_sites[2] = site + errthold;
	}


	if ((ungap_err4 <= errthold) && ungap_err4 == return_sites_error[3])
	{
		return_sites[3] = site + errthold;
	}


	/**
	if (strcmp(text, "TTTTAGGTTTAAGTAATTTTTTTGTTTTAGTTTTTTAAGTAGTTGGGGTTATAGGTGTATGTTATTATNTTTAGTTAATTTTTT") == 0
	&&
	strcmp(pattern1, "TAGCTCACTGAAGCCTCAGATGATCCTCCCACCTCAGCCTCCCAAGTAGCTGGGGCTACAGGTGCATGCTACCACGACTGGCTAATTTTAATATTTTTA") == 0)
	{
	fprintf(stderr, "\n\n\n\n*****************************\ntext: %s\n", text);
	fprintf(stderr, "text: %s\n", pattern1);


	fprintf(stderr, "err1: %d\n", err1);
	fprintf(stderr, "err2: %d\n", err2);
	fprintf(stderr, "err3: %d\n", err3);
	fprintf(stderr, "err4: %d\n", err4);

	fprintf(stderr, "return_sites[0]: %d\n", return_sites[0]);
	fprintf(stderr, "return_sites[1]: %d\n", return_sites[1]);
	fprintf(stderr, "return_sites[2]: %d\n", return_sites[2]);
	fprintf(stderr, "return_sites[3]: %d\n", return_sites[3]);

	fprintf(stderr, "*****************************\n\n\n\n\n");
	}
	**/




	return 1;

}



/*************************************************************************************************************************************************************************/




/***************************************************2���ȶԼ�sse**************************************************/
inline int BS_Reserve_Banded_BPM_2_SSE_only(char *pattern1, char *pattern2, int p_length, char *text, int t_length,
	int* return_sites, unsigned int* return_sites_error, unsigned short errthold, __m128i* Peq_SSE)

{

	memset(return_sites, -1, sizeof(int)* 2);
	memset(return_sites_error, -1, sizeof(unsigned int)* 2);

	Word Peq[256][2];
	int band_length = (errthold << 1) + 1;


	int i;

	Word tmp_Peq_1 = 1;


	memset(Peq['A'], 0, sizeof(Word)* 2);
	memset(Peq['C'], 0, sizeof(Word)* 2);
	memset(Peq['G'], 0, sizeof(Word)* 2);
	memset(Peq['T'], 0, sizeof(Word)* 2);

	for (i = 0; i<band_length; i++)
	{
		Peq[pattern1[i]][0] = Peq[pattern1[i]][0] | tmp_Peq_1;
		Peq[pattern2[i]][1] = Peq[pattern2[i]][1] | tmp_Peq_1;

		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq_SSE['A'] = _mm_set_epi64x(Peq['A'][1], Peq['A'][0]);
	Peq_SSE['C'] = _mm_set_epi64x(Peq['C'][1], Peq['C'][0]);
	Peq_SSE['G'] = _mm_set_epi64x(Peq['G'][1], Peq['G'][0]);
	Peq_SSE['T'] = _mm_set_epi64x(Peq['T'][1], Peq['T'][0]);

	Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);



	Word Mask = ((Word)1 << (errthold << 1));

	__m128i Mask1 = _mm_set_epi64x(0, Mask);
	__m128i Mask2 = _mm_set_epi64x(Mask, 0);



	__m128i VP = _mm_setzero_si128();
	__m128i VN = _mm_setzero_si128();
	__m128i X = _mm_setzero_si128();
	__m128i D0 = _mm_setzero_si128();
	__m128i HN = _mm_setzero_si128();
	__m128i HP = _mm_setzero_si128();
	__m128i tmp_process;
	__m128i tmp_process1;



	__m128i Err_2 = _mm_setzero_si128();
	__m128i err_mask = _mm_set_epi64x(1, 1);
	///for_not li quan shi 1
	__m128i for_not = _mm_set1_epi32(-1);
	__m128i err_arry;
	__m128i cmp_result;




	int i_bd = (errthold << 1);
	int last_high = (errthold << 1);
	int t_length_1 = t_length - 1;
	int err1;
	int err2;




	__m128i pre_end = _mm_set_epi64x(last_high + errthold, last_high + errthold);


	i = 0;
	while (i<t_length_1)
	{


		///X = Peq[text[i]] | VN;
		X = _mm_or_si128(Peq_SSE[text[i]], VN);



		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
		///X&VP
		tmp_process1 = _mm_and_si128(X, VP);
		///(VP + (X&VP))
		tmp_process = _mm_add_epi64(tmp_process1, VP);
		///((VP + (X&VP)) ^ VP)
		tmp_process = _mm_xor_si128(tmp_process, VP);
		///((VP + (X&VP)) ^ VP) | X
		D0 = _mm_or_si128(tmp_process, X);
		/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

		///HN = VP&D0;
		HN = _mm_and_si128(D0, VP);

		///HP = VN | ~(VP | D0);
		tmp_process = _mm_or_si128(D0, VP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		HP = _mm_or_si128(tmp_process, VN);


		///X = D0 >> 1;
		X = _mm_srli_epi64(D0, 1);
		///VN = X&HP;
		VN = _mm_and_si128(X, HP);
		///VP = HN | ~(X | HP);
		tmp_process = _mm_or_si128(X, HP);
		tmp_process = _mm_andnot_si128(tmp_process, for_not);
		VP = _mm_or_si128(HN, tmp_process);

		///D0&err_mask
		err_arry = _mm_and_si128(D0, err_mask);
		Err_2 = _mm_add_epi64(Err_2, err_mask);
		Err_2 = _mm_sub_epi64(Err_2, err_arry);



		Peq_SSE['A'] = _mm_srli_epi64(Peq_SSE['A'], 1);
		Peq_SSE['T'] = _mm_srli_epi64(Peq_SSE['T'], 1);
		Peq_SSE['G'] = _mm_srli_epi64(Peq_SSE['G'], 1);
		Peq_SSE['C'] = _mm_srli_epi64(Peq_SSE['C'], 1);

		++i;
		++i_bd;

		Peq_SSE[pattern1[i_bd]] = _mm_or_si128(Mask1, Peq_SSE[pattern1[i_bd]]);
		Peq_SSE[pattern2[i_bd]] = _mm_or_si128(Mask2, Peq_SSE[pattern2[i_bd]]);



		Peq_SSE['T'] = _mm_or_si128(Peq_SSE['T'], Peq_SSE['C']);

	}







	///X = Peq[text[i]] | VN;
	X = _mm_or_si128(Peq_SSE[text[i]], VN);

	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/
	///X&VP
	tmp_process1 = _mm_and_si128(X, VP);
	///(VP + (X&VP))
	tmp_process = _mm_add_epi64(tmp_process1, VP);
	///((VP + (X&VP)) ^ VP)
	tmp_process = _mm_xor_si128(tmp_process, VP);
	///((VP + (X&VP)) ^ VP) | X
	D0 = _mm_or_si128(tmp_process, X);
	/*************D0 = ((VP + (X&VP)) ^ VP) | X*********************/

	///HN = VP&D0;
	HN = _mm_and_si128(D0, VP);

	///HP = VN | ~(VP | D0);
	tmp_process = _mm_or_si128(D0, VP);
	tmp_process = _mm_andnot_si128(tmp_process, for_not);
	HP = _mm_or_si128(tmp_process, VN);


	///X = D0 >> 1;
	X = _mm_srli_epi64(D0, 1);
	///VN = X&HP;
	VN = _mm_and_si128(X, HP);
	///VP = HN | ~(X | HP);
	tmp_process = _mm_or_si128(X, HP);
	tmp_process = _mm_andnot_si128(tmp_process, for_not);
	VP = _mm_or_si128(HN, tmp_process);

	///D0&err_mask
	err_arry = _mm_and_si128(D0, err_mask);
	Err_2 = _mm_add_epi64(Err_2, err_mask);
	Err_2 = _mm_sub_epi64(Err_2, err_arry);

	///shi ji shang zhe ge zhi hen xiao d
	cmp_result = _mm_cmpgt_epi64(Err_2, pre_end);

	///jian zhi
	if (_mm_extract_epi64(cmp_result, 0) && _mm_extract_epi64(cmp_result, 1))
		return 1;








	int site = t_length - 1;
	err1 = _mm_extract_epi64(Err_2, 0);
	err2 = _mm_extract_epi64(Err_2, 1);




	if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
	{
		return_sites[0] = site;
		return_sites_error[0] = err1;
	}
	if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
	{
		return_sites[1] = site;
		return_sites_error[1] = err2;
	}






	i = 0;



	while (i<errthold)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm_srli_epi64(VP, i);
		tmp_process = _mm_and_si128(tmp_process, err_mask);
		Err_2 = _mm_add_epi64(Err_2, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm_srli_epi64(VN, i);
		tmp_process1 = _mm_and_si128(tmp_process1, err_mask);
		Err_2 = _mm_sub_epi64(Err_2, tmp_process1);
		++i;

		err1 = _mm_extract_epi64(Err_2, 0);
		err2 = _mm_extract_epi64(Err_2, 1);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}

	}





	unsigned int ungap_err1;
	unsigned int ungap_err2;

	ungap_err1 = err1;
	ungap_err2 = err2;





	while (i<last_high)
	{
		///err = err + ((VP >> i)&(Word)1);
		tmp_process = _mm_srli_epi64(VP, i);
		tmp_process = _mm_and_si128(tmp_process, err_mask);
		Err_2 = _mm_add_epi64(Err_2, tmp_process);

		///err = err - ((VN >> i)&(Word)1);
		tmp_process1 = _mm_srli_epi64(VN, i);
		tmp_process1 = _mm_and_si128(tmp_process1, err_mask);
		Err_2 = _mm_sub_epi64(Err_2, tmp_process1);
		++i;

		err1 = _mm_extract_epi64(Err_2, 0);
		err2 = _mm_extract_epi64(Err_2, 1);


		if ((err1 <= errthold) && (err1 <= return_sites_error[0]))
		{
			return_sites[0] = site + i;
			return_sites_error[0] = err1;
		}
		if ((err2 <= errthold) && (err2 <= return_sites_error[1]))
		{
			return_sites[1] = site + i;
			return_sites_error[1] = err2;
		}
	}






	if ((ungap_err1 <= errthold) && ungap_err1 == return_sites_error[0])
	{
		return_sites[0] = site + errthold;
	}


	if ((ungap_err2 <= errthold) && ungap_err2 == return_sites_error[1])
	{
		return_sites[1] = site + errthold;
	}



	return 1;

}



/*************************************************************************************************************************************************************************/









#endif
