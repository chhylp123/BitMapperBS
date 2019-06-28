// This is a demo program for showing how to call SACA_K.

#include<stdio.h>
#include<stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "bwt.h"
#include<stdint.h>
#include<ctype.h>
#include <unistd.h>
///#include <nmmintrin.h>

bwt_index bitmapper_index_params;



SA_block build_sa_block;

unsigned int debug_2 = 0;

unsigned int debug_total = 0;




void dec_bit_bwt(bwt_string_type input)
{
	unsigned int bit_2[64];
	int i = 0;
	for (i = 0; i < 64; i++)
	{

		bwt_string_type tmp = (bwt_string_type)((bwt_string_type)input&(bwt_string_type)1);

		input = input >> 1;

		if (tmp == (bwt_string_type)0)
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


void SACA_K(unsigned char *s, unsigned int *SA, unsigned int n,
           unsigned int K, unsigned int m, int level);

inline void bwt_direct_get_sa_interval_long_back_up_1
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int j_sp, j_ep;
	unsigned int l;
	unsigned int i = 0;

	unsigned int last_sp, last_ep;

	unsigned int actually_line_sp, actually_line_ep;

	unsigned int ans_sp, ans_ep, flag;


	SA_flag_string_type tmp_SA_pop_count;

	flag = 0;


	l = sp;
	///这是除以224，乘以256，最后再除以64，这就是确定属于哪个block
	actually_line_sp = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;
	///这是对224取余，这就是确定block内部还需要扫描的位数
	last_sp = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;



	l = ep;
	///这是除以224，乘以256，最后再除以64，这就是确定属于哪个block
	actually_line_ep = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;
	///这是除以224，乘以256，最后再除以64，这就是确定属于哪个block
	last_ep = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;




	///这说明sp和ep在同一个block内部,则count操作可以部分重叠
	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}








	ans_sp = SA_flag[actually_line_sp] >> bitmapper_index_params.SA_counter_shift_length;

	if (last_sp != 0)
	{
		///主要last本身已经加过SA_counter_length了, 所以这里j=0, 否则就应该把j = SA_counter_length
		j_sp = 0;

		tmp_SA_pop_count = SA_flag[actually_line_sp] & bitmapper_index_params.SA_pop_count_mode;

		while (j_sp + bitmapper_index_params.SA_flag_warp_number <= last_sp)
		{

			ans_sp = ans_sp + __builtin_popcountll(tmp_SA_pop_count);

			j_sp = j_sp + bitmapper_index_params.SA_flag_warp_number;

			actually_line_sp++;
			tmp_SA_pop_count = SA_flag[actually_line_sp];

		}

		///这个判断条件貌似是sp和ep都属于同一个block
		///所以ep的值貌似可以直接继承
		if (flag == 1)
		{
			j_ep = j_sp;

			ans_ep = ans_sp;

			actually_line_ep = actually_line_sp;

			flag = 2;
		}


		last_sp = last_sp - j_sp;


		if (last_sp != 0)
		{
			ans_sp = ans_sp +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_sp))));
		}


	}




	(*start) = ans_sp;





	if (last_ep != 0)
	{

		if (flag==2)
		{

			tmp_SA_pop_count = SA_flag[actually_line_ep];

			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}
		}
		else
		{


			ans_ep = SA_flag[actually_line_ep] >> bitmapper_index_params.SA_counter_shift_length;

			///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
			j_ep = 0;

			tmp_SA_pop_count = SA_flag[actually_line_ep] & bitmapper_index_params.SA_pop_count_mode;


			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}

		}







	}
	else
	{
		ans_ep = SA_flag[actually_line_ep] >> bitmapper_index_params.SA_counter_shift_length;
	}






	(*length) = ans_ep - (*start);

}


inline void bwt_direct_get_sa_interval_long
(SA_flag_string_type* SA_flag,
 bitmapper_bs_iter sp,
 bitmapper_bs_iter ep,
 bitmapper_bs_iter* start,
 bitmapper_bs_iter* length)
{
	bitmapper_bs_iter last_sp, last_ep;

	bitmapper_bs_iter actually_line_sp, actually_line_ep;

	bitmapper_bs_iter ans_sp, ans_ep;

	bitmapper_bs_iter j_sp, j_ep;

	bitmapper_bs_iter flag = 0;

	SA_flag_string_type tmp_SA_pop_count;


    actually_line_sp = ((sp >> 8 ) * bitmapper_index_params.acctuall_SA_flag_gap) >> 6;
	last_sp = sp & 255;

	actually_line_ep = ((ep >> 8 ) * bitmapper_index_params.acctuall_SA_flag_gap) >> 6;
	last_ep = ep & 255;


	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}



	ans_sp = SA_flag[actually_line_sp];


	if (last_sp != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		j_sp = 0;

		actually_line_sp++;
		tmp_SA_pop_count = SA_flag[actually_line_sp];




		while (j_sp + bitmapper_index_params.SA_flag_warp_number <= last_sp)
		{

			ans_sp = ans_sp + __builtin_popcountll(tmp_SA_pop_count);

			j_sp = j_sp + bitmapper_index_params.SA_flag_warp_number;

			actually_line_sp++;
			tmp_SA_pop_count = SA_flag[actually_line_sp];

		}


		if (flag == 1)
		{
			ans_ep = ans_sp;

			j_ep = j_sp;

			actually_line_ep = actually_line_sp;

			flag = 2;
		}



		last_sp = last_sp - j_sp;


		if (last_sp != 0)
		{
			ans_sp = ans_sp +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_sp))));
		}


	}




	(*start) = ans_sp;










	if (last_ep != 0)
	{

		if (flag == 2)
		{

			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}



		}
		else
		{
			ans_ep = SA_flag[actually_line_ep];
			actually_line_ep++;


			j_ep = 0;
			tmp_SA_pop_count = SA_flag[actually_line_ep];




			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}

		}






	}
	else
	{
		ans_ep = SA_flag[actually_line_ep];
	}




	(*length) = ans_ep - (*start);

}




inline void bwt_direct_get_sa_interval_long_back_up
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	SA_flag_string_type tmp_SA_pop_count;


	l = sp;

	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
	*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;


	///这里SA要改
	///last = l % compress_SA_flag;
	///last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;



	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}




	(*start) = ans;

	l = ep;



	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
	*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/
	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;



	///这里SA要改
	///last = l % compress_SA_flag;
	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;

	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}




	(*length) = ans - (*start);

}




inline unsigned int bwt_get_sa_restrict_steps_less_than_4
(SA_flag_string_type* SA_flag, unsigned int line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;


	if (line == bitmapper_index_params.shapline) return 0;
	unsigned int actually_line;
	unsigned int last;

	/**
	///这里SA要改
	unsigned int actually_line = ((l / bitmapper_index_params.compress_SA_flag)
	*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;


	///这里SA要改
	///unsigned int last = l % compress_SA_flag;
	unsigned int last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;


	actually_line = actually_line + last / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;

	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	actually_line = actually_line + (last >> 6);




	/**
	while (((SA_flag[actually_line] << (last % bitmapper_index_params.SA_flag_warp_number))
	&bitmapper_index_params.mode_SA_flag) == 0)
	**/
	while (((SA_flag[actually_line] << (last & 63))&bitmapper_index_params.mode_SA_flag) == 0)
	{

		if (i >= need_step)
		{
			return 0;
		}



		/**
		///似乎只有这里要改，因为这里涉及到取BWT了
		///这里做了修改
		///actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;
		actually_line = ((l / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.acctuall_bwt_gap) + l % bitmapper_index_params.compress_occ + 32;
		delta = (bwt[actually_line / bitmapper_index_params.bwt_warp_number]
		>> ((bitmapper_index_params.bwt_warp_number - actually_line%bitmapper_index_params.bwt_warp_number - 1) * 2))
		&(bwt_string_type)3;
		**/

		actually_line = ((l >> 8)*bitmapper_index_params.acctuall_bwt_gap) + (l & 255) + 32;

		/**
		delta = (bwt[(actually_line >> 5)]
		>> ((bitmapper_index_params.bwt_warp_number - (actually_line & 31) - 1) * 2))
		&(bwt_string_type)3;
		**/

		delta = (bwt[(actually_line >> 5)]
			>> ((bitmapper_index_params.bwt_warp_number - (actually_line & 31) - 1) << 1))
			&(bwt_string_type)3;







		l = find_occ_fm_index(l, delta, bwt, high_occ_table);


		i++;

		if (l == bitmapper_index_params.shapline)
		{
			(*accessed_sa) = i;

			return i;
		}




		/**
		///这里SA要改
		///actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;
		actually_line = ((l / bitmapper_index_params.compress_SA_flag)*bitmapper_index_params.acctuall_SA_flag_gap)
		/ bitmapper_index_params.SA_flag_warp_number;

		///这里SA要改
		///last = l % compress_SA_flag;
		last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

		actually_line = actually_line + last / bitmapper_index_params.SA_flag_warp_number;

		**/
		actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;

		last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

		actually_line = actually_line + (last >> 6);

	}

	/**
	actually_line
	= ((l / bitmapper_index_params.compress_SA_flag)*
	bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/

	actually_line = actually_line - (last >> 6);


	///这里SA要改
	///unsigned int ans = SA_flag[actually_line];
	unsigned int ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;


	///这里SA要改
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;
		SA_flag_string_type tmp_SA_pop_count;
		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}

	(*accessed_sa) = sa[ans] + i;

	return 1;
}






inline void bwt_accesss_SA_cur_more_than_3(
    unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter occ,
	bitmapper_bs_iter* tmp_SA_length,
	bitmapper_bs_iter need_step)
{

	bitmapper_bs_iter t;
	for (t = sp; t<ep; t++)
	{
		///这里SA要改
		if (bwt_get_sa_restrict_steps_more_than_3(SA_flag, t, bwt, high_occ_table, sa,
			need_step, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;
		}
	}

}




inline void bwt_accesss_SA_cur_less_than_4(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length, unsigned int need_step)
{

	unsigned int t;
	for (t = sp; t<ep; t++)
	{
		///这里SA要改
		if (bwt_get_sa_restrict_steps_less_than_4(SA_flag, t, bwt, high_occ_table, sa,
			need_step, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;
		}

	}

}


void read_n_uint40_from_file(uint40* tab, long length, char* fname)
{
	FILE *f = fopen(fname, "r");
	size_t fread_ret = fread(tab, sizeof(uint40), length, f);
	if ((long)fread_ret != length) {
		fprintf(stderr, "\nError: fread: fread_ret=%lld, length=%lld\n",
			fread_ret, length);
		exit(EXIT_FAILURE);
	}
	fclose(f);
}

void pure_read_n_uint40_from_file(uint40* tab, long length, bitmapper_bs_iter start)
{
	fseek(build_sa_block.f_sa, start*sizeof(uint40), SEEK_SET);

	size_t fread_ret = fread(tab, sizeof(uint40), length, build_sa_block.f_sa);
}


void pSAscan_indenpendent_get_sa_fromFILE(uint40 **sa, unsigned int text_length, char *refer)
{

	FILE* tmp_SA = fopen("tmp_ref.tmp", "w");

	fwrite(refer, sizeof(char), text_length, tmp_SA);

	fflush(tmp_SA);

	fclose(tmp_SA);


	int error=system("./psascan tmp_ref.tmp -m 8192");

	error = system("rm tmp_ref.tmp");

	////fprintf(stderr, "error=%d\n", error);

	uint40* SA = new uint40[text_length + 1];
	SA[0] = uint40(text_length, 0);

	read_n_uint40_from_file(SA+1, text_length, "tmp_ref.tmp.sa5");

	*sa = SA;

	fclose(tmp_SA);


}





uint40 get_sa_block(bitmapper_bs_iter iter)
{
	///fprintf(stderr, "iter=%llu, build_sa_block.start=%llu, build_sa_block.endl=%llu\n", iter, build_sa_block.start, build_sa_block.endl);


	if (iter >= build_sa_block.start
		&&
		iter <= build_sa_block.endl)
	{
		iter = iter - build_sa_block.start;
		return build_sa_block.sa_block[iter];
	}
	else
	{
		if (iter == 0)
		{
			build_sa_block.sa_block[0] = uint40(build_sa_block.text_length);
			build_sa_block.start = 0;
			build_sa_block.endl = build_sa_block.start + build_sa_block.block_size - 1;
			pure_read_n_uint40_from_file(build_sa_block.sa_block+1, build_sa_block.block_size - 1, 0);
		}
		else
		{
			build_sa_block.start = iter;
			build_sa_block.endl = build_sa_block.start + build_sa_block.block_size - 1;
			pure_read_n_uint40_from_file(build_sa_block.sa_block, build_sa_block.block_size, iter - 1);
		}

		iter = iter - build_sa_block.start;
		return build_sa_block.sa_block[iter];

	}

}

void init_SA_block(bitmapper_bs_iter text_length)
{
	build_sa_block.block_size = 1000000;
	build_sa_block.text_length = text_length;
	build_sa_block.sa_block = new uint40[build_sa_block.block_size];
	build_sa_block.f_sa = fopen("tmp_ref.tmp.sa5", "r");


	bitmapper_bs_iter file_length = 0;
	fseek(build_sa_block.f_sa, 0, SEEK_END);
	file_length = ftell(build_sa_block.f_sa);
	fseek(build_sa_block.f_sa, 0, SEEK_SET);
	file_length = file_length / sizeof(uint40);

	if (text_length == file_length)
	{
		fprintf(stderr, "text_length =%llu, file_length=%llu\n", text_length, file_length);
		fprintf(stderr, "text_length = file_length...sucess...\n");
		fflush(stderr);
	}
	else
	{
		fprintf(stderr, "text_length =%llu, file_length=%llu\n", text_length, file_length);
		fprintf(stderr, "ERROR: text_length != file_length...exit...\n");
		fflush(stderr);
		exit(0);
	}



	build_sa_block.sa_block[0] = uint40(text_length);
	build_sa_block.start = 0;
	build_sa_block.endl = build_sa_block.start + build_sa_block.block_size - 1;
	pure_read_n_uint40_from_file(build_sa_block.sa_block+1, build_sa_block.block_size-1, 0);


}


void indenpendent_get_sa_fromFILE(unsigned int **sa, unsigned int cc, char *refer)
{
	unsigned int i;
	unsigned int n;
	n = cc-1;

	fprintf(stderr, "r_n=%u\n", n);

	n++; // append the virtual sentinel
	fprintf(stderr, "Allocating input and output space: %u bytes = %.2lf MB", 5 * n, (double)5 * n / 1024 / 1024);
	unsigned char *s_ch = new unsigned char[n];
	//unsigned char *s_ch;
	unsigned int *SA = new unsigned int[n];
	//unsigned int *SA;
	///SA = (unsigned int*)malloc(sizeof(unsigned int)*n+10);
	if (s_ch == NULL || SA == NULL) {
		delete[] s_ch; delete[] SA;
		fprintf(stderr, "\nInsufficient memory, exit!");
		return;
	}


	// read the string into buffer.
	fprintf(stderr, "\nReading input string...");
	fseek(stdin, 0, SEEK_SET);
	for(i=0;i<n;i++) s_ch[i]=refer[i];
	// set the virtual sentinel
	s_ch[n - 1] = 0;








	clock_t start, finish;
	double  duration;
	start = clock();

	fprintf(stderr, "\nConstructing the suffix array...");
	SACA_K(s_ch, SA, n, 256, n, 0);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;

	fprintf(stderr, "\nSize: %u bytes, Time: %5.3f seconds\n", n - 1, duration);

	SA[0] = n - 1;
	printf("sa_n=%d\n",n);
	*sa = SA;
}






void new_version_pSAscan_build_sa(bitmapper_bs_iter text_length, char **refer)
{

	FILE* tmp_SA = fopen("tmp_ref.tmp", "w");

	fwrite(*refer, sizeof(char), text_length, tmp_SA);

	fflush(tmp_SA);

	fclose(tmp_SA);

	free(*refer);

	int error;
	
	error = system("rm tmp_ref.tmp.sa5");


	if (access("psascan", F_OK) == 0)
	{
		///文件存在
		fprintf(stderr, "the binary of psascan exists...\n");
		error = system("./psascan tmp_ref.tmp -m 8192");
	}
	else
	{
		///文件不存在
		fprintf(stderr, "the binary of psascan does not exist...\n");
		error = system("psascan tmp_ref.tmp -m 8192");
	}

	


	tmp_SA = fopen("tmp_ref.tmp", "r");
	*refer = (char*)malloc(sizeof(char)*(text_length + 1));
	fread(*refer, sizeof(char), text_length, tmp_SA);
	fclose(tmp_SA);




	error = system("rm tmp_ref.tmp");

	init_SA_block(text_length);

}




inline void string_match_scan(char* pattern, bitmapper_bs_iter pattern_length, bitmapper_bs_iter length, char *refer)
{

	char bit_to_car[3];

	bitmapper_bs_iter number_occ = 0;

	bit_to_car[0] = 'G';
	bit_to_car[1] = 'T';
	bit_to_car[2] = 'A';


	bitmapper_bs_iter i = 0;
	for (i = 0; i < length - pattern_length + 1; i++)
	{
		bitmapper_bs_iter j = 0;
		bitmapper_bs_iter unmatched = 0;
		for (j = 0; j < pattern_length; j++)
		{

			if (i == 209714077)
			{
				fprintf(stderr, "bit_to_car[refer[i + j]]: %c \n", bit_to_car[refer[i + j]]);
			}

			if (pattern[j] != bit_to_car[refer[i + j]])
			{
				unmatched = 1;
				break;

			}

		}

		if (unmatched==0)
		{
			number_occ++;
			fprintf(stderr, "scan: match site: %llu \n", i);
		}
	}

	fprintf(stderr, "scan: number_occ: %llu \n", number_occ);
}

unsigned int indenpendent_creadte_index(bitmapper_bs_iter text_length, char** input_refer,
	bitmapper_bs_iter compress_sa, char* filename)
{

	bitmapper_bs_iter nacgt[5];
	unsigned int ctoi[256];

	char filenames[100], filenameo[100], filenameb[100];
	char *refer = (*input_refer);
	bwt_string_type *bwt;
	uint40 *sa;
	char ch;
	bitmapper_bs_iter i, j;
	///这里要改
	FILE *f1, *f2, *fs, *fb, *fo;



	/**
    ctoi['A'] = 0;
	ctoi['G'] = 2;
	ctoi['T'] = 1;
	ctoi['a'] = 0;
	ctoi['g'] = 2;
	ctoi['t'] = 1;
	**/

	ctoi['A'] = 2;
	ctoi['G'] = 0;
	ctoi['T'] = 1;
	ctoi['a'] = 2;
	ctoi['g'] = 0;
	ctoi['t'] = 1;




	///this is the number of occ line
	bitmapper_bs_iter occ_line_number = (text_length + 1) / bitmapper_index_params.compress_occ + 1;
	bitmapper_bs_iter high_occ_line_number  = (text_length + 1) / bitmapper_index_params.high_compress_occ + 1;

	///this is the number of character which could be saved in one bwt_string_type
	///bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	///this is the number of bwt_string_type representing bwt string (does not include occ)
	bitmapper_bs_iter bwt_length = (text_length + 1) / bitmapper_index_params.bwt_warp_number + 1;

	bitmapper_bs_iter occ_byte_length
		= (occ_line_number * 2 * sizeof(unsigned short))/ (sizeof(bwt_string_type)) + 1;



	
	bitmapper_bs_iter SA_flag_length = text_length + 1;
	bitmapper_bs_iter SA_flag_occ_length = (text_length + 1) / bitmapper_index_params.compress_SA_flag + 1;

	bitmapper_bs_iter SA_flag_byte_length = (SA_flag_length + SA_flag_occ_length * SA_counter_length + 1)
		/ bitmapper_index_params.SA_flag_warp_number + 1;



	







    
	///这里要改
	high_occ_table_type* high_occ_table;
	bitmapper_bs_iter high_occ_table_length = ((text_length + 1) / bitmapper_index_params.idependent_high_compress_occ + 1) * 2 + 1;
	


	
	





	bitmapper_index_params.SA_length = text_length + 1;
	///这个应该不需要吧
	for (i = 0; i < text_length; i++)
	{


		ch = refer[i];

		if(
            ch == 'A' ||
            ch == 'G' ||
            ch == 'T' ||
            ch == 'a' ||
            ch == 'g' ||
            ch == 't'
          )
        {
            refer[i] = ctoi[ch];
        }
        else
        {
            fprintf(stderr, "Text includes character which dose not belong to {A, G, T, a, g, t}! FMtree will exit ...\n");
            return 1;
        }


	}




	strcpy(filename + strlen(filename), ".index");
	strcpy(filenames, filename);
	strcpy(filenameo, filename);
	strcpy(filenameb, filename);
	strcpy(filenames + strlen(filename), ".sa");
	strcpy(filenameo + strlen(filename), ".occ");
	strcpy(filenameb + strlen(filename), ".bwt");
	f2 = fopen(filename, "w");
	fs = fopen(filenames, "w");
	fb = fopen(filenameb, "w");
	fo = fopen(filenameo, "w");



	///这个应该不需要吧
	///indenpendent_get_sa_fromFILE(&sa, SA_length, refer);
	////pSAscan_indenpendent_get_sa_fromFILE(&sa, text_length, refer);


	///printf("first SA has been generated!\n");
	new_version_pSAscan_build_sa(text_length, &refer);



	printf("SA has been generated!\n");













	bwt = (bwt_string_type *)malloc(sizeof(bwt_string_type)*(bwt_length + occ_byte_length + 1));
	SA_flag_string_type* SA_flag = (SA_flag_string_type *)malloc(sizeof(SA_flag_string_type)*SA_flag_byte_length);
	high_occ_table = (high_occ_table_type*)malloc(sizeof(high_occ_table_type)*high_occ_table_length);


	///这里SA要改
	///for (i = 0; i < SA_flag_length + SA_occ_byte_length + 1; i++)
	for (i = 0; i < SA_flag_byte_length; i++)
	{
		SA_flag[i] = (SA_flag_string_type)0;
	}

	for (i = 0; i < bwt_length + occ_byte_length + 1; i++)
	{
		bwt[i] = (bwt_string_type)0;
	}


	fprintf(stderr, "bwt_length=%llu,occ_byte_length=%llu, bwt_length+occ_byte_length=%llu\n",
		bwt_length, occ_byte_length, bwt_length + occ_byte_length);























    for (j = 0; j <= 4; j++)
		nacgt[j] = 0;//nacgt[0][j] = 0;


	bitmapper_index_params.bwt_length = 0;///bitmapper_bs_iter bwt_iterater = 0;
	i = 0;
	bitmapper_bs_iter shift_length = 0;
	bwt_string_type tmp_bwt_0 = (bwt_string_type)0;
	bwt_string_type tmp_bwt_1 = (bwt_string_type)0;


    ///bwt的第一个元素存的是两个32-bit的counter值，这个肯定是0
	///但是第二个元素存的是6个8-bit的counter值，这个肯定不是0
	///这个要注意
	///bitmapper_bs_iter occ_bwt_number = 2;
	bitmapper_bs_iter occ_bwt_number = 1;


	for (bitmapper_index_params.bwt_length = 0; bitmapper_index_params.bwt_length < occ_bwt_number; bitmapper_index_params.bwt_length++)///for (bwt_iterater = 0; bwt_iterater < occ_bwt_number; bwt_iterater++)
	{
		bwt[bitmapper_index_params.bwt_length] = (bwt_string_type)0;///bwt[bwt_iterater] = (bwt_string_type)0;
	}


	///这里要改
	bitmapper_index_params.high_occ_table_length = 0;
	for (bitmapper_index_params.high_occ_table_length = 0; bitmapper_index_params.high_occ_table_length < 2; bitmapper_index_params.high_occ_table_length++)
	{
		high_occ_table[bitmapper_index_params.high_occ_table_length]=0;
	}
	


	i = 0;

	bitmapper_bs_iter bwt_iter_i = 0;///long long bwt_iter_i = 0;
	bitmapper_bs_iter high_occ_iter = 0;
	bitmapper_bs_iter low_occ_iter = 0;





	while (1)
	{


		tmp_bwt_0 = (bwt_string_type)0;
		tmp_bwt_1 = (bwt_string_type)0;

		if (i >= bitmapper_index_params.SA_length)
		{
			break;
		}


		if (uint64_t(get_sa_block(i)) != 0)
		{
			ch = refer[uint64_t(get_sa_block(i)) - 1];
		}
		else
		{
			bitmapper_index_params.shapline = i;

			i++;

			continue;
		}



		///shift_length = (bitmapper_index_params.bwt_warp_number - bwt_iter_i % bitmapper_index_params.bwt_warp_number - 1) * 2;
        shift_length = 64 - bwt_iter_i % 64 - 1;


		/**
		ctoi['A'] = 0;
		ctoi['G'] = 2;
		ctoi['T'] = 1;
		ctoi['a'] = 0;
		ctoi['g'] = 2;
		ctoi['t'] = 1;
		**/

        ///tmp_bwt_0实际上存的是ch的低位，tmp_bwt_1实际上存的是ch的高位
		tmp_bwt_0 = (bwt_string_type)(ch);
		tmp_bwt_0 = tmp_bwt_0 & (bwt_string_type)1;
		tmp_bwt_0 = tmp_bwt_0 << shift_length;

		tmp_bwt_1 = (bwt_string_type)(ch);
		tmp_bwt_1 = tmp_bwt_1 >> 1;
		tmp_bwt_1 = tmp_bwt_1 & (bwt_string_type)1;
		tmp_bwt_1 = tmp_bwt_1 << shift_length;





		///bwt[bitmapper_index_params.bwt_length] = bwt[bitmapper_index_params.bwt_length] | tmp_bwt;
		///实际是低位存在前一个64-bit里，高位存在后一个64-bit里
        bwt[bitmapper_index_params.bwt_length]
            = bwt[bitmapper_index_params.bwt_length] | tmp_bwt_0;

        bwt[bitmapper_index_params.bwt_length + 1]
            = bwt[bitmapper_index_params.bwt_length + 1] | tmp_bwt_1;




		///tmp_bwt = (bwt_string_type)0;
		tmp_bwt_0 = (bwt_string_type)0;
		tmp_bwt_1 = (bwt_string_type)0;


		///nacgt[0]肯定是$
		nacgt[ch + 1]++;



		i++;
		bwt_iter_i++;


		if (bwt_iter_i % 64 == 0)
		{
			///bitmapper_index_params.bwt_length++;///bwt_iterater++;
			///注意这里是+2啊
			///因为ch的低位和高位被分别存在了第一个64-bit和第二个64-bit中
			bitmapper_index_params.bwt_length = bitmapper_index_params.bwt_length + 2;
		}


		///我们实际上只存T和G的counter值
		if (bwt_iter_i % bitmapper_index_params.idependent_high_compress_occ == 0)
		{


			high_occ_table[bitmapper_index_params.high_occ_table_length] = nacgt[2];
			bitmapper_index_params.high_occ_table_length++;



			high_occ_table[bitmapper_index_params.high_occ_table_length] = nacgt[3];
			bitmapper_index_params.high_occ_table_length++;
		}



		///我们实际上只存T和G的counter值
		if (bwt_iter_i % bitmapper_index_params.high_compress_occ == 0)
		{
			////注意这里实际上就是要+1而已
			bitmapper_index_params.bwt_length++;
		}





		if (bwt_iter_i % bitmapper_index_params.compress_occ == 0)
		{

			///这个实际上是先取两个40-bit的counter值
			high_occ_iter = ((bwt_iter_i / bitmapper_index_params.high_compress_occ) * bitmapper_index_params.acctuall_bwt_gap)
				/ bitmapper_index_params.bwt_warp_number;


			low_occ_iter = (bwt_iter_i % bitmapper_index_params.high_compress_occ) / bitmapper_index_params.compress_occ;



			bitmapper_bs_iter tmp_counter_high32 = (bwt_string_type)0;
			bitmapper_bs_iter tmp_counter_low8 = (bwt_string_type)0;
			bitmapper_bs_iter tmp_counter_40 = (bwt_string_type)0;
			bitmapper_bs_iter mode_low8 = (bwt_string_type)255;
			bitmapper_bs_iter mode_low32 = ((bwt_string_type)-1) >> 32;


			if (low_occ_iter == 0)
			{
				/***********************************************字符1******************************************************/
				tmp_counter_40 = high_occ_table[(bwt_iter_i / bitmapper_index_params.idependent_high_compress_occ) * 2];
				///T实际上被转成了数字1，所以T的counter值实际是存在了nacgt[2]中
				tmp_counter_40 = nacgt[2] - tmp_counter_40;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				tmp_counter_40 = tmp_counter_40 << 48;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				bwt[high_occ_iter] = bwt[high_occ_iter] | tmp_counter_40;
				/***********************************************字符1******************************************************/

				/***********************************************字符2******************************************************/
				tmp_counter_40 = high_occ_table[(bwt_iter_i / bitmapper_index_params.idependent_high_compress_occ) * 2 + 1];
				tmp_counter_40 = nacgt[3] - tmp_counter_40;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				tmp_counter_40 = tmp_counter_40 << 32;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				bwt[high_occ_iter] = bwt[high_occ_iter] | tmp_counter_40;
				/***********************************************字符2******************************************************/
			}
			else if (low_occ_iter == 1)
			{

				/***********************************************字符1******************************************************/
				tmp_counter_40 = high_occ_table[(bwt_iter_i / bitmapper_index_params.idependent_high_compress_occ) * 2];
				///T实际上被转成了数字1，所以T的counter值实际是存在了nacgt[2]中
				tmp_counter_40 = nacgt[2] - tmp_counter_40;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				tmp_counter_40 = tmp_counter_40 << 16;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				bwt[high_occ_iter] = bwt[high_occ_iter] | tmp_counter_40;
				/***********************************************字符1******************************************************/

				/***********************************************字符2******************************************************/
				tmp_counter_40 = high_occ_table[(bwt_iter_i / bitmapper_index_params.idependent_high_compress_occ) * 2 + 1];
				tmp_counter_40 = nacgt[3] - tmp_counter_40;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				tmp_counter_40 = tmp_counter_40;
				/************************这里是不同low_occ_iter所不同的地方******************************/
				bwt[high_occ_iter] = bwt[high_occ_iter] | tmp_counter_40;
				/***********************************************字符2******************************************************/

			}

		}


		tmp_bwt_0 = (bwt_string_type)0;
		tmp_bwt_1 = (bwt_string_type)0;

	}




	/**
	///if ( bwt_iter_i%bwt_warp_number != 0 && bwt_iter_i % compress_occ != 0)
	if (bwt_iter_i%bitmapper_index_params.bwt_warp_number != 0 && bwt_iter_i % bitmapper_index_params.compress_occ != 0)
	{
		bitmapper_index_params.bwt_length++;//bwt_iterater++;
	}
	**/
	bitmapper_index_params.bwt_length = bitmapper_index_params.bwt_length + 2;

	printf("Occ and BWT have been built!\n");



	fprintf(stderr, "write SA_length=%u, shapline=%u\n", bitmapper_index_params.SA_length, bitmapper_index_params.shapline);








	SA_flag_string_type tmp_SA_flag = (SA_flag_string_type)0;

	bitmapper_index_params.SA_flag_iterater = 0;///unsigned int SA_flag_iterater = 0;

	bitmapper_index_params.sparse_suffix_array_length = 0;///unsigned int SA_number = 0;

	///这里SA要改
	bitmapper_bs_iter number_of_SA_flag_bits = 0;///unsigned int number_of_SA_flag_bits = 0;



	i = 0;

	SA_flag[bitmapper_index_params.SA_flag_iterater] = (SA_flag_string_type)0;


	///这里SA要改
	///注意这里SA_flag_iterater并不能++
	///SA_flag_iterater++;
	number_of_SA_flag_bits = number_of_SA_flag_bits + SA_counter_length;
	if (number_of_SA_flag_bits % bitmapper_index_params.SA_flag_warp_number == 0)
    {
        bitmapper_index_params.SA_flag_iterater++;
    }

	while (1)
	{

		tmp_SA_flag = (SA_flag_string_type)0;

		if (i >= bitmapper_index_params.SA_length)
		{
			break;
		}



		///if (uint64_t(sa[i]) % compress_sa == 0)  xxx
		if (uint64_t(get_sa_block(i)) % compress_sa == 0)
		{
			tmp_SA_flag = (SA_flag_string_type)1;
			bitmapper_index_params.sparse_suffix_array_length++;///SA_number++;

		}
		else
		{
			tmp_SA_flag = (SA_flag_string_type)0;
		}



		///这里SA要改
		///shift_length = SA_flag_warp_number - i % SA_flag_warp_number - 1;
		shift_length = bitmapper_index_params.SA_flag_warp_number -
            number_of_SA_flag_bits % bitmapper_index_params.SA_flag_warp_number - 1;

		tmp_SA_flag = tmp_SA_flag << shift_length;



		SA_flag[bitmapper_index_params.SA_flag_iterater] = SA_flag[bitmapper_index_params.SA_flag_iterater] | tmp_SA_flag;


		tmp_SA_flag = (SA_flag_string_type)0;



		///这里SA要改
		number_of_SA_flag_bits++;

		i++;

		if (number_of_SA_flag_bits%bitmapper_index_params.SA_flag_warp_number == 0)
		{
			bitmapper_index_params.SA_flag_iterater++;
		}




		if (i % bitmapper_index_params.compress_SA_flag == 0)
		{

			if (number_of_SA_flag_bits % bitmapper_index_params.SA_flag_warp_number != 0)
			{
				fprintf(stderr, "ERROR! \n");
			}



            SA_flag[bitmapper_index_params.SA_flag_iterater] = bitmapper_index_params.sparse_suffix_array_length;



            /**
            SA_flag[bitmapper_index_params.SA_flag_iterater]
                = (SA_flag_string_type)(SA_flag[bitmapper_index_params.SA_flag_iterater]
                                        << bitmapper_index_params.SA_counter_shift_length);
            **/







			number_of_SA_flag_bits = number_of_SA_flag_bits + SA_counter_length;



			if (number_of_SA_flag_bits % bitmapper_index_params.SA_flag_warp_number == 0)
            {
                bitmapper_index_params.SA_flag_iterater++;
            }

		}


		tmp_SA_flag = (SA_flag_string_type)0;


	}


	///这里SA要改
	/**
	if (i%SA_flag_warp_number != 0 &&
	i % compress_SA_flag != 0)
	{
	SA_flag_iterater++;
	}
	**/
	if (number_of_SA_flag_bits%bitmapper_index_params.SA_flag_warp_number != 0)
	{
		bitmapper_index_params.SA_flag_iterater++;
	}

	///这里SA要改,为了防止后面计算时溢出
	bitmapper_index_params.SA_flag_iterater++;

	fprintf(stderr, "SA_flag has been built\n");
	fflush(stderr);


	
	///这里要改
	fwrite(&bitmapper_index_params.high_occ_table_length, sizeof(bitmapper_index_params.high_occ_table_length), 1, fo);
	fwrite(high_occ_table, sizeof(high_occ_table_type), bitmapper_index_params.high_occ_table_length, fo);
	printf("Occ has been writed!\n");
	


	///fwrite(&bwt_iterater, sizeof(unsigned int), 1, fb);
	fwrite(&bitmapper_index_params.bwt_length, sizeof(bitmapper_bs_iter), 1, fb);///fwrite(&bwt_iterater, sizeof(bitmapper_bs_iter), 1, fb);
	fwrite(bwt, sizeof(bwt_string_type), bitmapper_index_params.bwt_length, fb);///fwrite(bwt, sizeof(bwt_string_type), bwt_iterater, fb);
	printf("BWT has been writed!\n");

	fwrite(&bitmapper_index_params.SA_length, sizeof(bitmapper_index_params.SA_length), 1, f2);
	fwrite(&bitmapper_index_params.shapline, sizeof(bitmapper_index_params.shapline), 1, f2);

	printf("cc&sharp_line has been writed!\n");



    nacgt[0] = 1;
    fwrite(&nacgt[0], sizeof(nacgt[0]), 1, f2);

    for (j = 1; j <= 4; j++)
    {
        nacgt[j] = nacgt[j] + nacgt[j - 1];
        fwrite(&nacgt[j], sizeof(nacgt[j]), 1, f2);
    }




    fwrite(&bitmapper_index_params.sparse_suffix_array_length, sizeof(bitmapper_index_params.sparse_suffix_array_length), 1, fs);///fwrite(&SA_number, sizeof(SA_number), 1, fs);





	///unsigned int ijkijkijkijk = 0;

	i = 0;

	unsigned int tmp_site = 0;






	while (1)
	{

		if (i >= bitmapper_index_params.SA_length)
		{
			break;
		}


		///if (uint64_t(sa[i]) % compress_sa == 0)  xxx
		if (uint64_t(get_sa_block(i)) % compress_sa == 0)
		{

			///if (compress_sa >= 4)
			///if (compress_sa >= 0)
			if (compress_sa >= 4)
			{
				///如果前一位是$，那么在sa必然是0，这个绝壁是不符合要求
				//所以不用记录，只是要判断是不是为0，为0就丢弃
				//不过这东西应该是被保存成010000000000000000000这个数
				//tmp_site = (unsigned int)0;
				////这个应该不需要吧
				////ch = refer[sa[i]] - 1;
				///if (uint64_t(sa[i]) != 0)  xxx
				if (uint64_t(get_sa_block(i)) != 0)
				{
					///ch = refer[uint64_t(sa[i]) - 1];  xxx
					ch = refer[uint64_t(get_sa_block(i)) - 1];
				}
				else
				{
					///没什么意义，仅仅为了和上一个版本保持一致...
					ch = 1;
				}



				tmp_site = (unsigned int)(ch);///tmp_site = (unsigned int)abs(ch);

				tmp_site = tmp_site << 30;

				///tmp_site = tmp_site | (uint64_t(sa[i]) / compress_sa);  xxx
				tmp_site = tmp_site | (uint64_t(get_sa_block(i)) / compress_sa);


				fwrite(&tmp_site, sizeof(tmp_site), 1, fs);
			}
			else
			{

				///fwrite(&sa[i], sizeof(sa[i]), 1, fs);  xxx
				tmp_site = uint64_t(get_sa_block(i));
				fwrite(&tmp_site, sizeof(tmp_site), 1, fs);

			}





			///ijkijkijkijk++;
		}

		i++;

	}

	fprintf(stderr, "sparse_suffix_array_length=%llu\n", bitmapper_index_params.sparse_suffix_array_length);

	fprintf(stderr, "SA_flag_iterater=%llu\n", bitmapper_index_params.SA_flag_iterater);


	fwrite(&bitmapper_index_params.SA_flag_iterater, sizeof(bitmapper_index_params.SA_flag_iterater), 1, fs);
	fwrite(SA_flag, sizeof(SA_flag_string_type), bitmapper_index_params.SA_flag_iterater, fs);





	fwrite(&compress_sa, sizeof(unsigned int), 1, f2);
	///fwrite(&compress_occ, sizeof(unsigned int), 1, f2);
	fwrite(&bitmapper_index_params.compress_occ, sizeof(unsigned int), 1, f2);

	///这里要改
	fwrite(&bitmapper_index_params.high_compress_occ, sizeof(unsigned int), 1, f2);



	bitmapper_bs_iter hash_table_i = 0;

	char pattern[17];

	char bit_to_car[3];


	/**
	ctoi['A'] = 2;
	ctoi['G'] = 0;
	ctoi['T'] = 1;
	ctoi['a'] = 2;
	ctoi['g'] = 0;
	ctoi['t'] = 1;
	**/

	bit_to_car[0] = 'G';
	bit_to_car[1] = 'T';
	bit_to_car[2] = 'A';



	

	/************************build 16-mer hash table******************************/


	bitmapper_index_params.hash_table_16_mer_size = (bitmapper_bs_iter)(pow(3.0, 16.0)) + 1;

	bitmapper_index_params.hash_table_16_mer_high_32 = (unsigned int*)malloc(sizeof(unsigned int)
		* bitmapper_index_params.hash_table_16_mer_size);

	bitmapper_index_params.hash_table_16_mer_low_8 = (unsigned char*)malloc(sizeof(unsigned char)
		* bitmapper_index_params.hash_table_16_mer_size);




	


	bitmapper_index_params.bwt = bwt;
	bitmapper_index_params.high_occ_table = high_occ_table;

    bitmapper_index_params.nacgt[0] = nacgt[0];
    bitmapper_index_params.nacgt[1] = nacgt[1];
    bitmapper_index_params.nacgt[2] = nacgt[2];
    bitmapper_index_params.nacgt[3] = nacgt[3];
    bitmapper_index_params.nacgt[4] = nacgt[4];

    bitmapper_index_params.ctoi['A'] = 2;
	bitmapper_index_params.ctoi['G'] = 0;
	bitmapper_index_params.ctoi['T'] = 1;
	bitmapper_index_params.ctoi['a'] = 2;
	bitmapper_index_params.ctoi['g'] = 0;
	bitmapper_index_params.ctoi['t'] = 1;

	





	bitmapper_bs_iter top, bot, pre_top, pre_bot;
	bitmapper_bs_iter mode_low8 = (bwt_string_type)255;
	unsigned int difference;
	///bitmapper_bs_iter tmp_bot=1, tmp_top;

	///hash表的起始应该是1,因为有个$
	bitmapper_index_params.hash_table_16_mer_high_32[0] = 0;
    bitmapper_index_params.hash_table_16_mer_low_8[0] = 1;



	


	for (hash_table_i = 0; hash_table_i < bitmapper_index_params.hash_table_16_mer_size - 1; hash_table_i++)
	{

		difference = 0;

		bitmapper_bs_iter convert_i = 0;
		bitmapper_bs_iter convert_tmp = hash_table_i;

		for (convert_i = 0; convert_i < 16; convert_i++)
		{
			pattern[15 - convert_i] = bit_to_car[convert_tmp % 3];
			convert_tmp = convert_tmp / 3;
		}

		pattern[16] = '\0';


		bitmapper_bs_iter num_occ = count(pattern, 16, &top, &bot, &pre_top, &pre_bot);

		/**
		bitmapper_bs_iter get_top = (bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i] << 8)
			| bitmapper_index_params.hash_table_16_mer_low_8[hash_table_i];
		**/
		/**
		bitmapper_bs_iter get_top = bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i];
		get_top = get_top << 8;
		get_top = get_top | bitmapper_index_params.hash_table_16_mer_low_8[hash_table_i];
		**/


		bitmapper_bs_iter tmp_high, tmp_low;
		///这一步是取低28位
		tmp_high = bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i] & bitmapper_index_params.hash_table_mode_get;
		tmp_high = tmp_high << 8;
		tmp_low = bitmapper_index_params.hash_table_16_mer_low_8[hash_table_i];
		tmp_high = tmp_high | tmp_low;
		bitmapper_bs_iter get_top = tmp_high;



        if(num_occ == 0)
        {
			top = get_top;
            bot = top;
        }
		else
		{
			///注意top有可能比get_top大.....
			///这其实就是那个hash表从上到下有可能不连续...
			///不连续的原因是因为存在xxxxx$这种字符串....
			if (top != get_top)
			{
				
				difference = top - get_top;
				difference = difference << bitmapper_index_params.hash_table_shift_length;

				/**
				fprintf(stderr, "****************************\n");
				fprintf(stderr, "hash_table_i=%llu\n****************************\n", hash_table_i);
				fprintf(stderr, "top=%llu, get_top=%llu\n", top, get_top);

				fprintf(stderr, "text_length=%llu, SA[get_top] =%llu \n", 
					text_length, uint64_t(get_sa_block(get_top)));

				
				

				fprintf(stderr, "****************************\n\n\n");
				**/

			}
		}

		/**
		if (hash_table_i == 7218558)
		{
			fprintf(stderr, "hash_table_i=%llu\n", hash_table_i);
			fprintf(stderr, "top=%llu, bot=%llu\n", top, bot);
			fprintf(stderr, "num_occ=%llu\n", num_occ);
		}
		**/

		bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i] = top >> 8;
		bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i]
			= bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i] | difference;
		bitmapper_index_params.hash_table_16_mer_low_8[hash_table_i] = top & mode_low8;

		bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i + 1] = bot >> 8;
		bitmapper_index_params.hash_table_16_mer_low_8[hash_table_i + 1] = bot & mode_low8;

				
	}



	/*************************************debug*******************************************/


	/**
	for (hash_table_i = 0; hash_table_i < bitmapper_index_params.hash_table_16_mer_size - 1; hash_table_i++)
	{

		bitmapper_bs_iter convert_i = 0;
		bitmapper_bs_iter convert_tmp = hash_table_i;

		for (convert_i = 0; convert_i < 16; convert_i++)
		{
			pattern[15 - convert_i] = bit_to_car[convert_tmp % 3];
			convert_tmp = convert_tmp / 3;
		}

		pattern[16] = '\0';



		bitmapper_bs_iter num_occ = count(pattern, 16, &top, &bot, &pre_top, &pre_bot);

		bitmapper_bs_iter get_top, get_bot;


		query_16_mer_hash_table(hash_table_i, &get_top, &get_bot);

		if(num_occ==0 && get_bot - get_top ==0)
        {
            continue;
        }

        if (get_top != top || get_bot != bot)
		{
			
			fprintf(stderr, "################################\n", hash_table_i);
			fprintf(stderr, "hash_table_i=%llu\n", hash_table_i);
			fprintf(stderr, "top=%llu, bot=%llu\n", top, bot);
			fprintf(stderr, "get_top=%llu, get_bot=%llu\n", get_top, get_bot);
			

		}

		bitmapper_bs_iter hash_value = get_3_letter_hash_value(pattern, 16);


		if (hash_value != hash_table_i)
		{
			fprintf(stderr, "hash error\n");
			fprintf(stderr, "hash_table_i = %llu\n", hash_table_i);

		}


	}
	**/







	

	/*************************************debug*******************************************/



	fwrite(&bitmapper_index_params.hash_table_16_mer_size, sizeof(bitmapper_bs_iter), 1, fb);
	fwrite(bitmapper_index_params.hash_table_16_mer_high_32,
		sizeof(unsigned int), bitmapper_index_params.hash_table_16_mer_size, fb);
	fwrite(bitmapper_index_params.hash_table_16_mer_low_8,
		sizeof(unsigned char), bitmapper_index_params.hash_table_16_mer_size, fb);



	printf("hash table has been writed!\n");



	/************************build 16-mer hash table******************************/






	fflush(fo);
	fflush(f2);
	fflush(fs);
	fflush(fb);








	/*******************************debug block*******************************************/

	/**

	bitmapper_index_params.bwt = bwt;
	bitmapper_index_params.high_occ_table = high_occ_table;


	bitmapper_index_params.nacgt[0] = nacgt[0];
	bitmapper_index_params.nacgt[1] = nacgt[1];
	bitmapper_index_params.nacgt[2] = nacgt[2];
	bitmapper_index_params.nacgt[3] = nacgt[3];
	bitmapper_index_params.nacgt[4] = nacgt[4];


	i = 0;

	nacgt[0] = 0;
	nacgt[1] = 0;
	nacgt[2] = 0;
	nacgt[3] = 0;
	nacgt[4] = 0;


	while (1)
	{


		if (i >= bitmapper_index_params.SA_length)
		{
			break;
		}


		if (uint64_t(get_sa_block(i)) != 0)
		{
			ch = refer[uint64_t(get_sa_block(i)) - 1];
		}
		else
		{
			bitmapper_index_params.shapline = i;

			i++;

			continue;
		}


		///nacgt[0]肯定是$
		nacgt[ch + 1]++;

		i++;






		bitmapper_bs_iter get_value;
		high_occ_table_type* tmp_high_occ_table = NULL;

		get_value =
			find_occ_fm_index(i, 1, bwt, tmp_high_occ_table) - bitmapper_index_params.nacgt[1];

		if (get_value != nacgt[2])
		{
			fprintf(stderr, "i=%llu, get_value=%llu, nacgt[2]=%llu, find_occ_fm_index T is error ...\n",
				i, get_value, nacgt[2]);
			return 1;
		}



		get_value =
			find_occ_fm_index(i, 2, bwt, tmp_high_occ_table) - bitmapper_index_params.nacgt[2];

		if (get_value != nacgt[3])
		{
			fprintf(stderr, "i=%llu, get_value=%llu, nacgt[3]=%llu, find_occ_fm_index G is error ...\n",
				i, get_value, nacgt[3]);
			return 1;
		}


		get_value =
			find_occ_fm_index(i, 0, bwt, tmp_high_occ_table) - bitmapper_index_params.nacgt[0];

		if (get_value != nacgt[1])
		{
			fprintf(stderr, "i=%llu, get_value=%llu, nacgt[1]=%llu, find_occ_fm_index A is error ...\n",
				i, get_value, nacgt[1]);
			return 1;
		}




	}




	fprintf(stderr, "bitmapper_index_params.nacgt[0] = %llu \n", bitmapper_index_params.nacgt[0]);
	fprintf(stderr, "bitmapper_index_params.nacgt[1] = %llu \n", bitmapper_index_params.nacgt[1]);
	fprintf(stderr, "bitmapper_index_params.nacgt[2] = %llu \n", bitmapper_index_params.nacgt[2]);
	fprintf(stderr, "bitmapper_index_params.nacgt[3] = %llu \n", bitmapper_index_params.nacgt[3]);
	fprintf(stderr, "bitmapper_index_params.nacgt[4] = %llu \n", bitmapper_index_params.nacgt[4]);

	fprintf(stderr, "nacgt[0] = %llu \n", nacgt[0]);
	fprintf(stderr, "nacgt[1] = %llu \n", nacgt[1]);
	fprintf(stderr, "nacgt[2] = %llu \n", nacgt[2]);
	fprintf(stderr, "nacgt[3] = %llu \n", nacgt[3]);
	fprintf(stderr, "nacgt[4] = %llu \n", nacgt[4]);


	**/






	/*******************************debug block*******************************************/










	///return 1;







	fprintf(stderr, "Sucess!\n");
	fflush(stderr);



	

	fclose(fo);
	fclose(f2);
	fclose(fs);
	fclose(fb);
































	





















































	return 0;


}



unsigned int init_bitmapper_index_params()
{
	bitmapper_bs_iter i;

	for (i = 0; i < 256; i++)
	{
		bitmapper_index_params.ctoi[i] = 3;
	}


	/**
	bitmapper_index_params.ctoi['A'] = 0;
	bitmapper_index_params.ctoi['G'] = 2;
	bitmapper_index_params.ctoi['T'] = 1;
	bitmapper_index_params.ctoi['a'] = 0;
	bitmapper_index_params.ctoi['g'] = 2;
	bitmapper_index_params.ctoi['t'] = 1;
	**/


	bitmapper_index_params.ctoi['A'] = 2;
	bitmapper_index_params.ctoi['G'] = 0;
	bitmapper_index_params.ctoi['T'] = 1;
	bitmapper_index_params.ctoi['a'] = 2;
	bitmapper_index_params.ctoi['g'] = 0;
	bitmapper_index_params.ctoi['t'] = 1;



	bitmapper_index_params.cut_thr = 1;

	bitmapper_index_params.tree_nodes_number = (pow(4, bitmapper_index_params.compress_sa) - 1) / 3;


	/**
	bitmapper_index_params.sp_tree = (unsigned int*)malloc(sizeof(unsigned int)* bitmapper_index_params.tree_nodes_number);
	bitmapper_index_params.ep_tree = (unsigned int*)malloc(sizeof(unsigned int)* bitmapper_index_params.tree_nodes_number);
	**/



	bitmapper_index_params.FMtree_queue
        = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)* bitmapper_index_params.tree_nodes_number*3);
		///= (unsigned int*)malloc(sizeof(unsigned int)* bitmapper_index_params.tree_nodes_number*3);



	bitmapper_index_params.FMtree_queue_start_point = 0;
	bitmapper_index_params.FMtree_queue_end_point = 0;





















	for (i = 0; i < 65536; i++)
	{
		bitmapper_index_params.hash_count[i] = __builtin_popcountll(i);
	}







}







unsigned int init_locate_queue_muti_thread(bwt_locate_queue* get_queue)
{

	get_queue->FMtree_queue
		= (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)* bitmapper_index_params.tree_nodes_number * 3);

	get_queue->FMtree_queue_start_point = 0;
	get_queue->FMtree_queue_end_point = 0;

}



unsigned int load_index(char* filename_prefix)
{
	char filename[1000], filename1[1000], filenames[1000], filenameo[1000], filenameb[1000];
	long long i, j, t;

	FILE *f1, *f2, *fs, *fb, *fout, *fo;
	///FILE *f1, *f2, *fs, *fb, *fout;



	bitmapper_index_params.pop_count_mode[0] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[0]
			= (bwt_string_type)(bitmapper_index_params.pop_count_mode[0] << 2) | (bwt_string_type)3;
	}


	bitmapper_index_params.pop_count_mode[1] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[1] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[1] << 2) | (bwt_string_type)2;
	}


	bitmapper_index_params.pop_count_mode[2] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[2] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[2] << 2) | (bwt_string_type)1;
	}

	bitmapper_index_params.pop_count_mode[3] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[3] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[3] << 2) | (bwt_string_type)0;
	}

	bitmapper_index_params.mode_high = (bwt_string_type)0;
	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.mode_high =
			(bwt_string_type)(bitmapper_index_params.mode_high << 2) | (bwt_string_type)2;
	}


	bitmapper_index_params.mode_low = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.mode_low =
			(bwt_string_type)(bitmapper_index_params.mode_low << 2) | (bwt_string_type)1;
	}




	strcpy(filename, filename_prefix);
	strcpy(filename + strlen(filename), ".index");
	strcpy(filenames, filename);
	strcpy(filenameo, filename);
	strcpy(filenameb, filename);
	strcpy(filenames + strlen(filename), ".sa");
	strcpy(filenameo + strlen(filename), ".occ");
	strcpy(filenameb + strlen(filename), ".bwt");


	///这里要改
	
	fo = fopen(filenameo, "r");
	if (fo == NULL)
	{
		fprintf(stderr, "Failed to open %s! FMtree will exit ...\n", filenameo);
		return 1;
	}
	

	fs = fopen(filenames, "r");
	if (fs == NULL)
	{
		fprintf(stderr, "Failed to open %s! FMtree will exit ...\n", filenames);
		return 1;
	}
	fb = fopen(filenameb, "r");
	if (fb == NULL)
	{
		fprintf(stderr, "Failed to open %s! FMtree will exit ...\n", filenameb);
		return 1;
	}
	strcpy(filename1, filename);
	f2 = fopen(filename1, "r");

	if (f2 == NULL)
	{
		fprintf(stderr, "Failed to open %s! FMtree will exit ...\n", filename1);
		return 1;
	}


	fread(&bitmapper_index_params.SA_length, sizeof(bitmapper_index_params.SA_length), 1, f2);
	fread(&bitmapper_index_params.shapline, sizeof(bitmapper_index_params.shapline), 1, f2);

	fprintf(stderr, "shapline=%llu\n", bitmapper_index_params.shapline);



	for (j = 0; j <= 4; j++)
		fread(&bitmapper_index_params.nacgt[j], sizeof(bitmapper_index_params.nacgt[j]), 1, f2);




	fread(&bitmapper_index_params.compress_sa, sizeof(bitmapper_index_params.compress_sa), 1, f2);
	fread(&bitmapper_index_params.compress_occ, sizeof(bitmapper_index_params.compress_occ), 1, f2);
	fread(&bitmapper_index_params.high_compress_occ, sizeof(bitmapper_index_params.high_compress_occ), 1, f2);


	bitmapper_index_params.mode_4[0] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[1] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[2] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[3] = (bwt_string_type)0;


	fread(&bitmapper_index_params.bwt_length, sizeof(bitmapper_index_params.bwt_length), 1, fb);
	bitmapper_index_params.bwt = (bwt_string_type *)malloc(sizeof(bwt_string_type)*bitmapper_index_params.bwt_length);
	fread(bitmapper_index_params.bwt, sizeof(bwt_string_type), bitmapper_index_params.bwt_length, fb);
	fprintf(stderr,"BWT has been loaded!\n");

	fread(&bitmapper_index_params.hash_table_16_mer_size, sizeof(bitmapper_bs_iter), 1, fb);


	bitmapper_index_params.hash_table_16_mer_high_32 = (unsigned int*)malloc(sizeof(unsigned int)
		* bitmapper_index_params.hash_table_16_mer_size);

	bitmapper_index_params.hash_table_16_mer_low_8 = (unsigned char*)malloc(sizeof(unsigned char)
		* bitmapper_index_params.hash_table_16_mer_size);


	fread(bitmapper_index_params.hash_table_16_mer_high_32,
		sizeof(unsigned int), bitmapper_index_params.hash_table_16_mer_size, fb);
	fread(bitmapper_index_params.hash_table_16_mer_low_8,
		sizeof(unsigned char), bitmapper_index_params.hash_table_16_mer_size, fb);

	fprintf(stderr,"hash table has been loaded!\n");



	fprintf(stderr,"SA_length=%u\n", bitmapper_index_params.SA_length);
	fread(&bitmapper_index_params.sparse_suffix_array_length,
		sizeof(bitmapper_index_params.sparse_suffix_array_length), 1, fs);
	fprintf(stderr,"sparse_suffix_array_length=%u\n", bitmapper_index_params.sparse_suffix_array_length);
	bitmapper_index_params.sa
		= (unsigned int *)malloc(sizeof(unsigned int)*(bitmapper_index_params.sparse_suffix_array_length));
	fread(bitmapper_index_params.sa, sizeof(unsigned int), bitmapper_index_params.sparse_suffix_array_length, fs);








	fread(&bitmapper_index_params.SA_flag_iterater,
		sizeof(bitmapper_index_params.SA_flag_iterater), 1, fs);

	fprintf(stderr, "SA_flag_iterater=%llu\n", bitmapper_index_params.SA_flag_iterater);

	bitmapper_index_params.SA_flag =
		(SA_flag_string_type*)malloc(sizeof(SA_flag_string_type)*bitmapper_index_params.SA_flag_iterater);

	fread(bitmapper_index_params.SA_flag, sizeof(SA_flag_string_type),
		bitmapper_index_params.SA_flag_iterater, fs);

	
	fread(&bitmapper_index_params.high_occ_table_length, sizeof(bitmapper_index_params.high_occ_table_length), 1, fo);
	fprintf(stderr, "high_occ_table_length=%llu\n", bitmapper_index_params.high_occ_table_length);
	bitmapper_index_params.high_occ_table = (high_occ_table_type*)malloc(sizeof(high_occ_table_type)
		*bitmapper_index_params.high_occ_table_length);
	fread(bitmapper_index_params.high_occ_table, sizeof(high_occ_table_type),
		bitmapper_index_params.high_occ_table_length, fo);
	








	fclose(f2);
	fclose(fs);
	fclose(fb);
	fclose(fo);

	init_bitmapper_index_params();

























	/*************************************debug*******************************************/

	/**
	bitmapper_bs_iter hash_table_i = 0;

	char pattern[17];

	char bit_to_car[3];

	bit_to_car[0] = 'G';
	bit_to_car[1] = 'T';
	bit_to_car[2] = 'A';

	bitmapper_bs_iter top, bot, pre_top, pre_bot;

	for (hash_table_i = 0; hash_table_i < bitmapper_index_params.hash_table_16_mer_size - 1; hash_table_i++)
	{

		bitmapper_bs_iter convert_i = 0;
		bitmapper_bs_iter convert_tmp = hash_table_i;

		for (convert_i = 0; convert_i < 16; convert_i++)
		{
			pattern[15 - convert_i] = bit_to_car[convert_tmp % 3];
			convert_tmp = convert_tmp / 3;
		}

		pattern[16] = '\0';



		bitmapper_bs_iter num_occ = count(pattern, 16, &top, &bot, &pre_top, &pre_bot);

		bitmapper_bs_iter get_top, get_bot;


		query_16_mer_hash_table(hash_table_i, &get_top, &get_bot);

		if (num_occ == 0 && get_bot - get_top == 0)
		{
			continue;
		}

		if (get_top != top || get_bot != bot)
		{

			fprintf(stderr, "################################\n", hash_table_i);
			fprintf(stderr, "hash_table_i=%llu\n", hash_table_i);
			fprintf(stderr, "top=%llu, bot=%llu\n", top, bot);
			fprintf(stderr, "get_top=%llu, get_bot=%llu\n", get_top, get_bot);

			fprintf(stderr, "hash_table_i=%llu, high_32=%llu, low8=%llu\n", 
				hash_table_i, 
				bitmapper_index_params.hash_table_16_mer_high_32[hash_table_i], 
				bitmapper_index_params.hash_table_16_mer_low_8[hash_table_i]);




		}

		bitmapper_bs_iter hash_value = get_3_letter_hash_value(pattern, 16);


		if (hash_value != hash_table_i)
		{
			fprintf(stderr, "hash error\n");
			fprintf(stderr, "hash_table_i = %llu\n", hash_table_i);

		}


	}


	fprintf(stderr, "Debug sucess!\n");
	**/
	

	/*************************************debug*******************************************/
















	return 1;



}


























inline unsigned int bwt_get_sa_restrict_zero_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	unsigned int last_wods;

	SA_flag_string_type tmp_SA_pop_count;

	l = sp;

	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
		*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;


	///这里SA要改
	///last = l % compress_SA_flag;
	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;




	///这里SA要改
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}


	///(*start) = sa[ans];
	(*start) = (sa[ans] & bitmapper_index_params.SA_header_mode)*bitmapper_index_params.compress_sa;


	return 1;

}




inline unsigned int bwt_get_sa_restrict_zero_steps
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	unsigned int last_wods;

	SA_flag_string_type tmp_SA_pop_count;


	l = sp;

	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
		*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;



	///这里SA要改
	///last = l % compress_SA_flag;
	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;


	///这里SA要改
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}


	(*start) = sa[ans];

	return 1;

}


inline bitmapper_bs_iter bwt_accesss_pre_SA(
    unsigned int *sa, SA_flag_string_type* SA_flag,
    bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter* tmp_SA_length,
    bitmapper_bs_iter delta)
{
	bitmapper_bs_iter SA_start, SA_length, t, tmp_SA, tmp_SA_tail;

	///if (ep - sp>1)
	{
		///这里SA要改
		bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);

		/**
		if (SA_length>150000)
		{
			return 0;
		}
		**/


		for (t = 0; t<SA_length; t++)
		{
			tmp_SA = sa[t + SA_start];
			tmp_SA_tail = (tmp_SA & bitmapper_index_params.SA_header_mode)
				*bitmapper_index_params.compress_sa;



			if (tmp_SA_tail != 0 &&
				((bitmapper_index_params.SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)++] = tmp_SA_tail - 1;
			}


			///locates[(*tmp_SA_length)++] = sa[t + SA_start] & SA_header_mode - 1;

		}
	}/**
	else if (ep - sp == 1)   ///说实话程序是不可能到这里的
	{




		///这里SA要改
		if (bwt_get_sa_restrict_zero_steps(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{


			tmp_SA = locates[(*tmp_SA_length)];
			tmp_SA_tail = (tmp_SA & bitmapper_index_params.SA_header_mode)*bitmapper_index_params.compress_sa;



			if (tmp_SA_tail != 0 &&
				((bitmapper_index_params.SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)] = tmp_SA_tail - 1;
				(*tmp_SA_length)++;
			}



			///locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] & SA_header_mode - 1;
			///(*tmp_SA_length)++;

		}
	}
	**/

	return 1;

}



inline void bwt_accesss_SA_more_than_3(
    unsigned int *sa, SA_flag_string_type* SA_flag,
    bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter occ,
	bitmapper_bs_iter* tmp_SA_length)
{
	bitmapper_bs_iter SA_start, SA_length, t;

	///这里SA要改
	bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);

	for (t = 0; t<SA_length; t++)
	{
		locates[(*tmp_SA_length)++] = (bitmapper_bs_iter)(sa[t + SA_start] & bitmapper_index_params.SA_header_mode)
			*bitmapper_index_params.compress_sa + occ;
	}


}

/**
inline void bwt_accesss_SA_less_than_4(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length)
{
	unsigned int SA_start, SA_length, t;

	///这里SA要改
	bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);

	for (t = 0; t<SA_length; t++)
	{
		locates[(*tmp_SA_length)++] = sa[t + SA_start] + occ;
	}


}
**/

/**
inline void bwt_accesss_SA_more_than_3_back_up(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length)
{
	unsigned int SA_start, SA_length, t;

	if (ep - sp>1)
	{


		///这里SA要改
		bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);




		for (t = 0; t<SA_length; t++)
		{
			locates[(*tmp_SA_length)++] = (sa[t + SA_start] & bitmapper_index_params.SA_header_mode)
				*bitmapper_index_params.compress_sa + occ;
		}
	}
	else if (ep - sp == 1)
	{
		///这里SA要改
		if (bwt_get_sa_restrict_zero_steps_more_than_3(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;

		}

	}

}
**/


inline void bwt_find_occ_all_sp_ep_optimal_back(
    bitmapper_bs_iter sp,
    bitmapper_bs_iter ep,
    bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	bitmapper_bs_iter* ans_A_sp,
	bitmapper_bs_iter* ans_C_sp,
	bitmapper_bs_iter* ans_G_sp,
	bitmapper_bs_iter* ans_T_sp,
	bitmapper_bs_iter* ans_A_ep,
	bitmapper_bs_iter* ans_C_ep,
	bitmapper_bs_iter* ans_G_ep,
	bitmapper_bs_iter* ans_T_ep)
{


	if (sp > bitmapper_index_params.shapline)
		sp--;
	if (ep > bitmapper_index_params.shapline)
		ep--;



	bitmapper_bs_iter i_sp;
	bitmapper_bs_iter i_ep;

	bitmapper_bs_iter flag = 0;

	bwt_string_type P_tmp, P_A, P_B;

	bitmapper_bs_iter tmp_ans_C_sp, tmp_ans_T_sp;
	bitmapper_bs_iter tmp_ans_C_ep, tmp_ans_T_ep;
	bitmapper_bs_iter rank_1_sp, rank_1_ep;

	rank_1_sp = 0;
	rank_1_ep = 0;

	tmp_ans_C_sp = 0;
	tmp_ans_T_sp = 0;

	tmp_ans_C_ep = 0;
	tmp_ans_T_ep = 0;


	bitmapper_bs_iter high_occ_table_line = (sp >> 16) << 2;


	(*ans_A_sp) = high_occ_table[high_occ_table_line];
	(*ans_C_sp) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_sp) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_sp) = high_occ_table[high_occ_table_line + 3];




	high_occ_table_line = (ep >> 16) << 2;
	(*ans_A_ep) = high_occ_table[high_occ_table_line];
	(*ans_C_ep) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_ep) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_ep) = high_occ_table[high_occ_table_line + 3];
	///从概率上来讲，ep和sp的superblock在一起的概率还是很大的,所以一开始一起查cache命中率还是不低的



	bitmapper_bs_iter actually_line_sp = ((sp >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;


	bitmapper_bs_iter actually_line_ep = ((ep >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;


	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}



	///sp的block
	(*ans_A_sp) += (bwt[actually_line_sp] >> 48);

	(*ans_C_sp) += (bwt[actually_line_sp] >> 32) & bitmapper_index_params.mode_16;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 16) & bitmapper_index_params.mode_16;

	(*ans_T_sp) += bwt[actually_line_sp] & bitmapper_index_params.mode_16;

	actually_line_sp++;


	bitmapper_bs_iter need_line_sp = sp & 255;
	bitmapper_bs_iter need_line_ep = ep & 255;


	if (need_line_sp != 0)
	{

		i_sp = 0;

		while (i_sp + bitmapper_index_params.bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);

			///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);



			i_sp = i_sp + bitmapper_index_params.bwt_warp_number;

		}


		///这个判断条件貌似是sp和ep都属于同一个block
		///所以ep的值貌似可以直接继承
		if (flag == 1)
		{
			i_ep = i_sp;

			rank_1_ep = rank_1_sp;

			tmp_ans_C_ep = tmp_ans_C_sp;
		    (*ans_C_ep) = (*ans_C_sp);

			(*ans_G_ep) = (*ans_G_sp);

			tmp_ans_T_ep = tmp_ans_T_sp;
			(*ans_T_ep) = (*ans_T_sp);


			actually_line_ep = actually_line_sp;

			flag = 2;
		}


		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];

			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) << 1);


			P_tmp = P_tmp & P_B;

			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);


			///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);


		}



		(*ans_C_sp) = (*ans_C_sp) + tmp_ans_C_sp;
		(*ans_G_sp) = (*ans_G_sp) + rank_1_sp - tmp_ans_C_sp - (tmp_ans_T_sp << 1);
		(*ans_T_sp) = (*ans_T_sp) + tmp_ans_T_sp;

		/**
		if ((bitmapper_index_params.shapline >=
			((sp >> 8) << 8))
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;


		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>bitmapper_index_params.shapline)
		{
			(*ans_A_sp)--;
		}
		**/


		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


	}



	(*ans_A_sp) += bitmapper_index_params.nacgt[0];
	(*ans_C_sp) += bitmapper_index_params.nacgt[1];
	(*ans_G_sp) += bitmapper_index_params.nacgt[2];
	(*ans_T_sp) += bitmapper_index_params.nacgt[3];


















	if (need_line_ep != 0)
	{

		///如果进到这个循环，则说明sp和ep属于一个block
		///那么直接接着count就好, 这个cache命中率还是蛮高的
		if (flag == 2)
		{



			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];

				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);

				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);


				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];

				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);

				P_tmp = P_tmp & P_B;


				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);


				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);




			}


			(*ans_C_ep) = (*ans_C_ep) + tmp_ans_C_ep;
			(*ans_G_ep) = (*ans_G_ep) + rank_1_ep - tmp_ans_C_ep - (tmp_ans_T_ep << 1);
			(*ans_T_ep) = (*ans_T_ep) + tmp_ans_T_ep;


			/**
			if ((bitmapper_index_params.shapline >=
				((ep >> 8) << 8))
				&& (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}
			**/

			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);









		}
		else
		{

			///这里就是sp和ep在不同的block，这样的话这就是一次cache损失
			(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

			(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

			(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

			(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

			actually_line_ep++;



			i_ep = 0;

			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];

				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);


				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);


				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];

				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);



				P_tmp = P_tmp & P_B;

				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);


				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);

			}


			(*ans_C_ep) = (*ans_C_ep) + tmp_ans_C_ep;
			(*ans_G_ep) = (*ans_G_ep) + rank_1_ep - tmp_ans_C_ep - (tmp_ans_T_ep<<1);
			(*ans_T_ep) = (*ans_T_ep) + tmp_ans_T_ep;


			/**
			if ((bitmapper_index_params.shapline >= ((ep >> 8) << 8)) && (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}

			**/

			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);

		}




	}
	else
	{
		(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

		(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

		(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

		(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

		actually_line_ep++;
	}



	(*ans_A_ep) += bitmapper_index_params.nacgt[0];
	(*ans_C_ep) += bitmapper_index_params.nacgt[1];
	(*ans_G_ep) += bitmapper_index_params.nacgt[2];
	(*ans_T_ep) += bitmapper_index_params.nacgt[3];








	return;
}


inline void bwt_find_occ_all_sp_ep_optimal(
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	bitmapper_bs_iter* ans_0_sp,
	bitmapper_bs_iter* ans_1_sp,
	bitmapper_bs_iter* ans_2_sp,
	bitmapper_bs_iter* ans_0_ep,
	bitmapper_bs_iter* ans_1_ep,
	bitmapper_bs_iter* ans_2_ep)
{

	/*************************************calculate sp******************************************************/

	if (sp > bitmapper_index_params.shapline)
		sp--;



	bitmapper_bs_iter actually_line_sp;
	bitmapper_bs_iter high_occ_actually_line_sp;
	bitmapper_bs_iter occ_value_sp_1;
	bitmapper_bs_iter occ_value_sp_2;
	bitmapper_bs_iter high_occ_value_sp_1;
	bitmapper_bs_iter high_occ_value_sp_2;


	get_occ_value_sp_all(sp, &actually_line_sp, &high_occ_actually_line_sp, &occ_value_sp_1, &occ_value_sp_2,
		&high_occ_value_sp_1, &high_occ_value_sp_2);

	(*ans_1_sp) = occ_value_sp_1;
	(*ans_2_sp) = occ_value_sp_2;




	bitmapper_bs_iter need_sp = (sp & 63);

	bwt_string_type P_A;

	if (likely(need_sp != 0))
	{

		/***************************************delta = 1***************************************************/
		P_A = (bwt[actually_line_sp] >> (64 - need_sp));
		(*ans_1_sp) += __builtin_popcountll(P_A);
		/***************************************delta = 1***************************************************/



		/***************************************delta = 2***************************************************/
		P_A = (bwt[actually_line_sp + 1] >> (64 - need_sp));
		(*ans_2_sp) += __builtin_popcountll(P_A);
		/***************************************delta = 2***************************************************/


	}

	/***************************************delta = 0***************************************************/
	(*ans_0_sp) = sp - (*ans_1_sp) - (*ans_2_sp);
	/***************************************delta = 0***************************************************/



	(*ans_0_sp) += bitmapper_index_params.nacgt[0];
	(*ans_1_sp) += bitmapper_index_params.nacgt[1];
	(*ans_2_sp) += bitmapper_index_params.nacgt[2];



	/*************************************calculate sp******************************************************/






	/*************************************calculate ep******************************************************/


	if (ep > bitmapper_index_params.shapline)
		ep--;

	bitmapper_bs_iter actually_line_ep;
	bitmapper_bs_iter occ_value_ep_1;
	bitmapper_bs_iter occ_value_ep_2;

	get_occ_value_ep_all(ep, actually_line_sp, high_occ_actually_line_sp, occ_value_sp_1, occ_value_sp_2,
		high_occ_value_sp_1, high_occ_value_sp_2, &actually_line_ep, &occ_value_ep_1, &occ_value_ep_2);


	(*ans_1_ep) = occ_value_ep_1;
	(*ans_2_ep) = occ_value_ep_2;

	bitmapper_bs_iter need_ep = (ep & 63);

	if (likely(need_ep != 0))
	{
		/***************************************delta = 1***************************************************/
		P_A = (bwt[actually_line_ep] >> (64 - need_ep));
		(*ans_1_ep) += __builtin_popcountll(P_A);
		/***************************************delta = 1***************************************************/

		/***************************************delta = 2***************************************************/
		P_A = (bwt[actually_line_ep + 1] >> (64 - need_ep));
		(*ans_2_ep) += __builtin_popcountll(P_A);
		/***************************************delta = 2***************************************************/
	}

	/***************************************delta = 0***************************************************/
	(*ans_0_ep) = ep - (*ans_1_ep) - (*ans_2_ep);
	/***************************************delta = 0***************************************************/



	(*ans_0_ep) += bitmapper_index_params.nacgt[0];
	(*ans_1_ep) += bitmapper_index_params.nacgt[1];
	(*ans_2_ep) += bitmapper_index_params.nacgt[2];

}





inline void bwt_find_occ_all_sp_ep_optimal_small(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep)
{


	unsigned int i_sp;

	unsigned int flag = 0;

	bwt_string_type P_tmp, P_A, P_B;

	unsigned int tmp_ans_C_sp, tmp_ans_T_sp;
	unsigned int rank_1_sp = 0;


	tmp_ans_C_sp = 0;
	tmp_ans_T_sp = 0;


	unsigned int high_occ_table_line = (sp >> 16) << 2;


	(*ans_A_sp) = high_occ_table[high_occ_table_line];
	(*ans_C_sp) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_sp) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_sp) = high_occ_table[high_occ_table_line + 3];


	unsigned int actually_line_sp = ((sp >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;



	///sp的block
	(*ans_A_sp) += (bwt[actually_line_sp] >> 48);

	(*ans_C_sp) += (bwt[actually_line_sp] >> 32) & bitmapper_index_params.mode_16;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 16) & bitmapper_index_params.mode_16;

	(*ans_T_sp) += bwt[actually_line_sp] & bitmapper_index_params.mode_16;

	actually_line_sp++;


	unsigned int need_line_sp = sp & 255;


	if (need_line_sp != 0)
	{

		i_sp = 0;

		while (i_sp + bitmapper_index_params.bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);

			///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);



			i_sp = i_sp + bitmapper_index_params.bwt_warp_number;

		}

		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];

			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) << 1);


			P_tmp = P_tmp & P_B;

			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);


			///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);


		}



		(*ans_C_sp) = (*ans_C_sp) + tmp_ans_C_sp;
		(*ans_G_sp) = (*ans_G_sp) + rank_1_sp - tmp_ans_C_sp - (tmp_ans_T_sp << 1);
		(*ans_T_sp) = (*ans_T_sp) + tmp_ans_T_sp;


		if ((bitmapper_index_params.shapline >=
			((sp >> 8) << 8))
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;


		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>bitmapper_index_params.shapline)
		{
			(*ans_A_sp)--;
		}


	}



	(*ans_A_sp) += bitmapper_index_params.nacgt[0];
	(*ans_C_sp) += bitmapper_index_params.nacgt[1];
	(*ans_G_sp) += bitmapper_index_params.nacgt[2];
	(*ans_T_sp) += bitmapper_index_params.nacgt[3];







	unsigned int delta;

	unsigned char count_ep[4];
	count_ep[0] = 0;
	count_ep[1] = 0;
	count_ep[2] = 0;
	count_ep[3] = 0;



	while (sp<ep)
	{
		actually_line_sp = ((sp >> 8)*bitmapper_index_params.acctuall_bwt_gap) + (sp & 255) + 32;


		delta = (bwt[(actually_line_sp >> 5)]
			>> ((bitmapper_index_params.bwt_warp_number - (actually_line_sp & 31) - 1) << 1))
			&(bwt_string_type)3;

		count_ep[delta]++;

		sp++;
	}


	(*ans_A_ep) = (*ans_A_sp) + count_ep[0];
	(*ans_C_ep) = (*ans_C_sp) + count_ep[1];
	(*ans_G_ep) = (*ans_G_sp) + count_ep[2];
	(*ans_T_ep) = (*ans_T_sp) + count_ep[3];

	return;
}




inline void bwt_find_occ_all_sp_ep_optimal_back_up(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep)
{


	unsigned int i_sp;
	unsigned int i_ep;

	unsigned int flag = 0;

	bwt_string_type P_tmp, P_A, P_B;



	///unsigned int high_occ_table_line = (sp / bitmapper_index_params.high_compress_occ) * 4;
	///unsigned int high_occ_table_line = (sp / bitmapper_index_params.high_compress_occ) << 2;
	unsigned int high_occ_table_line = (sp >> 16) << 2;


	(*ans_A_sp) = high_occ_table[high_occ_table_line];
	(*ans_C_sp) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_sp) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_sp) = high_occ_table[high_occ_table_line + 3];



	///high_occ_table_line = (ep / bitmapper_index_params.high_compress_occ) * 4;
	///high_occ_table_line = (ep / bitmapper_index_params.high_compress_occ) << 2;
	high_occ_table_line = (ep >> 16) << 2;
	(*ans_A_ep) = high_occ_table[high_occ_table_line];
	(*ans_C_ep) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_ep) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_ep) = high_occ_table[high_occ_table_line + 3];
	///从概率上来讲，ep和sp的superblock在一起的概率还是很大的,所以一开始一起查cache命中率还是不低的



	/**
	unsigned int actually_line_sp = ((sp / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.acctuall_bwt_gap) / bitmapper_index_params.bwt_warp_number;
	**/

	unsigned int actually_line_sp = ((sp >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;



	/**
	unsigned int actually_line_ep = ((ep / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.acctuall_bwt_gap) / bitmapper_index_params.bwt_warp_number;
	**/

	unsigned int actually_line_ep = ((ep >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;










	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}






















	///sp的block
	(*ans_A_sp) += (bwt[actually_line_sp] >> 48);

	(*ans_C_sp) += (bwt[actually_line_sp] >> 32) & bitmapper_index_params.mode_16;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 16) & bitmapper_index_params.mode_16;

	(*ans_T_sp) += bwt[actually_line_sp] & bitmapper_index_params.mode_16;

	actually_line_sp++;

	/**
	unsigned int need_line_sp = sp % bitmapper_index_params.compress_occ;
	unsigned int need_line_ep = ep % bitmapper_index_params.compress_occ;
	**/
	unsigned int need_line_sp = sp & 255;
	unsigned int need_line_ep = ep & 255;


	if (need_line_sp != 0)
	{




		i_sp = 0;

		while (i_sp + bitmapper_index_params.bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);



			/**
			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);
			**/
			///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);





			i_sp = i_sp + bitmapper_index_params.bwt_warp_number;

		}


		///这个判断条件貌似是sp和ep都属于同一个block
		///所以ep的值貌似可以直接继承
		if (flag == 1)
		{
			i_ep = i_sp;


			(*ans_C_ep) = (*ans_C_sp);

			(*ans_G_ep) = (*ans_G_sp);

			(*ans_T_ep) = (*ans_T_sp);

			actually_line_ep = actually_line_sp;

			flag = 2;
		}


		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];
			/**
			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) * 2);
			**/

			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) << 1);


			P_tmp = P_tmp & P_B;



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);

			/**
			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);
			**/
			///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);






		}




		///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);



		/**
		if ((bitmapper_index_params.shapline >=
			(sp / bitmapper_index_params.compress_occ)*bitmapper_index_params.compress_occ)
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;
		**/


		if ((bitmapper_index_params.shapline >=
			((sp >> 8) << 8))
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;


		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>bitmapper_index_params.shapline)
		{
			(*ans_A_sp)--;
		}


	}



	(*ans_A_sp) += bitmapper_index_params.nacgt[0];
	(*ans_C_sp) += bitmapper_index_params.nacgt[1];
	(*ans_G_sp) += bitmapper_index_params.nacgt[2];
	(*ans_T_sp) += bitmapper_index_params.nacgt[3];


















	if (need_line_ep != 0)
	{

		///如果进到这个循环，则说明sp和ep属于一个block
		///那么直接接着count就好, 这个cache命中率还是蛮高的
		if (flag == 2)
		{



			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);






				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];

				///P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) * 2);

				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);

				P_tmp = P_tmp & P_B;



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);





			}




			///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);



			/**
			if ((bitmapper_index_params.shapline >=
				(ep / bitmapper_index_params.compress_occ)*bitmapper_index_params.compress_occ)
				&& (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;
			**/
			if ((bitmapper_index_params.shapline >=
				((ep >> 8) << 8))
				&& (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}












		}
		else
		{

			///这里就是sp和ep在不同的block，这样的话这就是一次cache损失
			(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

			(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

			(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

			(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

			actually_line_ep++;



			i_ep = 0;

			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);

				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);






				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];


				///P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) * 2);
				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);



				P_tmp = P_tmp & P_B;



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///因为任何数亦或0等于没进行任何操作，所以这里和上面有点不一样，不需要亦或
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);





			}


			/**
			if ((bitmapper_index_params.shapline >= (ep / bitmapper_index_params.compress_occ)
				*bitmapper_index_params.compress_occ) && (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;
			**/
			if ((bitmapper_index_params.shapline >= ((ep >> 8)<<8)) && (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}

		}




	}
	else
	{
		(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

		(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

		(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

		(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

		actually_line_ep++;
	}



	(*ans_A_ep) += bitmapper_index_params.nacgt[0];
	(*ans_C_ep) += bitmapper_index_params.nacgt[1];
	(*ans_G_ep) += bitmapper_index_params.nacgt[2];
	(*ans_T_ep) += bitmapper_index_params.nacgt[3];








	return;
}



inline bitmapper_bs_iter enqueue_FMtree(bitmapper_bs_iter sp, bitmapper_bs_iter ep, bitmapper_bs_iter layer)
{
	bitmapper_bs_iter queue_point = bitmapper_index_params.FMtree_queue_end_point * 3;
	bitmapper_index_params.FMtree_queue[queue_point++] = sp;
	bitmapper_index_params.FMtree_queue[queue_point++] = ep;
	bitmapper_index_params.FMtree_queue[queue_point] = layer;


	bitmapper_index_params.FMtree_queue_end_point++;
}

inline bitmapper_bs_iter dequeue_FMtree(bitmapper_bs_iter* sp, bitmapper_bs_iter* ep, bitmapper_bs_iter* layer)

{

	bitmapper_bs_iter queue_point = bitmapper_index_params.FMtree_queue_start_point * 3;
	*sp =
		bitmapper_index_params.FMtree_queue[queue_point++];
	*ep =
		bitmapper_index_params.FMtree_queue[queue_point++];


	*layer =
		bitmapper_index_params.FMtree_queue[queue_point];


	bitmapper_index_params.FMtree_queue_start_point++;
}

inline bitmapper_bs_iter queue_length_FMtree()
{
	return (bitmapper_index_params.FMtree_queue_end_point - bitmapper_index_params.FMtree_queue_start_point);
}

inline bitmapper_bs_iter empty_queue_FMtree()
{
	bitmapper_index_params.FMtree_queue_start_point = 0;
	bitmapper_index_params.FMtree_queue_end_point = 0;
}







bitmapper_bs_iter locate(char* pattern,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter sp1,
	bitmapper_bs_iter ep1,
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter length_read,
	bitmapper_bs_iter* occurrences)
{

	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter number_of_hits = ep - sp;

	bitmapper_bs_iter current_sp;
	bitmapper_bs_iter current_ep;
	bitmapper_bs_iter current_layer;

	/**
	bitmapper_bs_iter sp_A;
	bitmapper_bs_iter sp_C;
	bitmapper_bs_iter sp_G;
	bitmapper_bs_iter sp_T;


	bitmapper_bs_iter ep_A;
	bitmapper_bs_iter ep_C;
	bitmapper_bs_iter ep_G;
	bitmapper_bs_iter ep_T;
	**/
	bitmapper_bs_iter sp_0;
	bitmapper_bs_iter sp_1;
	bitmapper_bs_iter sp_2;


	bitmapper_bs_iter ep_0;
	bitmapper_bs_iter ep_1;
	bitmapper_bs_iter ep_2;


	///因为用了early leaf node calculation
	bitmapper_bs_iter tree_height = bitmapper_index_params.compress_sa - 1;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;



	///如果ep-sp在阈值之下，则所有位置one-to-one的拿到，这样就可以直接返回了
	if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
	{

		///注意这里最多是向下扩展bitmapper_index_params.compress_sa - current_layer - 1
		///也就是D-1步
		///和后面的循环里的结果不一样
		///后面的循环里因为已经经过了early leaf node calculation, 所以应该是
		///tree_height - current_layer - 1
		bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length,
			bitmapper_index_params.compress_sa - current_layer - 1);

		(*occurrences) = tmp_SA_length;

		return 1;


	}
	else  ///这就是正常展开，拿到sampled的位置
	{

		bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length);

		///如果已经拿到了所有的位置，则直接返回
		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}

	}


	///early leaf node calculation
	if (length_read >= 2)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta);



		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}


	}



	/**
	///early leaf node calculation
	if (length_read >= 2)
	{

	bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

	///这里不要改
	///这里SA要改
	if (bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
	sp1, ep1,
	&tmp_SA_length, bitmapper_index_params.delta))
	{
	if (tmp_SA_length == number_of_hits)
	{
	(*occurrences) = tmp_SA_length;

	return 1;
	}
	}
	else
	{
	tree_height = bitmapper_index_params.compress_sa;
	}




	}
	**/





	///队列置空
	empty_queue_FMtree();

	///程序能执行到这里，说明必须要根节点要扩展了
	///说明树的第一层和最后一层已经遍历完了，而且没有找到所有的位置
	///这个if判断说明整棵树大于两层即D>2
	if (current_layer + 1 < tree_height)
	{

		bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			&sp_0,
			&sp_1,
			&sp_2,
			&ep_0,
			&ep_1,
			&ep_2);





		if (ep_0 != sp_0)
		{
			enqueue_FMtree(sp_0, ep_0, current_layer + 1);
		}

		if (ep_1 != sp_1)
		{
			enqueue_FMtree(sp_1, ep_1, current_layer + 1);
		}

		if (ep_2 != sp_2)
		{
			enqueue_FMtree(sp_2, ep_2, current_layer + 1);
		}

	}








	while (queue_length_FMtree() != 0)
	{

		////节点出队
		dequeue_FMtree(&current_sp, &current_ep, &current_layer);

		///如果ep-sp在阈值之下，则所有位置one-to-one的拿到
		if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
		{

			///注意这里最多是向下扩展tree_height - current_layer - 1
			///而不是bitmapper_index_params.compress_sa - current_layer - 1
			///也就是D-2步
			///和前面的循环里的结果不一样
			///因为这里已经经过了early leaf node calculation, 所以应该是
			bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length,
				tree_height - current_layer - 1);


			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

			///如果在这个阈值之下,则其不可能再分支了
			continue;




		}
		else  ///这就是正常展开，拿到sampled的位置
		{

			bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length);

			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

		}




		if (current_layer + 1 < tree_height)
		{


			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_0,
				&sp_1,
				&sp_2,
				&ep_0,
				&ep_1,
				&ep_2);




			if (ep_0 != sp_0)
			{

				enqueue_FMtree(sp_0, ep_0, current_layer + 1);
			}

			if (ep_1 != sp_1)
			{
				enqueue_FMtree(sp_1, ep_1, current_layer + 1);
			}

			if (ep_2 != sp_2)
			{
				enqueue_FMtree(sp_2, ep_2, current_layer + 1);
			}

		}



	}





	(*occurrences) = tmp_SA_length;
}







inline bitmapper_bs_iter enqueue_FMtree_muti_thread
(bitmapper_bs_iter sp, bitmapper_bs_iter ep, bitmapper_bs_iter layer, bwt_locate_queue* get_queue)
{
	bitmapper_bs_iter queue_point = get_queue->FMtree_queue_end_point * 3;
	get_queue->FMtree_queue[queue_point++] = sp;
	get_queue->FMtree_queue[queue_point++] = ep;
	get_queue->FMtree_queue[queue_point] = layer;


	get_queue->FMtree_queue_end_point++;
}

inline bitmapper_bs_iter dequeue_FMtree_muti_thread
(bitmapper_bs_iter* sp, bitmapper_bs_iter* ep, bitmapper_bs_iter* layer, bwt_locate_queue* get_queue)

{

	bitmapper_bs_iter queue_point = get_queue->FMtree_queue_start_point * 3;
	*sp =
		get_queue->FMtree_queue[queue_point++];
	*ep =
		get_queue->FMtree_queue[queue_point++];


	*layer =
		get_queue->FMtree_queue[queue_point];


	get_queue->FMtree_queue_start_point++;
}

inline bitmapper_bs_iter queue_length_FMtree_muti_thread(bwt_locate_queue* get_queue)
{
	return (get_queue->FMtree_queue_end_point - get_queue->FMtree_queue_start_point);
}

inline bitmapper_bs_iter empty_queue_FMtree_muti_thread(bwt_locate_queue* get_queue)
{
	get_queue->FMtree_queue_start_point = 0;
	get_queue->FMtree_queue_end_point = 0;
}




bitmapper_bs_iter locate_muti_thread(char* pattern,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter sp1,
	bitmapper_bs_iter ep1,
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter length_read,
	bitmapper_bs_iter* occurrences,
	bwt_locate_queue* get_queue)
{

	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter number_of_hits = ep - sp;

	bitmapper_bs_iter current_sp;
	bitmapper_bs_iter current_ep;
	bitmapper_bs_iter current_layer;

	bitmapper_bs_iter sp_0;
	bitmapper_bs_iter sp_1;
	bitmapper_bs_iter sp_2;


	bitmapper_bs_iter ep_0;
	bitmapper_bs_iter ep_1;
	bitmapper_bs_iter ep_2;


	///因为用了early leaf node calculation
	bitmapper_bs_iter tree_height = bitmapper_index_params.compress_sa - 1;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;



	///如果ep-sp在阈值之下，则所有位置one-to-one的拿到，这样就可以直接返回了
	if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
	{

		///注意这里最多是向下扩展bitmapper_index_params.compress_sa - current_layer - 1
		///也就是D-1步
		///和后面的循环里的结果不一样
		///后面的循环里因为已经经过了early leaf node calculation, 所以应该是
		///tree_height - current_layer - 1
		bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length,
			bitmapper_index_params.compress_sa - current_layer - 1);

		(*occurrences) = tmp_SA_length;

		return 1;


	}
	else  ///这就是正常展开，拿到sampled的位置
	{

		bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length);

		///如果已经拿到了所有的位置，则直接返回
		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}

	}


	///early leaf node calculation
	if (length_read >= 2)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta);



		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}


	}



	/**
	///early leaf node calculation
	if (length_read >= 2)
	{

	bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

	///这里不要改
	///这里SA要改
	if (bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
	sp1, ep1,
	&tmp_SA_length, bitmapper_index_params.delta))
	{
	if (tmp_SA_length == number_of_hits)
	{
	(*occurrences) = tmp_SA_length;

	return 1;
	}
	}
	else
	{
	tree_height = bitmapper_index_params.compress_sa;
	}




	}
	**/





	///队列置空
	empty_queue_FMtree_muti_thread(get_queue);

	///程序能执行到这里，说明必须要根节点要扩展了
	///说明树的第一层和最后一层已经遍历完了，而且没有找到所有的位置
	///这个if判断说明整棵树大于两层即D>2
	if (current_layer + 1 < tree_height)
	{

		bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			&sp_0,
			&sp_1,
			&sp_2,
			&ep_0,
			&ep_1,
			&ep_2);





		if (ep_0 != sp_0)
		{
			enqueue_FMtree_muti_thread(sp_0, ep_0, current_layer + 1, get_queue);
		}

		if (ep_1 != sp_1)
		{
			enqueue_FMtree_muti_thread(sp_1, ep_1, current_layer + 1, get_queue);
		}

		if (ep_2 != sp_2)
		{
			enqueue_FMtree_muti_thread(sp_2, ep_2, current_layer + 1, get_queue);
		}

	}








	while (queue_length_FMtree_muti_thread(get_queue) != 0)
	{

		////节点出队
		dequeue_FMtree_muti_thread(&current_sp, &current_ep, &current_layer, get_queue);

		///如果ep-sp在阈值之下，则所有位置one-to-one的拿到
		if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
		{

			///注意这里最多是向下扩展tree_height - current_layer - 1
			///而不是bitmapper_index_params.compress_sa - current_layer - 1
			///也就是D-2步
			///和前面的循环里的结果不一样
			///因为这里已经经过了early leaf node calculation, 所以应该是
			bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length,
				tree_height - current_layer - 1);


			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

			///如果在这个阈值之下,则其不可能再分支了
			continue;




		}
		else  ///这就是正常展开，拿到sampled的位置
		{

			bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length);

			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

		}




		if (current_layer + 1 < tree_height)
		{


			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_0,
				&sp_1,
				&sp_2,
				&ep_0,
				&ep_1,
				&ep_2);




			if (ep_0 != sp_0)
			{

				enqueue_FMtree_muti_thread(sp_0, ep_0, current_layer + 1, get_queue);
			}

			if (ep_1 != sp_1)
			{
				enqueue_FMtree_muti_thread(sp_1, ep_1, current_layer + 1, get_queue);
			}

			if (ep_2 != sp_2)
			{
				enqueue_FMtree_muti_thread(sp_2, ep_2, current_layer + 1, get_queue);
			}

		}



	}





	(*occurrences) = tmp_SA_length;
}











bitmapper_bs_iter locate_one_by_one(char* pattern,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter sp1,
	bitmapper_bs_iter ep1,
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter length_read,
	bitmapper_bs_iter* occurrences)
{

	bitmapper_bs_iter tmp_SA_length = 0;
	bitmapper_bs_iter number_of_hits = ep - sp;

	bitmapper_bs_iter current_sp;
	bitmapper_bs_iter current_ep;
	bitmapper_bs_iter current_layer;


	bitmapper_bs_iter sp_0;
	bitmapper_bs_iter sp_1;
	bitmapper_bs_iter sp_2;


	bitmapper_bs_iter ep_0;
	bitmapper_bs_iter ep_1;
	bitmapper_bs_iter ep_2;


	///因为用了early leaf node calculation
	bitmapper_bs_iter tree_height = bitmapper_index_params.compress_sa - 1;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;



	bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
		bitmapper_index_params.SA_flag,
		bitmapper_index_params.bwt,
		bitmapper_index_params.high_occ_table,
		locates,
		current_sp,
		current_ep,
		current_layer,
		&tmp_SA_length,
		bitmapper_index_params.compress_sa - current_layer - 1);

	(*occurrences) = tmp_SA_length;

	return 1;

}








/**
unsigned int locate_debug(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences)
{

	unsigned int tmp_SA_length = 0;
	unsigned int number_of_hits = ep - sp;

	unsigned int current_sp;
	unsigned int current_ep;
	unsigned int current_layer;

	unsigned int sp_A;
	unsigned int sp_C;
	unsigned int sp_G;
	unsigned int sp_T;


	unsigned int ep_A;
	unsigned int ep_C;
	unsigned int ep_G;
	unsigned int ep_T;



	///因为用了early leaf node calculation
	unsigned int tree_height = bitmapper_index_params.compress_sa - 1;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;

	///如果ep-sp在阈值之下，则所有位置one-to-one的拿到，这样就可以直接返回了
	if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
	{

		///注意这里最多是向下扩展bitmapper_index_params.compress_sa - current_layer - 1
		///也就是D-1步
		///和后面的循环里的结果不一样
		///后面的循环里因为已经经过了early leaf node calculation, 所以应该是
		///tree_height - current_layer - 1
		bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length,
			bitmapper_index_params.compress_sa - current_layer - 1);

		(*occurrences) = tmp_SA_length;

		return 1;


	}
	else  ///这就是正常展开，拿到sampled的位置
	{

		bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length);

		///如果已经拿到了所有的位置，则直接返回
		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}

	}


	///early leaf node calculation
	if (length_read >= 2)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta);



		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}


	}






	///队列置空
	empty_queue_FMtree();

	///程序能执行到这里，说明必须要根节点要扩展了
	///说明树的第一层和最后一层已经遍历完了，而且没有找到所有的位置
	///这个if判断说明整棵树大于两层即D>2
	if (current_layer + 1 < tree_height)
	{
		debug_total++;
		if (current_ep - current_sp < 3
			&&
			current_ep - current_sp > 1)
		{
			debug_2++;
		}


		bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			&sp_A,
			&sp_C,
			&sp_G,
			&sp_T,
			&ep_A,
			&ep_C,
			&ep_G,
			&ep_T);



		if (ep_A != sp_A)
		{
			enqueue_FMtree(sp_A, ep_A, current_layer + 1);
		}

		if (ep_C != sp_C)
		{
			enqueue_FMtree(sp_C, ep_C, current_layer + 1);
		}

		if (ep_G != sp_G)
		{
			enqueue_FMtree(sp_G, ep_G, current_layer + 1);
		}
		if (ep_T != sp_T)
		{
			enqueue_FMtree(sp_T, ep_T, current_layer + 1);
		}
	}








	while (queue_length_FMtree() != 0)
	{

		////节点出队
		dequeue_FMtree(&current_sp, &current_ep, &current_layer);

		///如果ep-sp在阈值之下，则所有位置one-to-one的拿到
		if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
		{

			///注意这里最多是向下扩展tree_height - current_layer - 1
			///而不是bitmapper_index_params.compress_sa - current_layer - 1
			///也就是D-2步
			///和前面的循环里的结果不一样
			///因为这里已经经过了early leaf node calculation, 所以应该是
			bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length,
				tree_height - current_layer - 1);


			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

			///如果在这个阈值之下,则其不可能再分支了
			continue;




		}
		else  ///这就是正常展开，拿到sampled的位置
		{

			bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length);

			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

		}




		if (current_layer + 1 < tree_height)
		{
			debug_total++;
			if (current_ep - current_sp < 3
				&&
				current_ep - current_sp > 1)
			{
				debug_2++;
			}

			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_A,
				&sp_C,
				&sp_G,
				&sp_T,
				&ep_A,
				&ep_C,
				&ep_G,
				&ep_T);


			if (ep_A != sp_A)
			{

				enqueue_FMtree(sp_A, ep_A, current_layer + 1);
			}

			if (ep_C != sp_C)
			{
				enqueue_FMtree(sp_C, ep_C, current_layer + 1);
			}

			if (ep_G != sp_G)
			{
				enqueue_FMtree(sp_G, ep_G, current_layer + 1);
			}
			if (ep_T != sp_T)
			{
				enqueue_FMtree(sp_T, ep_T, current_layer + 1);
			}
		}



	}





	(*occurrences) = tmp_SA_length;
}

**/
/**
unsigned int locate_less_than_4(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences)
{

	unsigned int tmp_SA_length = 0;
	unsigned int number_of_hits = ep - sp;

	unsigned int current_sp;
	unsigned int current_ep;
	unsigned int current_layer;

	unsigned int sp_A;
	unsigned int sp_C;
	unsigned int sp_G;
	unsigned int sp_T;


	unsigned int ep_A;
	unsigned int ep_C;
	unsigned int ep_G;
	unsigned int ep_T;



	///因为不用early leaf node calculation
	unsigned int tree_height = bitmapper_index_params.compress_sa;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;

	///如果ep-sp在阈值之下，则所有位置one-to-one的拿到，这样就可以直接返回了
	if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
	{

		bwt_accesss_SA_cur_less_than_4(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length,
			tree_height - current_layer - 1);

		(*occurrences) = tmp_SA_length;

		return 1;


	}
	else  ///这就是正常展开，拿到sampled的位置
	{

		bwt_accesss_SA_less_than_4(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length);

		///如果已经拿到了所有的位置，则直接返回
		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}

	}





	///队列置空
	empty_queue_FMtree();


	if (current_layer + 1 < tree_height)
	{
		bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			&sp_A,
			&sp_C,
			&sp_G,
			&sp_T,
			&ep_A,
			&ep_C,
			&ep_G,
			&ep_T);



		if (ep_A != sp_A)
		{
			enqueue_FMtree(sp_A, ep_A, current_layer + 1);
		}

		if (ep_C != sp_C)
		{
			enqueue_FMtree(sp_C, ep_C, current_layer + 1);
		}

		if (ep_G != sp_G)
		{
			enqueue_FMtree(sp_G, ep_G, current_layer + 1);
		}
		if (ep_T != sp_T)
		{
			enqueue_FMtree(sp_T, ep_T, current_layer + 1);
		}
	}








	while (queue_length_FMtree() != 0)
	{

		////节点出队
		dequeue_FMtree(&current_sp, &current_ep, &current_layer);

		///如果ep-sp在阈值之下，则所有位置one-to-one的拿到
		if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
		{

			bwt_accesss_SA_cur_less_than_4(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length,
				tree_height - current_layer - 1);


			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

			///如果在这个阈值之下,则其不可能再分支了
			continue;




		}
		else  ///这就是正常展开，拿到sampled的位置
		{

			bwt_accesss_SA_less_than_4(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length);

			///如果已经拿到了所有的位置，则直接返回
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

		}




		if (current_layer + 1 < tree_height)
		{
			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_A,
				&sp_C,
				&sp_G,
				&sp_T,
				&ep_A,
				&ep_C,
				&ep_G,
				&ep_T);


			if (ep_A != sp_A)
			{

				enqueue_FMtree(sp_A, ep_A, current_layer + 1);
			}

			if (ep_C != sp_C)
			{
				enqueue_FMtree(sp_C, ep_C, current_layer + 1);
			}

			if (ep_G != sp_G)
			{
				enqueue_FMtree(sp_G, ep_G, current_layer + 1);
			}
			if (ep_T != sp_T)
			{
				enqueue_FMtree(sp_T, ep_T, current_layer + 1);
			}
		}



	}





	(*occurrences) = tmp_SA_length;
}
**/

/**
unsigned int locate(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences)
{

	unsigned int tmp_SA_length = 0;
	unsigned int number_of_hits = ep - sp;
	unsigned int tree_layers_i = 0;
	unsigned int father_sp;
	unsigned int father_ep;

	unsigned int occ = 0;


	bitmapper_index_params.tree_index = 0;

	bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index] = sp;
	bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] = ep;

	bitmapper_index_params.tree_layers = 0;

	bitmapper_index_params.need_step = bitmapper_index_params.compress_sa - bitmapper_index_params.tree_layers - 1;

	if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] -
		bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index] <= bitmapper_index_params.cut_thr)
	{
		///这个函数内部要改
		///这里SA要改
		bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
			locates,
			bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
			bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
			&tmp_SA_length, bitmapper_index_params.need_step);

		bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
			= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] = (unsigned int)-1;
	}
	else
	{
		///这里不要改
		///这里SA要改
		bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
			bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
			&tmp_SA_length);
	}




	if (length_read >= 2 &&
		tmp_SA_length != number_of_hits)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		///这里不要改
		///这里SA要改
		bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta);


	}






	bitmapper_index_params.tree_layer_length = 1;







	bitmapper_index_params.tree_layer_length = 1;

	bitmapper_index_params.tree_index = 1;

	for (bitmapper_index_params.tree_layers = 1;
		bitmapper_index_params.tree_layers < bitmapper_index_params.compress_sa - 1;
		bitmapper_index_params.tree_layers++)
	{


		if (tmp_SA_length == number_of_hits)
		{
			break;
		}


		//need_step = compress_sa - tree_layers - 1;
		bitmapper_index_params.need_step
			= bitmapper_index_params.compress_sa - bitmapper_index_params.tree_layers - 2;

		occ = bitmapper_index_params.tree_layers;


		for (tree_layers_i = 0; tree_layers_i < bitmapper_index_params.tree_layer_length; tree_layers_i++)
		{

			father_sp = bitmapper_index_params.sp_tree
				[(bitmapper_index_params.tree_index - 1) / 4];
			father_ep = bitmapper_index_params.ep_tree
				[(bitmapper_index_params.tree_index - 1) / 4];





			if (father_sp == (unsigned int)-1)
			{
				bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
					= bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1]
					= bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2]
					= bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3]
					= (unsigned int)-1;
			}
			else
			{
				///这里要改
				bwt_find_occ_all_sp_ep_optimal(father_sp, father_ep, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]);





				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
					<= bitmapper_index_params.cut_thr)
				{
					///这里要改
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);




					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index]);

					fflush(stderr);



					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] = (unsigned int)-1;
				}
				else
				{
					///这里不要改
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
						&tmp_SA_length);



					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index]);

					fflush(stderr);

				}




				if (tmp_SA_length == number_of_hits)
				{

					break;
				}


				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1]
					<= bitmapper_index_params.cut_thr)
				{

					///这里要改
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);



					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1]);

					fflush(stderr);


					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1] = (unsigned int)-1;
				}
				else
				{
					///这里不要改
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1], occ,
						&tmp_SA_length);


					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1]);

					fflush(stderr);

				}




				if (tmp_SA_length == number_of_hits)
				{
					break;
				}


				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2]
					<= bitmapper_index_params.cut_thr)
				{

					///这里要改
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);




					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2]);

					fflush(stderr);




					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2] = (unsigned int)-1;
				}
				else
				{
					///这里不要改
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2], occ,
						&tmp_SA_length);


					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2]);

					fflush(stderr);

				}




				if (tmp_SA_length == number_of_hits)
				{
					break;
				}

				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3]
					<= bitmapper_index_params.cut_thr)
				{

					///这里要改
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);


					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]);

					fflush(stderr);


					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]
						= (unsigned int)-1;
				}
				else
				{
					///这里不要改
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3], occ,
						&tmp_SA_length);



					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]);

					fflush(stderr);

				}



				if (tmp_SA_length == number_of_hits)
				{
					break;
				}

			}

			bitmapper_index_params.tree_index = bitmapper_index_params.tree_index + 4;
		}

		bitmapper_index_params.tree_layer_length = bitmapper_index_params.tree_layer_length * 4;
	}

	(*occurrences) = tmp_SA_length;
}

**/

void debug_information()
{
	fprintf(stderr, "debug_2 = %llu \n", debug_2);
	fprintf(stderr, "debug_total = %llu \n", debug_total);
}
