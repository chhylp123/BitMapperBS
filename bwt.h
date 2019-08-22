
#ifndef __BWT__
#define __BWT__


#include "uint40.h"

///#define SA_counter_length 32
#define SA_counter_length 64
typedef uint64_t bwt_string_type;
typedef uint64_t SA_flag_string_type;
///typedef uint32_t high_occ_table_type;
typedef uint64_t high_occ_table_type;
typedef uint64_t bitmapper_bs_iter;


#define likely(x) __builtin_expect (!!(x), 1)
#define unlikely(x) __builtin_expect (!!(x), 0)

typedef struct bwt_locate_queue
{
	bitmapper_bs_iter* FMtree_queue;
	bitmapper_bs_iter FMtree_queue_start_point;
	bitmapper_bs_iter FMtree_queue_end_point;


	bwt_locate_queue():FMtree_queue_start_point(0), FMtree_queue_end_point(0)
	{};

} bwt_locate_queue;



typedef struct bwt_index
{

	unsigned char hash_count[65536];



	unsigned int *sa;
	SA_flag_string_type* SA_flag;
	bwt_string_type *bwt;
	high_occ_table_type* high_occ_table;




	///unsigned int SA_length;
	bitmapper_bs_iter SA_length;
	///unsigned int bwt_length;
	bitmapper_bs_iter bwt_length;
	bitmapper_bs_iter sparse_suffix_array_length;///unsigned int sparse_suffix_array_length;
	bitmapper_bs_iter SA_flag_iterater;///unsigned int SA_flag_iterater;
	///unsigned int high_occ_table_length;
	bitmapper_bs_iter high_occ_table_length;


	///unsigned int shapline;
	bitmapper_bs_iter shapline;


	unsigned int compress_sa;
	///unsigned int compress_occ = 256;
	unsigned int compress_occ;
	///unsigned int high_compress_occ = 65536;
	///unsigned int high_compress_occ = 256;
	unsigned int high_compress_occ;
	unsigned int idependent_high_compress_occ;
	///unsigned int compress_SA_flag = 224;
	unsigned int compress_SA_flag;
	unsigned int bwt_warp_number;
	unsigned int SA_flag_warp_number;
	uint64_t tmp_SA_flag;
	uint64_t* long_SA_flag;
	unsigned int single_occ_bwt;
	unsigned int occ_words;
	///unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned short)* 4 * 8 / 2;
	///unsigned int acctuall_bwt_gap = high_compress_occ + 64;
	unsigned int acctuall_bwt_gap;
	unsigned int acctuall_SA_flag_gap;
	unsigned int SA_counter_shift_length;
	bwt_string_type mode_4[4];
	bwt_string_type mode;
	bwt_string_type mode_high_1;
	bwt_string_type mode_16;
	bwt_string_type mode_32;
	///bitmapper_bs_iter mode_low8 = (bwt_string_type)255;
	bitmapper_bs_iter mode_low32;
	bwt_string_type mode_high;
	bwt_string_type mode_low;
	SA_flag_string_type mode_SA_flag;
	bwt_string_type mode_bwt_delta;
	SA_flag_string_type SA_pop_count_mode;
	bwt_string_type pop_count_mode[4];
	unsigned int bwt_count_hash_table_bit;
	unsigned int text_length;
	unsigned int na, nc, ng, nt;
	bitmapper_bs_iter nacgt[5];
	unsigned int bwt_step;
	unsigned int ctoi[256];
	char itoc[4];

	unsigned int total_gap;
	unsigned int num_r, len_r;
	unsigned int SA_header_mode;
	unsigned int SA_header_mode_reverse;
	FILE* _rg_fp;


	int delta;



	unsigned int tree_nodes_number;
	unsigned int* sp_tree;
	unsigned int* ep_tree;
	unsigned int tree_index;
	unsigned int tree_layers;
	unsigned int need_step;
	unsigned int cut_thr;
	unsigned int tree_layer_length;
	bitmapper_bs_iter* FMtree_queue;///unsigned int* FMtree_queue;
	///char* FMtree_queue_int;

	bitmapper_bs_iter FMtree_queue_start_point;
	bitmapper_bs_iter FMtree_queue_end_point;


	bitmapper_bs_iter c_0_times;
	bitmapper_bs_iter c_1_times;
	bitmapper_bs_iter c_2_times;


	unsigned int* hash_table_16_mer_high_32;
	unsigned int hash_table_shift_length;
	unsigned int hash_table_mode;
	unsigned int hash_table_mode_get;
	unsigned char* hash_table_16_mer_low_8;
	bitmapper_bs_iter hash_table_16_mer_size;


	bwt_index() :compress_sa(8), compress_occ(64), high_compress_occ(128), idependent_high_compress_occ(65536),
		compress_SA_flag(256), bwt_warp_number(sizeof(bwt_string_type)* 8 / 2), SA_flag_warp_number(sizeof(SA_flag_string_type)* 8),
		tmp_SA_flag((uint64_t)0), long_SA_flag(NULL), single_occ_bwt(sizeof(bwt_string_type) / sizeof(unsigned short)),
		occ_words(4 / single_occ_bwt), acctuall_bwt_gap(high_compress_occ + 32), acctuall_SA_flag_gap(compress_SA_flag + SA_counter_length),
		SA_counter_shift_length(sizeof(SA_flag_string_type)* 8 - SA_counter_length), mode((bwt_string_type)-1), mode_high_1((bwt_string_type)1 << (SA_flag_warp_number - 1)),
		mode_16((bwt_string_type)65535), mode_32(((bwt_string_type)-1) >> 32), mode_low32 (((bwt_string_type)-1) >> 32),
		mode_SA_flag((SA_flag_string_type)(((SA_flag_string_type)-1) << (sizeof(SA_flag_string_type)* 8 - 1))),
		mode_bwt_delta ((bwt_string_type)(((bwt_string_type)-1) << (sizeof(bwt_string_type)* 8 - 1))),
		SA_pop_count_mode ((SA_flag_string_type)(((SA_flag_string_type)-1) >> SA_counter_length)),
		bwt_count_hash_table_bit(16), bwt_step(1), itoc({ 'A', 'C', 'G', 'T' }),
		num_r(10000), len_r(12),SA_header_mode((unsigned int)(((unsigned int)-1) >> 2)),
		SA_header_mode_reverse((unsigned int)(((unsigned int)-1) << 30)), cut_thr(1),
		tree_layer_length(1), FMtree_queue_start_point(0), FMtree_queue_end_point(0),
		c_0_times(0), c_1_times(0), c_2_times(0),
		hash_table_shift_length(28),hash_table_mode(((unsigned int)-1) << hash_table_shift_length),
		hash_table_mode_get(((unsigned int)-1) >> (32 - hash_table_shift_length))
		{};



} bwt_index;




typedef struct SA_block
{
	bitmapper_bs_iter start, endl, block_size, text_length;
	uint40* sa_block;
	FILE *f_sa;


}SA_block;

extern bwt_index bitmapper_index_params;


//void dec_bit_bwt(bwt_string_type input);


/**
unsigned int locate(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences);
**/


bitmapper_bs_iter locate_one_by_one(char* pattern,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter sp1,
	bitmapper_bs_iter ep1,
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter length_read,
	bitmapper_bs_iter* occurrences);

bitmapper_bs_iter locate(char* pattern, bitmapper_bs_iter sp, bitmapper_bs_iter ep,
	bitmapper_bs_iter sp1, bitmapper_bs_iter ep1, bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter length_read, bitmapper_bs_iter* occurrences);

unsigned int init_locate_queue_muti_thread(bwt_locate_queue* get_queue);

bitmapper_bs_iter locate_muti_thread(char* pattern,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter ep,
	bitmapper_bs_iter sp1,
	bitmapper_bs_iter ep1,
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter length_read,
	bitmapper_bs_iter* occurrences,
	bwt_locate_queue* get_queue);
/**
unsigned int locate_debug(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences);

unsigned int locate_less_than_4(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences);
**/
unsigned int init_bitmapper_index_params();

unsigned int indenpendent_creadte_index(bitmapper_bs_iter text_length, char** input_refer,
	bitmapper_bs_iter compress_sa, char* filename);

/**
unsigned int count(char* pattern, unsigned int length,
	unsigned int* sp, unsigned int* ep, unsigned int* sp1, unsigned int* ep1);
**/



///unsigned int find_occ_fm_index(unsigned int line, int delta, bwt_string_type *bwt, high_occ_table_type* high_occ_table);
/**
void bwt_accesss_SA_cur_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length, unsigned int need_step);

**/

/**
unsigned int bwt_get_sa_restrict_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa);
**/

/**
void bwt_direct_get_sa_interval_long
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length);
**/
/**
void bwt_accesss_SA_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length);
**/

/**
unsigned int bwt_get_sa_restrict_zero_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start);
**/


/**
void bwt_accesss_pre_SA(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep,
	unsigned int* tmp_SA_length, unsigned int delta);
**/

/**
unsigned int bwt_get_sa_restrict_zero_steps
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start);
**/

/**
void bwt_find_occ_all_sp_ep_optimal(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep);
**/

unsigned int load_index(char* filename_prefix);


inline void query_16_mer_hash_table(bitmapper_bs_iter hash_key,
	bitmapper_bs_iter* sp, bitmapper_bs_iter* ep)
{
	bitmapper_bs_iter tmp_high, tmp_low;
	tmp_high = bitmapper_index_params.hash_table_16_mer_high_32[hash_key] & bitmapper_index_params.hash_table_mode_get;
	tmp_high = tmp_high << 8;
	tmp_low = bitmapper_index_params.hash_table_16_mer_low_8[hash_key];
	tmp_high = tmp_high | tmp_low;
	*sp = tmp_high;



	tmp_high = bitmapper_index_params.hash_table_16_mer_high_32[hash_key + 1] & bitmapper_index_params.hash_table_mode_get;
	tmp_high = tmp_high << 8;
	tmp_low = bitmapper_index_params.hash_table_16_mer_low_8[hash_key + 1];
	tmp_high = tmp_high | tmp_low;
	*ep = tmp_high;

	tmp_high = bitmapper_index_params.hash_table_16_mer_high_32[hash_key + 1] & bitmapper_index_params.hash_table_mode;
	tmp_high = tmp_high >> bitmapper_index_params.hash_table_shift_length;
	*ep = *ep - tmp_high;

}


inline bitmapper_bs_iter get_3_letter_hash_value(char* pattern, unsigned int length)
{
	bitmapper_bs_iter hash_value = 0;

	bitmapper_bs_iter convert_i = 0;

	bitmapper_bs_iter delta;

	for (convert_i = 0; convert_i < length; convert_i++)
	{
		hash_value = hash_value * 3;
		delta = bitmapper_index_params.ctoi[pattern[convert_i]];
		if (delta > 2)
		{
			return ((bitmapper_bs_iter)-1);
		}
		hash_value = hash_value + delta;
		///hash_value = hash_value + bitmapper_index_params.ctoi[pattern[convert_i]];

	}

	return hash_value;

}

inline bitmapper_bs_iter get_3_letter_hash_value_N_site(char* pattern, unsigned int length, bitmapper_bs_iter* N_site)
{
	bitmapper_bs_iter hash_value = 0;

	bitmapper_bs_iter convert_i = 0;

	bitmapper_bs_iter delta;

	for (convert_i = 0; convert_i < length; convert_i++)
	{
		hash_value = hash_value * 3;
		delta = bitmapper_index_params.ctoi[pattern[convert_i]];
		if (delta > 2)
		{
			(*N_site) = convert_i;
			return ((bitmapper_bs_iter)-1);
		}
		hash_value = hash_value + delta;
		///hash_value = hash_value + bitmapper_index_params.ctoi[pattern[convert_i]];

	}

	return hash_value;

}






///注意字符0的high_occ_value都是没减过的值，而occ_value是减过的值
inline bitmapper_bs_iter get_occ_value_sp_back(bitmapper_bs_iter line, int delta,
                                          bitmapper_bs_iter* actually_line, bitmapper_bs_iter* high_occ_actually_line,
                                          bitmapper_bs_iter* occ_value, bitmapper_bs_iter* high_occ_value)
{

	///high_occ_iter = ((line / bitmapper_index_params.high_compress_occ)*bitmapper_index_params.acctuall_bwt_gap)/ bitmapper_index_params.bwt_warp_number;
	///bitmapper_bs_iter high_occ_iter = ((line >> 8) * bitmapper_index_params.acctuall_bwt_gap)>>5;
	bitmapper_bs_iter high_occ_iter = ((line >> 8) * 320) >> 5;


	///bitmapper_bs_iter low_occ_iter = (bwt_iter_i % bitmapper_index_params.high_compress_occ) / bitmapper_index_params.compress_occ;
	bitmapper_bs_iter low_occ_iter = (line & 255) >> 6;



	(*actually_line) = high_occ_iter + 2 + (low_occ_iter << 1);
	(*high_occ_actually_line) = high_occ_iter;




	bitmapper_bs_iter tmp_counter_high32 = (bwt_string_type)0;
	bitmapper_bs_iter tmp_counter_low8 = (bwt_string_type)0;
	bitmapper_bs_iter tmp_counter_40 = (bwt_string_type)0;



	///1代表是字符T
	if (delta == 1)
	{
		///先取T的32-bit counter,因为T的东西都存在高32-bit
		tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] >> 32;
		///再取T的8-bit counter
		tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 56;
		tmp_counter_40 = tmp_counter_high32 << 8;
		tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;

		(*high_occ_value) = tmp_counter_40;

		///low_occ_iter=0根本不用处理
		/**
		if(low_occ_iter==1)
        {
            tmp_counter_40 = tmp_counter_40 +
                ((bitmapper_index_params.bwt[high_occ_iter + 1] >> 40) & bitmapper_index_params.mode_low8);
        }else if(low_occ_iter==2)
        {
            tmp_counter_40 = tmp_counter_40 +
                ((bitmapper_index_params.bwt[high_occ_iter + 1] >> 24) & bitmapper_index_params.mode_low8);
        }else if(low_occ_iter==3)
        {
            tmp_counter_40 = tmp_counter_40 +
                ((bitmapper_index_params.bwt[high_occ_iter + 1] >> 8) & bitmapper_index_params.mode_low8);
        }
		**/

		if (low_occ_iter != 0)
		{
			/**
			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> (56 - (low_occ_iter << 4))) & bitmapper_index_params.mode_low8);
			**/

			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> (56 - (low_occ_iter << 4))) & 255);

		}


        (*occ_value) = tmp_counter_40;

		///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);

        return tmp_counter_40;

	}
	else if (delta == 2)  ///2代表是字符G
	{


	    ///先取G的32-bit counter
        tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] & bitmapper_index_params.mode_low32;
        ///再取G的8-bit counter
        tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 48;
        ///tmp_counter_low8 = tmp_counter_low8 & bitmapper_index_params.mode_low8;
		tmp_counter_low8 = tmp_counter_low8 & 255;
        ///最后就是完整的40-bit counter
        tmp_counter_40 = tmp_counter_high32 << 8;
        tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;

        (*high_occ_value) = tmp_counter_40;


		/**
        ///low_occ_iter=0根本不用处理
		if(low_occ_iter==1)
        {
            tmp_counter_40 = tmp_counter_40 +
                ((bitmapper_index_params.bwt[high_occ_iter + 1] >> 32) & bitmapper_index_params.mode_low8);
        }else if(low_occ_iter==2)
        {
            tmp_counter_40 = tmp_counter_40 +
                ((bitmapper_index_params.bwt[high_occ_iter + 1] >> 16) & bitmapper_index_params.mode_low8);
        }else if(low_occ_iter==3)
        {
            tmp_counter_40 = tmp_counter_40 +
                ((bitmapper_index_params.bwt[high_occ_iter + 1]) & bitmapper_index_params.mode_low8);
        }
		**/

		if (low_occ_iter != 0)
		{
			/**
			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4))) & bitmapper_index_params.mode_low8);
			**/

			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4))) & 255);
		}


		(*occ_value) = tmp_counter_40;

		///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);

        return tmp_counter_40;
	}
	else if (delta == 0)  ///0代表是字符A
	{

		/****************************T的high occ value*********************************************/
		///先取T的32-bit counter,因为T的东西都存在高32-bit
		tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] >> 32;
		///再取T的8-bit counter
		tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 56;
		tmp_counter_40 = tmp_counter_high32 << 8;
		tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;
		/****************************T的high occ value*********************************************/


		/****************************G的high occ value*********************************************/
		///先取G的32-bit counter
		tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] & bitmapper_index_params.mode_low32;
		///再取G的8-bit counter
		tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 48;
		///tmp_counter_low8 = tmp_counter_low8 & bitmapper_index_params.mode_low8;
		tmp_counter_low8 = tmp_counter_low8 & 255;
		///最后就是完整的40-bit counter
		tmp_counter_high32 = tmp_counter_high32 << 8;
		tmp_counter_high32 = tmp_counter_high32 | tmp_counter_low8;
		/****************************G的high occ value*********************************************/

		/****************************A+G的high occ value*********************************************/
		tmp_counter_40 = tmp_counter_40 + tmp_counter_high32;
		/****************************A+G的high occ value*********************************************/

        ///注意字符0的high_occ_value是没减过的值
		(*high_occ_value) = tmp_counter_40;

		/**
		///low_occ_iter=0根本不用处理
		if (low_occ_iter == 1)
		{
			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> 40) & bitmapper_index_params.mode_low8) +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> 32) & bitmapper_index_params.mode_low8);
		}
		else if (low_occ_iter == 2)
		{
			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> 24) & bitmapper_index_params.mode_low8) +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> 16) & bitmapper_index_params.mode_low8);
		}
		else if (low_occ_iter == 3)
		{
			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> 8) & bitmapper_index_params.mode_low8) +
				((bitmapper_index_params.bwt[high_occ_iter + 1]) & bitmapper_index_params.mode_low8);
		}
		**/






		if (low_occ_iter != 0)
		{
			tmp_counter_high32 = (bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4)));
			/**
			tmp_counter_40 = tmp_counter_40 + (tmp_counter_high32 & bitmapper_index_params.mode_low8) +
				((tmp_counter_high32 >> 8) & bitmapper_index_params.mode_low8);
			**/
			tmp_counter_40 = tmp_counter_40 + (tmp_counter_high32 & 255) + ((tmp_counter_high32 >> 8) & 255);
		}





        ///注意字符0的occ_value是没减过的值
        (*occ_value) = ((line >> 6) << 6) - tmp_counter_40;



		///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);


		///return ((line >> 6) << 6) - tmp_counter_40;
		return tmp_counter_40;


	}
}





///注意字符0的high_occ_value_sp都是没减过的值，而occ_value_sp是减过的值
inline bitmapper_bs_iter get_occ_value_ep_back(bitmapper_bs_iter line,
                                          int delta,
                                          bitmapper_bs_iter actually_line_sp,
                                          bitmapper_bs_iter high_occ_actually_line_sp,
                                          bitmapper_bs_iter occ_value_sp,
                                          bitmapper_bs_iter high_occ_value_sp,
                                          bitmapper_bs_iter* actually_line)
{


	bitmapper_bs_iter high_occ_iter = ((line >> 8) * 320) >> 5;
	bitmapper_bs_iter low_occ_iter = (line & 255) >> 6;
	(*actually_line) = high_occ_iter + 2 + (low_occ_iter << 1);

	bitmapper_bs_iter tmp_counter_high32 = (bwt_string_type)0;
	bitmapper_bs_iter tmp_counter_low8 = (bwt_string_type)0;
	bitmapper_bs_iter tmp_counter_40 = (bwt_string_type)0;


    ///这种情况下sp和ep在同一个block内，直接返回就好了
	if((*actually_line) == actually_line_sp)
    {
        return occ_value_sp;
    }
    else if(high_occ_iter == high_occ_actually_line_sp) ///这种情况下sp和ep在同一个super block内，只需要计算小block的occ value
    {


        tmp_counter_40 = high_occ_value_sp;


        if (low_occ_iter != 0)
        {
            if (delta == 1)
            {

                tmp_counter_40 = tmp_counter_40 +
                        ((bitmapper_index_params.bwt[high_occ_iter + 1] >> (56 - (low_occ_iter << 4))) & 255);

                ///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);
                return tmp_counter_40;
            }
            else if (delta == 2)  ///2代表是字符G
            {
                tmp_counter_40 = tmp_counter_40 +
                        ((bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4))) & 255);

                ///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);
                return tmp_counter_40;
            }
            else if (delta == 0)  ///0代表是字符A
            {

                tmp_counter_high32 = (bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4)));
                tmp_counter_40 = tmp_counter_40 + (tmp_counter_high32 & 255) + ((tmp_counter_high32 >> 8) & 255);

                ///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);
                return ((line >> 6) << 6) - tmp_counter_40;
            }

        }
        else
        {
            if (delta == 0)  ///0代表是字符A
            {
                tmp_counter_40 = ((line >> 6) << 6) - tmp_counter_40;
            }

            return tmp_counter_40;
        }

    }
    else    ///这种情况下sp和ep也不在同一个super block内
    {


        ///1代表是字符T
        if (delta == 1)
        {
            ///先取T的32-bit counter,因为T的东西都存在高32-bit
            tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] >> 32;
            ///再取T的8-bit counter
            tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 56;
            tmp_counter_40 = tmp_counter_high32 << 8;
            tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;


            if (low_occ_iter != 0)
            {

                tmp_counter_40 = tmp_counter_40 +
                    ((bitmapper_index_params.bwt[high_occ_iter + 1] >> (56 - (low_occ_iter << 4))) & 255);
            }

            ///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);

            return tmp_counter_40;

        }
        else if (delta == 2)  ///2代表是字符G
        {


            ///先取G的32-bit counter
            tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] & bitmapper_index_params.mode_low32;
            ///再取G的8-bit counter
            tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 48;
            ///tmp_counter_low8 = tmp_counter_low8 & bitmapper_index_params.mode_low8;
            tmp_counter_low8 = tmp_counter_low8 & 255;
            ///最后就是完整的40-bit counter
            tmp_counter_40 = tmp_counter_high32 << 8;
            tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;


            if (low_occ_iter != 0)
            {

                tmp_counter_40 = tmp_counter_40 +
                    ((bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4))) & 255);
            }

            ///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);

            return tmp_counter_40;
        }
        else if (delta == 0)  ///0代表是字符A
        {

            /****************************T的high occ value*********************************************/
            ///先取T的32-bit counter,因为T的东西都存在高32-bit
            tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] >> 32;
            ///再取T的8-bit counter
            tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 56;
            tmp_counter_40 = tmp_counter_high32 << 8;
            tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;
            /****************************T的high occ value*********************************************/


            /****************************G的high occ value*********************************************/
            ///先取G的32-bit counter
            tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] & bitmapper_index_params.mode_low32;
            ///再取G的8-bit counter
            tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 48;
            ///tmp_counter_low8 = tmp_counter_low8 & bitmapper_index_params.mode_low8;
            tmp_counter_low8 = tmp_counter_low8 & 255;
            ///最后就是完整的40-bit counter
            tmp_counter_high32 = tmp_counter_high32 << 8;
            tmp_counter_high32 = tmp_counter_high32 | tmp_counter_low8;
            /****************************G的high occ value*********************************************/

            /****************************A+G的high occ value*********************************************/
            tmp_counter_40 = tmp_counter_40 + tmp_counter_high32;
            /****************************A+G的high occ value*********************************************/



            if (low_occ_iter != 0)
            {
                tmp_counter_high32 = (bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4)));
                tmp_counter_40 = tmp_counter_40 + (tmp_counter_high32 & 255) + ((tmp_counter_high32 >> 8) & 255);
            }


            ///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);


            return ((line >> 6) << 6) - tmp_counter_40;


        }

    }



}
































inline bitmapper_bs_iter get_occ_value_back(bitmapper_bs_iter line, int delta, bitmapper_bs_iter* actually_line)
{

	///high_occ_iter = ((line / bitmapper_index_params.high_compress_occ)*bitmapper_index_params.acctuall_bwt_gap)/ bitmapper_index_params.bwt_warp_number;
	///bitmapper_bs_iter high_occ_iter = ((line >> 8) * bitmapper_index_params.acctuall_bwt_gap)>>5;
	bitmapper_bs_iter high_occ_iter = ((line >> 8) * 320) >> 5;


	///bitmapper_bs_iter low_occ_iter = (bwt_iter_i % bitmapper_index_params.high_compress_occ) / bitmapper_index_params.compress_occ;
	bitmapper_bs_iter low_occ_iter = (line & 255) >> 6;



	(*actually_line) = high_occ_iter + 2 + (low_occ_iter << 1);




	bitmapper_bs_iter tmp_counter_high32 = (bwt_string_type)0;
	bitmapper_bs_iter tmp_counter_low8 = (bwt_string_type)0;
	bitmapper_bs_iter tmp_counter_40 = (bwt_string_type)0;



	///1代表是字符T
	if (delta == 1)
	{
		///先取T的32-bit counter,因为T的东西都存在高32-bit
		tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] >> 32;
		///再取T的8-bit counter
		tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 56;
		tmp_counter_40 = tmp_counter_high32 << 8;
		tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;

		///low_occ_iter=0根本不用处理
		/**
		if(low_occ_iter==1)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 40) & bitmapper_index_params.mode_low8);
		}else if(low_occ_iter==2)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 24) & bitmapper_index_params.mode_low8);
		}else if(low_occ_iter==3)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 8) & bitmapper_index_params.mode_low8);
		}
		**/

		if (low_occ_iter != 0)
		{
			/**
			tmp_counter_40 = tmp_counter_40 +
			((bitmapper_index_params.bwt[high_occ_iter + 1] >> (56 - (low_occ_iter << 4))) & bitmapper_index_params.mode_low8);
			**/

			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> (56 - (low_occ_iter << 4))) & 255);

		}




		///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);

		return tmp_counter_40;

	}
	else if (delta == 2)  ///2代表是字符G
	{


		///先取G的32-bit counter
		tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] & bitmapper_index_params.mode_low32;
		///再取G的8-bit counter
		tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 48;
		///tmp_counter_low8 = tmp_counter_low8 & bitmapper_index_params.mode_low8;
		tmp_counter_low8 = tmp_counter_low8 & 255;
		///最后就是完整的40-bit counter
		tmp_counter_40 = tmp_counter_high32 << 8;
		tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;


		/**
		///low_occ_iter=0根本不用处理
		if(low_occ_iter==1)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 32) & bitmapper_index_params.mode_low8);
		}else if(low_occ_iter==2)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 16) & bitmapper_index_params.mode_low8);
		}else if(low_occ_iter==3)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1]) & bitmapper_index_params.mode_low8);
		}
		**/

		if (low_occ_iter != 0)
		{
			/**
			tmp_counter_40 = tmp_counter_40 +
			((bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4))) & bitmapper_index_params.mode_low8);
			**/

			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4))) & 255);
		}


		///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);

		return tmp_counter_40;
	}
	else if (delta == 0)  ///0代表是字符A
	{

		/****************************T的high occ value*********************************************/
		///先取T的32-bit counter,因为T的东西都存在高32-bit
		tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] >> 32;
		///再取T的8-bit counter
		tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 56;
		tmp_counter_40 = tmp_counter_high32 << 8;
		tmp_counter_40 = tmp_counter_40 | tmp_counter_low8;
		/****************************T的high occ value*********************************************/


		/****************************G的high occ value*********************************************/
		///先取G的32-bit counter
		tmp_counter_high32 = bitmapper_index_params.bwt[high_occ_iter] & bitmapper_index_params.mode_low32;
		///再取G的8-bit counter
		tmp_counter_low8 = bitmapper_index_params.bwt[high_occ_iter + 1] >> 48;
		///tmp_counter_low8 = tmp_counter_low8 & bitmapper_index_params.mode_low8;
		tmp_counter_low8 = tmp_counter_low8 & 255;
		///最后就是完整的40-bit counter
		tmp_counter_high32 = tmp_counter_high32 << 8;
		tmp_counter_high32 = tmp_counter_high32 | tmp_counter_low8;
		/****************************G的high occ value*********************************************/

		/****************************A+G的high occ value*********************************************/
		tmp_counter_40 = tmp_counter_40 + tmp_counter_high32;
		/****************************A+G的high occ value*********************************************/



		/**
		///low_occ_iter=0根本不用处理
		if (low_occ_iter == 1)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 40) & bitmapper_index_params.mode_low8) +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 32) & bitmapper_index_params.mode_low8);
		}
		else if (low_occ_iter == 2)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 24) & bitmapper_index_params.mode_low8) +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 16) & bitmapper_index_params.mode_low8);
		}
		else if (low_occ_iter == 3)
		{
		tmp_counter_40 = tmp_counter_40 +
		((bitmapper_index_params.bwt[high_occ_iter + 1] >> 8) & bitmapper_index_params.mode_low8) +
		((bitmapper_index_params.bwt[high_occ_iter + 1]) & bitmapper_index_params.mode_low8);
		}
		**/






		if (low_occ_iter != 0)
		{
			tmp_counter_high32 = (bitmapper_index_params.bwt[high_occ_iter + 1] >> (48 - (low_occ_iter << 4)));
			/**
			tmp_counter_40 = tmp_counter_40 + (tmp_counter_high32 & bitmapper_index_params.mode_low8) +
			((tmp_counter_high32 >> 8) & bitmapper_index_params.mode_low8);
			**/
			tmp_counter_40 = tmp_counter_40 + (tmp_counter_high32 & 255) + ((tmp_counter_high32 >> 8) & 255);
		}










		///__builtin_prefetch(&(bitmapper_index_params.bwt[(*actually_line)]), 0, 1);


		return ((line >> 6) << 6) - tmp_counter_40;


	}
}







///inline bitmapper_bs_iter get_occ_value(bitmapper_bs_iter line, int delta)
inline bitmapper_bs_iter get_occ_value(bitmapper_bs_iter line, bitmapper_bs_iter delta, bitmapper_bs_iter* actually_line)
{

	///high_occ_iter = ((line / bitmapper_index_params.high_compress_occ)*bitmapper_index_params.acctuall_bwt_gap)/ bitmapper_index_params.bwt_warp_number;
	bitmapper_bs_iter high_occ_iter = ((line >> 7) * 160) >> 5;


	bitmapper_bs_iter low_occ_iter = (line & 127) >> 6;

	(*actually_line) = high_occ_iter + 1 + (low_occ_iter << 1);

	bitmapper_bs_iter tmp_counter_40 = (bwt_string_type)0;

	/**
	if (delta == 1)
	{
		tmp_counter_40 = bitmapper_index_params.high_occ_table[((line >> 16) << 1)];

		tmp_counter_40 = tmp_counter_40 + ((bitmapper_index_params.bwt[high_occ_iter] >> (48 - (low_occ_iter << 5))) & 65535);

		return tmp_counter_40;
	}
	else if (delta == 2)
	{
		tmp_counter_40 = bitmapper_index_params.high_occ_table[((line >> 16) << 1) + 1];

		tmp_counter_40 = tmp_counter_40 + ((bitmapper_index_params.bwt[high_occ_iter] >> (32 - (low_occ_iter << 5))) & 65535);

		return tmp_counter_40;

	}
	else if (delta == 0)
	{
		bitmapper_bs_iter tmp_iter = ((line >> 16) << 1);

		tmp_counter_40 = bitmapper_index_params.high_occ_table[tmp_iter]
			+ bitmapper_index_params.high_occ_table[tmp_iter + 1];

		tmp_iter = (bitmapper_index_params.bwt[high_occ_iter] >> (32 - (low_occ_iter << 5)));

		tmp_counter_40 = tmp_counter_40 + ((tmp_iter >> 16) & 65535) + (tmp_iter & 65535);

		return ((line >> 6) << 6) - tmp_counter_40;

	}
	**/


	if (delta != 0)
	{
		delta = delta >> 1;
		
		tmp_counter_40 = bitmapper_index_params.high_occ_table[((line >> 16) << 1) + delta];

		tmp_counter_40 = tmp_counter_40 + 
			((bitmapper_index_params.bwt[high_occ_iter] >> ((48 - (delta << 4)) - (low_occ_iter << 5))) & 65535);

		return tmp_counter_40;
	}
	else
	{
		bitmapper_bs_iter tmp_iter = ((line >> 16) << 1);

		tmp_counter_40 = bitmapper_index_params.high_occ_table[tmp_iter]
			+ bitmapper_index_params.high_occ_table[tmp_iter + 1];

		tmp_iter = (bitmapper_index_params.bwt[high_occ_iter] >> (32 - (low_occ_iter << 5)));

		tmp_counter_40 = tmp_counter_40 + ((tmp_iter >> 16) & 65535) + (tmp_iter & 65535);

		return ((line >> 6) << 6) - tmp_counter_40;

	}

}






///注意字符0的high_occ_value都是没减过的值，而occ_value是减过的值
inline bitmapper_bs_iter get_occ_value_sp(bitmapper_bs_iter line, int delta,
	bitmapper_bs_iter* actually_line, bitmapper_bs_iter* high_occ_actually_line,
	bitmapper_bs_iter* occ_value, bitmapper_bs_iter* high_occ_value)
{

	bitmapper_bs_iter high_occ_iter = ((line >> 7) * 160) >> 5;
	bitmapper_bs_iter low_occ_iter = (line & 127) >> 6;
	(*actually_line) = high_occ_iter + 1 + (low_occ_iter << 1);
	(*high_occ_actually_line) = high_occ_iter;



	bitmapper_bs_iter tmp_counter_40 = (bwt_string_type)0;

	if (delta != 0)
	{
		delta = delta >> 1;

		tmp_counter_40 = bitmapper_index_params.high_occ_table[((line >> 16) << 1) + delta];
		(*high_occ_value) = tmp_counter_40;



		tmp_counter_40 = tmp_counter_40 +
			((bitmapper_index_params.bwt[high_occ_iter] >> ((48 - (delta << 4)) - (low_occ_iter << 5))) & 65535);
		(*occ_value) = tmp_counter_40;

		return tmp_counter_40;
	}
	else
	{
		bitmapper_bs_iter tmp_iter = ((line >> 16) << 1);

		tmp_counter_40 = bitmapper_index_params.high_occ_table[tmp_iter]
			+ bitmapper_index_params.high_occ_table[tmp_iter + 1];
		(*high_occ_value) = tmp_counter_40;


		tmp_iter = (bitmapper_index_params.bwt[high_occ_iter] >> (32 - (low_occ_iter << 5)));
		tmp_counter_40 = tmp_counter_40 + ((tmp_iter >> 16) & 65535) + (tmp_iter & 65535);
		(*occ_value) = ((line >> 6) << 6) - tmp_counter_40;

		return tmp_counter_40;
		///return ((line >> 6) << 6) - tmp_counter_40;

	}

}




///只需要算字符1和字符2就行了，字符0的结果可以整体减一下就行
inline void get_occ_value_sp_all(bitmapper_bs_iter line, 
	bitmapper_bs_iter* actually_line, bitmapper_bs_iter* high_occ_actually_line,
	bitmapper_bs_iter* occ_value_1, bitmapper_bs_iter* occ_value_2, 
	bitmapper_bs_iter* high_occ_value_1, bitmapper_bs_iter* high_occ_value_2)
{

	bitmapper_bs_iter high_occ_iter = ((line >> 7) * 160) >> 5;
	bitmapper_bs_iter low_occ_iter = (line & 127) >> 6;
	(*actually_line) = high_occ_iter + 1 + (low_occ_iter << 1);
	(*high_occ_actually_line) = high_occ_iter;

	bitmapper_bs_iter tmp;

	bitmapper_bs_iter tmp_counter_1 = (bwt_string_type)0;
	bitmapper_bs_iter tmp_counter_2 = (bwt_string_type)0;

	tmp = ((line >> 16) << 1);

	tmp_counter_1 = bitmapper_index_params.high_occ_table[tmp];
	tmp_counter_2 = bitmapper_index_params.high_occ_table[tmp + 1];

	(*high_occ_value_1) = tmp_counter_1;
	(*high_occ_value_2) = tmp_counter_2;

	tmp = (low_occ_iter << 5);

	tmp_counter_1 = tmp_counter_1 +
		((bitmapper_index_params.bwt[high_occ_iter] >> (48 - tmp)) & 65535);

	tmp_counter_2 = tmp_counter_2 +
		((bitmapper_index_params.bwt[high_occ_iter] >> (32 - tmp)) & 65535);

	(*occ_value_1) = tmp_counter_1;
	(*occ_value_2) = tmp_counter_2;
}




///只需要算字符1和字符2就行了，字符0的结果可以整体减一下就行
inline void get_occ_value_ep_all(bitmapper_bs_iter line,
	bitmapper_bs_iter actually_line_sp,
	bitmapper_bs_iter high_occ_actually_line_sp,
	bitmapper_bs_iter occ_value_sp_1,
	bitmapper_bs_iter occ_value_sp_2,
	bitmapper_bs_iter high_occ_value_sp_1,
	bitmapper_bs_iter high_occ_value_sp_2,
	bitmapper_bs_iter* actually_line,
	bitmapper_bs_iter* occ_value_1, 
	bitmapper_bs_iter* occ_value_2)
{



	bitmapper_bs_iter high_occ_iter = ((line >> 7) * 160) >> 5;
	bitmapper_bs_iter low_occ_iter = (line & 127) >> 6;
	(*actually_line) = high_occ_iter + 1 + (low_occ_iter << 1);
	bitmapper_bs_iter tmp;





	///这种情况下sp和ep在同一个block内，直接返回就好了
	if ((*actually_line) == actually_line_sp)
	{
		(*occ_value_1) = occ_value_sp_1;
		(*occ_value_2) = occ_value_sp_2;
	}
	else if (high_occ_iter == high_occ_actually_line_sp) ///这种情况下sp和ep在同一个super block内，只需要计算小block的occ value
	{

		(*occ_value_1) = high_occ_value_sp_1;
		(*occ_value_2) = high_occ_value_sp_2;


		tmp = (low_occ_iter << 5);

		(*occ_value_1) = (*occ_value_1) +
			((bitmapper_index_params.bwt[high_occ_iter] >> (48 - tmp)) & 65535);

		(*occ_value_2) = (*occ_value_2) +
			((bitmapper_index_params.bwt[high_occ_iter] >> (32 - tmp)) & 65535);



	}
	else    ///这种情况下sp和ep也不在同一个super block内
	{

		tmp = ((line >> 16) << 1);

		(*occ_value_1) = bitmapper_index_params.high_occ_table[tmp];
		(*occ_value_2) = bitmapper_index_params.high_occ_table[tmp + 1];


		tmp = (low_occ_iter << 5);

		(*occ_value_1) = (*occ_value_1) +
			((bitmapper_index_params.bwt[high_occ_iter] >> (48 - tmp)) & 65535);

		(*occ_value_2) = (*occ_value_2) +
			((bitmapper_index_params.bwt[high_occ_iter] >> (32 - tmp)) & 65535);

	}



}






///注意字符0的high_occ_value_sp都是没减过的值，而occ_value_sp是减过的值
inline bitmapper_bs_iter get_occ_value_ep(bitmapper_bs_iter line,
	int delta,
	bitmapper_bs_iter actually_line_sp,
	bitmapper_bs_iter high_occ_actually_line_sp,
	bitmapper_bs_iter occ_value_sp,
	bitmapper_bs_iter high_occ_value_sp,
	bitmapper_bs_iter* actually_line)
{



	bitmapper_bs_iter high_occ_iter = ((line >> 7) * 160) >> 5;
	bitmapper_bs_iter low_occ_iter = (line & 127) >> 6;
	(*actually_line) = high_occ_iter + 1 + (low_occ_iter << 1);
	bitmapper_bs_iter tmp_counter_40 = (bwt_string_type)0;





	///这种情况下sp和ep在同一个block内，直接返回就好了
	if ((*actually_line) == actually_line_sp)
	{
		return occ_value_sp;
	}
	else if (high_occ_iter == high_occ_actually_line_sp) ///这种情况下sp和ep在同一个super block内，只需要计算小block的occ value
	{


		tmp_counter_40 = high_occ_value_sp;

		if (delta != 0)
		{
			delta = delta >> 1;

			///tmp_counter_40 = bitmapper_index_params.high_occ_table[((line >> 16) << 1) + delta];

			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter] >> ((48 - (delta << 4)) - (low_occ_iter << 5))) & 65535);

			return tmp_counter_40;
		}
		else
		{
			bitmapper_bs_iter tmp_iter;

			tmp_iter = (bitmapper_index_params.bwt[high_occ_iter] >> (32 - (low_occ_iter << 5)));

			tmp_counter_40 = tmp_counter_40 + ((tmp_iter >> 16) & 65535) + (tmp_iter & 65535);

			return ((line >> 6) << 6) - tmp_counter_40;

		}

	}
	else    ///这种情况下sp和ep也不在同一个super block内
	{

		if (delta != 0)
		{
			delta = delta >> 1;

			tmp_counter_40 = bitmapper_index_params.high_occ_table[((line >> 16) << 1) + delta];

			tmp_counter_40 = tmp_counter_40 +
				((bitmapper_index_params.bwt[high_occ_iter] >> ((48 - (delta << 4)) - (low_occ_iter << 5))) & 65535);

			return tmp_counter_40;
		}
		else
		{
			bitmapper_bs_iter tmp_iter = ((line >> 16) << 1);

			tmp_counter_40 = bitmapper_index_params.high_occ_table[tmp_iter]
				+ bitmapper_index_params.high_occ_table[tmp_iter + 1];

			tmp_iter = (bitmapper_index_params.bwt[high_occ_iter] >> (32 - (low_occ_iter << 5)));

			tmp_counter_40 = tmp_counter_40 + ((tmp_iter >> 16) & 65535) + (tmp_iter & 65535);

			return ((line >> 6) << 6) - tmp_counter_40;

		}
	}



}













/**
inline bitmapper_bs_iter popcount_hash_table(bitmapper_bs_iter value)
{
	bitmapper_bs_iter count = 0;

	count += bitmapper_index_params.hash_count[value & bitmapper_index_params.mode_16];
	count += bitmapper_index_params.hash_count[(value >> 16) & bitmapper_index_params.mode_16];
	count += bitmapper_index_params.hash_count[(value >> 32) & bitmapper_index_params.mode_16];
	count += bitmapper_index_params.hash_count[(value >> 48) & bitmapper_index_params.mode_16];

	return count;
}
**/

inline bitmapper_bs_iter find_occ_fm_index(bitmapper_bs_iter line, int delta, bwt_string_type *bwt, high_occ_table_type* high_occ_table)
{

	if (line > bitmapper_index_params.shapline)
		line--;

	bwt_string_type ans = bitmapper_index_params.nacgt[delta];

	bitmapper_bs_iter actually_line;

	/**
	fprintf(stderr, "*******************\nline=%llu, delta=%llu\n", line, delta);
	**/

	ans = ans + get_occ_value(line, delta, &actually_line);

	/**
	fprintf(stderr, "actually_line=%llu, get_occ_value=%llu\n", actually_line, ans - bitmapper_index_params.nacgt[delta]);

	dec_bit_bwt(bwt[actually_line]);
	dec_bit_bwt(bwt[actually_line + 1]);
	**/


	///__builtin_prefetch(&(bwt[actually_line]), 0, 1);


	bitmapper_bs_iter need_line = (line & 63);

	if (unlikely(need_line == 0))
		return ans;

    bwt_string_type P_A;


	/**
	///1代表是字符T
	if (delta == 1)
	{
        P_A = (bwt[actually_line] >> (64 - need_line));
        ans = ans + __builtin_popcountll(P_A);
		///ans = ans + popcount_hash_table(P_A);

		///bitmapper_index_params.c_1_times++;


	}
	else if (delta == 2)  ///2代表是字符G
	{
        P_A = (bwt[actually_line + 1] >> (64 - need_line));
        ans = ans + __builtin_popcountll(P_A);
		///ans = ans + popcount_hash_table(P_A);

		///bitmapper_index_params.c_2_times++;
	}
	else if (delta == 0)  ///0代表是字符A
	{
        P_A = bwt[actually_line] | bwt[actually_line + 1];
        P_A = ~P_A;
        P_A = (P_A >> (64 - need_line));
        ans = ans + __builtin_popcountll(P_A);
		///ans = ans + popcount_hash_table(P_A);

		///bitmapper_index_params.c_0_times++;
	}
	**/




	if (delta != 0)
	{
		P_A = (bwt[actually_line + (delta >> 1)] >> (64 - need_line));
		ans = ans + __builtin_popcountll(P_A);
		///ans = ans + popcount_hash_table(P_A);


	}
	else  ///0代表是字符A
	{
		P_A = bwt[actually_line] | bwt[actually_line + 1];
		P_A = ~P_A;
		P_A = (P_A >> (64 - need_line));
		ans = ans + __builtin_popcountll(P_A);
		///ans = ans + popcount_hash_table(P_A);
	}


	///__builtin_prefetch(&(bwt[ans]), 0, 1);

    return ans;

}







inline void find_occ_fm_index_combine
    (bitmapper_bs_iter sp, bitmapper_bs_iter ep,
     bitmapper_bs_iter* new_sp, bitmapper_bs_iter* new_ep,
     int delta, bwt_string_type *bwt)
{




    /*************************************calculate sp******************************************************/

	if (sp > bitmapper_index_params.shapline)
		sp--;


	bwt_string_type ans_sp = bitmapper_index_params.nacgt[delta];

	bitmapper_bs_iter actually_line_sp;
	bitmapper_bs_iter high_occ_actually_line_sp;
	bitmapper_bs_iter occ_value_sp;
	bitmapper_bs_iter high_occ_value_sp;


	///注意字符0的high_occ_value都是没减过的值，而occ_value是减过的值
    get_occ_value_sp(sp, delta, &actually_line_sp, &high_occ_actually_line_sp, &occ_value_sp, &high_occ_value_sp);


    ans_sp = ans_sp + occ_value_sp;


	///__builtin_prefetch(&(bwt[actually_line]), 0, 1);


	bitmapper_bs_iter need_sp = (sp & 63);

	bwt_string_type P_A;

	if (likely(need_sp != 0))
    {
        if (delta != 0)
        {
            P_A = (bwt[actually_line_sp + (delta >> 1)] >> (64 - need_sp));
            ans_sp = ans_sp + __builtin_popcountll(P_A);
            ///ans = ans + popcount_hash_table(P_A);


        }
        else  ///0代表是字符A
        {
            P_A = bwt[actually_line_sp] | bwt[actually_line_sp + 1];
            P_A = ~P_A;
            P_A = (P_A >> (64 - need_sp));
            ans_sp = ans_sp + __builtin_popcountll(P_A);
            ///ans = ans + popcount_hash_table(P_A);
        }

    }


    (*new_sp) = ans_sp;

    /*************************************calculate sp******************************************************/




    /*************************************calculate ep******************************************************/






    if (ep > bitmapper_index_params.shapline)
		ep--;


    bwt_string_type ans_ep = bitmapper_index_params.nacgt[delta];

    bitmapper_bs_iter actually_line_ep;

	///注意字符0的high_occ_value_sp都是没减过的值，而occ_value_sp是减过的值
	ans_ep = ans_ep + get_occ_value_ep(ep, delta, actually_line_sp, high_occ_actually_line_sp, occ_value_sp, high_occ_value_sp,
                                          &actually_line_ep);

    bitmapper_bs_iter need_ep = (ep & 63);


    if (likely(need_ep != 0))
    {
        if (delta != 0)
        {
            P_A = (bwt[actually_line_ep + (delta >> 1)] >> (64 - need_ep));
            ans_ep = ans_ep + __builtin_popcountll(P_A);
            ///ans = ans + popcount_hash_table(P_A);


        }
        else  ///0代表是字符A
        {
            P_A = bwt[actually_line_ep] | bwt[actually_line_ep + 1];
            P_A = ~P_A;
            P_A = (P_A >> (64 - need_ep));
            ans_ep = ans_ep + __builtin_popcountll(P_A);
            ///ans = ans + popcount_hash_table(P_A);
        }

    }


    (*new_ep) = ans_ep;



    /*************************************calculate ep******************************************************/




	///__builtin_prefetch(&(bwt[ans]), 0, 1);

    ///return ans_sp;

}



inline void find_occ_fm_index_combine_back
(bitmapper_bs_iter sp, bitmapper_bs_iter ep,
bitmapper_bs_iter* new_sp, bitmapper_bs_iter* new_ep,
int delta, bwt_string_type *bwt)
{




	/*************************************calculate sp******************************************************/

	if (sp > bitmapper_index_params.shapline)
		sp--;


	bwt_string_type ans_sp = bitmapper_index_params.nacgt[delta];

	bitmapper_bs_iter actually_line_sp;
	bitmapper_bs_iter high_occ_actually_line_sp;
	bitmapper_bs_iter occ_value_sp;
	bitmapper_bs_iter high_occ_value_sp;


	///注意字符0的high_occ_value都是没减过的值，而occ_value是减过的值
	get_occ_value_sp(sp, delta, &actually_line_sp, &high_occ_actually_line_sp, &occ_value_sp, &high_occ_value_sp);


	ans_sp = ans_sp + occ_value_sp;


	///__builtin_prefetch(&(bwt[actually_line]), 0, 1);


	bitmapper_bs_iter need_sp = (sp & 63);

	bwt_string_type P_A;

	if (likely(need_sp != 0))
	{
		if (delta != 0)
		{
			P_A = (bwt[actually_line_sp + (delta >> 1)] >> (64 - need_sp));
			ans_sp = ans_sp + __builtin_popcountll(P_A);
			///ans = ans + popcount_hash_table(P_A);


		}
		else  ///0代表是字符A
		{
			P_A = bwt[actually_line_sp] | bwt[actually_line_sp + 1];
			P_A = ~P_A;
			P_A = (P_A >> (64 - need_sp));
			ans_sp = ans_sp + __builtin_popcountll(P_A);
			///ans = ans + popcount_hash_table(P_A);
		}

	}


	(*new_sp) = ans_sp;

	/*************************************calculate sp******************************************************/




	/*************************************calculate ep******************************************************/






	if (ep > bitmapper_index_params.shapline)
		ep--;


	bwt_string_type ans_ep = bitmapper_index_params.nacgt[delta];

	bitmapper_bs_iter actually_line_ep;

	///注意字符0的high_occ_value_sp都是没减过的值，而occ_value_sp是减过的值
	ans_ep = ans_ep + get_occ_value_ep(ep, delta, actually_line_sp, high_occ_actually_line_sp, occ_value_sp, high_occ_value_sp,
		&actually_line_ep);

	bitmapper_bs_iter need_ep = (ep & 63);


	if (likely(need_ep != 0))
	{
		if (delta != 0)
		{
			P_A = (bwt[actually_line_ep + (delta >> 1)] >> (64 - need_ep));
			ans_ep = ans_ep + __builtin_popcountll(P_A);
			///ans = ans + popcount_hash_table(P_A);


		}
		else  ///0代表是字符A
		{
			P_A = bwt[actually_line_ep] | bwt[actually_line_ep + 1];
			P_A = ~P_A;
			P_A = (P_A >> (64 - need_ep));
			ans_ep = ans_ep + __builtin_popcountll(P_A);
			///ans = ans + popcount_hash_table(P_A);
		}

	}


	(*new_ep) = ans_ep;



	/*************************************calculate ep******************************************************/




	///__builtin_prefetch(&(bwt[ans]), 0, 1);

	///return ans_sp;

}







inline bitmapper_bs_iter count(char* pattern, bitmapper_bs_iter length,
	bitmapper_bs_iter* sp, bitmapper_bs_iter* ep, bitmapper_bs_iter* sp1, bitmapper_bs_iter* ep1)
{

	long long j = length - 1;

	bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];


	///if (bitmapper_index_params.delta > 3)
	///if (unlikely(bitmapper_index_params.delta > 2))
	if (bitmapper_index_params.delta > 2)
	{
		return 0;
	}

	bitmapper_bs_iter top = bitmapper_index_params.nacgt[bitmapper_index_params.delta];
	bitmapper_bs_iter bot = bitmapper_index_params.nacgt[bitmapper_index_params.delta + 1];


	for (j = j - 1; j >= 1; j--)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

		///if (bitmapper_index_params.delta > 3)
		///if (unlikely(bitmapper_index_params.delta > 2))
		if (bitmapper_index_params.delta > 2)
		{
			return 0;
		}


		top = find_occ_fm_index(top, bitmapper_index_params.delta,
			bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);
		bot = find_occ_fm_index(bot, bitmapper_index_params.delta,
			bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);

		if (bot <= top)
		{
			break;
		}

	}

	(*sp1) = top;
	(*ep1) = bot;


	for (; j >= 0; j--)
	{

		if (bot <= top)
		{
			break;
		}


		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];


		///if (bitmapper_index_params.delta > 3)
		///if (unlikely(bitmapper_index_params.delta > 2))
		if (bitmapper_index_params.delta > 2)
		{
			return 0;
		}



		top = find_occ_fm_index(top, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);
		bot = find_occ_fm_index(bot, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);


	}










	(*sp) = top;
	(*ep) = bot;


	if (bot <= top)
	{
		return 0;
	}
	else
	{


		return bot - top;
	}



}














inline bitmapper_bs_iter count_hash_table(char* pattern, bitmapper_bs_iter length,
	bitmapper_bs_iter* sp, bitmapper_bs_iter* ep, bitmapper_bs_iter* sp1, bitmapper_bs_iter* ep1)
{

	if (length < 17)
	{
		return 0;
	}

	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;



	////这个j还必须得是long long...
	long long j = length - 16;

	bitmapper_bs_iter hash_value = get_3_letter_hash_value(pattern + j, 16);


	if (hash_value == ((bitmapper_bs_iter)-1))
	{
		return 0;
	}

	query_16_mer_hash_table(hash_value, &top, &bot);

	if (bot <= top)
	{
		return 0;
	}










	for (j = j - 1; j >= 1; j--)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

		if (bitmapper_index_params.delta > 2)
		{
			return 0;
		}


		find_occ_fm_index_combine(top, bot, &top, &bot, bitmapper_index_params.delta, bitmapper_index_params.bwt);


		if (bot <= top)
		{
			break;
		}

	}

	(*sp1) = top;
	(*ep1) = bot;


	for (; j >= 0; j--)
	{

		if (bot <= top)
		{
			break;
		}


		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

		if (bitmapper_index_params.delta > 2)
		{
			return 0;
		}


        find_occ_fm_index_combine(top, bot, &top, &bot, bitmapper_index_params.delta,bitmapper_index_params.bwt);

	}


	(*sp) = top;
	(*ep) = bot;


	if (bot <= top)
	{
		return 0;
	}
	else
	{
		return bot - top;
	}




}




inline bitmapper_bs_iter count_backward_as_much(char* pattern, bitmapper_bs_iter length,
	bitmapper_bs_iter* sp, bitmapper_bs_iter* ep, bitmapper_bs_iter* sp1, bitmapper_bs_iter* ep1,
	bitmapper_bs_iter* match_length)
{

	(*match_length) = 0;

	///这个干脆改成18好了...
	if (length < 18)
	{
		return 0;
	}

	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter pre_pre_top;
	bitmapper_bs_iter pre_pre_bot;


	////这个j还必须得是long long...
	long long j = length - 16;

	bitmapper_bs_iter hash_value = get_3_letter_hash_value(pattern + j, 16);


	if (hash_value == ((bitmapper_bs_iter)-1))
	{
		return 0;
	}

	query_16_mer_hash_table(hash_value, &top, &bot);

	if (bot <= top)
	{
		return 0;
	}





	
	pre_top = (bitmapper_bs_iter)-1;
	pre_bot = (bitmapper_bs_iter)-1;
	

	for (j = j - 1; j >= 0; j--)
	{

		pre_pre_top = pre_top;
		pre_pre_bot = pre_bot;
		pre_top = top;
		pre_bot = bot;


		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

		if (bitmapper_index_params.delta > 2)
		{
			(*match_length) = length - j -1;
			bot = top;
			break;
		}


		find_occ_fm_index_combine(top, bot, &top, &bot, bitmapper_index_params.delta, bitmapper_index_params.bwt);


		if (bot <= top)
		{
			(*match_length) = length - j - 1;
			break;
		}

		

	}



	///这个说明遇到了不能再搜的情况
	if (bot <= top)
	{
		///在这种情况下, 当前的有效区间实际上不是top和bot
		///而是倒退一位的pre_top和pre_bot
		(*sp) = pre_top;
		(*ep) = pre_bot;

		///而有效的前驱区间实际上是pre_pre_top和pre_pre_bot
		(*sp1) = pre_pre_top;
		(*ep1) = pre_pre_bot;


		///请注意, 当(*match_length)是16时，这个pre_pre_top和pre_pre_bot均是-1
		///从逻辑上实际上有问题
		///不过我们一般对长为16的种子直接不要

		return (*ep) - (*sp);
	}
	else   ///这里说明搜到read的末尾了，仍然匹配
	{
		///首先匹配了这么多长度
		(*match_length) = length;


		///在这种情况下, 当前的有效区间就是top和bot
		(*sp) = top;
		(*ep) = bot;


		///有效的前驱区间是pre_top和pre_bot
		(*sp1) = pre_top;
		(*ep1) = pre_bot;

		return (*ep) - (*sp);
	}
}





inline bitmapper_bs_iter count_backward_as_much_1_terminate(char* pattern, bitmapper_bs_iter length,
	bitmapper_bs_iter* sp, bitmapper_bs_iter* ep, bitmapper_bs_iter* sp1, bitmapper_bs_iter* ep1,
	bitmapper_bs_iter* match_length)
{

	(*match_length) = 0;

	///这个干脆改成18好了...
	if (length < 18)
	{
		return 0;
	}

	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;
	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;
	bitmapper_bs_iter pre_pre_top;
	bitmapper_bs_iter pre_pre_bot;


	////这个j还必须得是long long...
	long long j = length - 16;

	bitmapper_bs_iter hash_value = get_3_letter_hash_value(pattern + j, 16);


	if (hash_value == ((bitmapper_bs_iter)-1))
	{
		return 0;
	}

	query_16_mer_hash_table(hash_value, &top, &bot);

	if (bot <= top)
	{
		return 0;
	}






	pre_top = (bitmapper_bs_iter)-1;
	pre_bot = (bitmapper_bs_iter)-1;


	for (j = j - 1; j >= 0; j--)
	{

		

		pre_pre_top = pre_top;
		pre_pre_bot = pre_bot;
		pre_top = top;
		pre_bot = bot;


		if (bot - top == 1)
		{
			(*match_length) = length - j - 1;
			break;
		}



		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

		if (bitmapper_index_params.delta > 2)
		{
			(*match_length) = length - j - 1;
			bot = top;
			break;
		}


		find_occ_fm_index_combine(top, bot, &top, &bot, bitmapper_index_params.delta, bitmapper_index_params.bwt);


		if (bot <= top)
		{
			(*match_length) = length - j - 1;
			break;
		}



	}



	///这个说明遇到了不能再搜的情况
	if (bot <= top)
	{
		///在这种情况下, 当前的有效区间实际上不是top和bot
		///而是倒退一位的pre_top和pre_bot
		(*sp) = pre_top;
		(*ep) = pre_bot;

		///而有效的前驱区间实际上是pre_pre_top和pre_pre_bot
		(*sp1) = pre_pre_top;
		(*ep1) = pre_pre_bot;


		///请注意, 当(*match_length)是16时，这个pre_pre_top和pre_pre_bot均是-1
		///从逻辑上实际上有问题
		///不过我们一般对长为16的种子直接不要

		return (*ep) - (*sp);
	}
	else   ///这里说明搜到read的末尾了，仍然匹配
	{
		///首先匹配了这么多长度
		///(*match_length) = length;
		(*match_length) = length - j - 1;

		///在这种情况下, 当前的有效区间就是top和bot
		(*sp) = top;
		(*ep) = bot;


		///有效的前驱区间是pre_top和pre_bot
		(*sp1) = pre_top;
		(*ep1) = pre_bot;

		return (*ep) - (*sp);
	}
}








inline bitmapper_bs_iter count_forward_as_mush(char* pattern, bitmapper_bs_iter length,
	bitmapper_bs_iter* sp, bitmapper_bs_iter* ep, bitmapper_bs_iter* sp1, bitmapper_bs_iter* ep1, 
	bitmapper_bs_iter* match_length)
{

	(*match_length) = 0;

	///这里干脆改成大于18好了....
	if (length < 18)
	{
		return 0;
	}

	bitmapper_bs_iter top;
	bitmapper_bs_iter bot;

	bitmapper_bs_iter pre_top;
	bitmapper_bs_iter pre_bot;

	bitmapper_bs_iter pre_pre_top;
	bitmapper_bs_iter pre_pre_bot;



	////这个j还必须得是long long...
	///long long j = length - 16;
	long long j = 0;

	bitmapper_bs_iter hash_value = get_3_letter_hash_value(pattern + j, 16);
	///N_site代表的是第一次出现N的位置
	///bitmapper_bs_iter hash_value = get_3_letter_hash_value_N_site(pattern + j, 16, N_site);

	if (hash_value == ((bitmapper_bs_iter)-1))
	{
		return 0;
	}

	query_16_mer_hash_table(hash_value, &top, &bot);

	if (bot <= top)
	{
		return 0;
	}

	pre_top = (bitmapper_bs_iter)-1;
	pre_bot = (bitmapper_bs_iter)-1;

	for (j = 16; j < length; j++)
	{

		pre_pre_top = pre_top;
		pre_pre_bot = pre_bot;
		pre_top = top;
		pre_bot = bot;



		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

		if (bitmapper_index_params.delta > 2)
		{
			(*match_length) = j;
			bot = top;
			break;
		}


		find_occ_fm_index_combine(top, bot, &top, &bot, bitmapper_index_params.delta, bitmapper_index_params.bwt);


		if (bot <= top)
		{
			(*match_length) = j;
			break;
		}

	}




	///这个说明遇到了不能再搜的情况
	if (bot <= top)
	{
		///在这种情况下, 当前的有效区间实际上不是top和bot
		///而是倒退一位的pre_top和pre_bot
		(*sp) = pre_top;
		(*ep) = pre_bot;

		///而有效的前驱区间实际上是pre_pre_top和pre_pre_bot
		(*sp1) = pre_pre_top;
		(*ep1) = pre_pre_bot;


		///请注意, 当(*match_length)是16时，这个pre_pre_top和pre_pre_bot均是-1
		///从逻辑上实际上有问题
		///不过我们一般对长为16的种子直接不要

		return (*ep) - (*sp);
	}
	else   ///这里说明搜到read的末尾了，仍然匹配
	{
		///首先匹配了这么多长度
		(*match_length) = length;


		///在这种情况下, 当前的有效区间就是top和bot
		(*sp) = top;
		(*ep) = bot;


		///有效的前驱区间是pre_top和pre_bot
		(*sp1) = pre_top;
		(*ep1) = pre_bot;

		return (*ep) - (*sp);
	}
}
















inline unsigned int count_less_than_4(char* pattern, unsigned int length,
	unsigned int* sp, unsigned int* ep)
{


	long long j = length - 1;

	bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

	if (bitmapper_index_params.delta > 3)
	{
		return 0;
	}

	unsigned int top = bitmapper_index_params.nacgt[bitmapper_index_params.delta];
	unsigned int bot = bitmapper_index_params.nacgt[bitmapper_index_params.delta + 1];



	j--;



	for (; j >= 0; j--)
	{
		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];


		if (bitmapper_index_params.delta > 3)
		{
			return 0;
		}



		top = find_occ_fm_index(top, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);
		bot = find_occ_fm_index(bot, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);

		if (bot <= top)
		{
			return 0;
		}

	}




	(*sp) = top;
	(*ep) = bot;


	return bot - top;

}



void debug_information();


inline bwt_string_type access_bwt_delta(bitmapper_bs_iter l)
{
	///这一溜是取BWT的,必须考虑$
	if (l > bitmapper_index_params.shapline)
		l--;


	bitmapper_bs_iter actually_line = (((l >> 7) * 160) >> 5) + 1 + (((l & 127) >> 6) << 1);

	bitmapper_bs_iter need_line = (l & 63);


	bwt_string_type delta;


	delta = (bitmapper_index_params.bwt[actually_line] << need_line) & bitmapper_index_params.mode_bwt_delta;

	if (delta)
	{
		return 1;
	}
	else
	{
		delta = (bitmapper_index_params.bwt[actually_line + 1] << need_line) & bitmapper_index_params.mode_bwt_delta;

		if (delta)
		{
			return 2;
		}
		else
		{
			return 0;
		}
	}
}

inline bitmapper_bs_iter bwt_get_sa_restrict_steps_more_than_3
(SA_flag_string_type* SA_flag,
bitmapper_bs_iter line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa,
bitmapper_bs_iter need_step,
bitmapper_bs_iter* accessed_sa)///unsigned int* accessed_sa)
{
	bitmapper_bs_iter l = line;
	bwt_string_type delta;
	bitmapper_bs_iter i = 0;


	if (line == bitmapper_index_params.shapline)
	{
		(*accessed_sa) = 0;
		return 1;
	}


	bitmapper_bs_iter actually_line;
	bitmapper_bs_iter last;

	///注意对SA flag而言, $并不影响其取法
	actually_line = ((l >> 8) * bitmapper_index_params.acctuall_SA_flag_gap) >> 6;
	last = l & 255;
	actually_line = actually_line + (last >> 6) + 1; ///注意这里要+1



	while (((SA_flag[actually_line] << (last & 63))&bitmapper_index_params.mode_SA_flag) == 0)
	{

		if (i >= need_step)
		{
			return 0;
		}


		delta = access_bwt_delta(l);

		l = find_occ_fm_index(l, delta, bwt, high_occ_table);


		i++;

		if (l == bitmapper_index_params.shapline)
		{
			(*accessed_sa) = i;

			return 1;
		}


		///注意对SA flag而言, $并不影响其取法
		actually_line = ((l >> 8) * bitmapper_index_params.acctuall_SA_flag_gap) >> 6;
		last = l & 255;
		actually_line = actually_line + (last >> 6) + 1; ///注意这里要+1

	}




	actually_line = actually_line - (last >> 6) - 1; ///注意这里要-1
	bitmapper_bs_iter ans = SA_flag[actually_line];




	///这里SA要改
	if (last != 0)
	{
		bitmapper_bs_iter j = 0;


		actually_line++;
		SA_flag_string_type tmp_SA_pop_count = SA_flag[actually_line];



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

	(*accessed_sa) = (bitmapper_bs_iter)(sa[ans] & bitmapper_index_params.SA_header_mode)
		*bitmapper_index_params.compress_sa + i;

	return 1;
}




inline void locate_one_position(
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter* tmp_SA_length)
{


		///这里SA要改
	if (bwt_get_sa_restrict_steps_more_than_3(bitmapper_index_params.SA_flag, sp, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table, bitmapper_index_params.sa,
		bitmapper_index_params.compress_sa - 1, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			//locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
		    locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)];
			(*tmp_SA_length)++;
		}


}


inline void locate_one_position_direct(
	bitmapper_bs_iter* locates,///unsigned int* locates,
	bitmapper_bs_iter sp,
	bitmapper_bs_iter* tmp_SA_length,
	bitmapper_bs_iter total_SA_length,
	bitmapper_bs_iter seed_length, bitmapper_bs_iter seed_offset)
{


	///这里SA要改
	if (bwt_get_sa_restrict_steps_more_than_3(bitmapper_index_params.SA_flag, sp, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table, bitmapper_index_params.sa,
		bitmapper_index_params.compress_sa - 1, (&(locates[(*tmp_SA_length)]))) == 1)
	{
		locates[(*tmp_SA_length)] = total_SA_length - locates[(*tmp_SA_length)] - seed_length - seed_offset;

		(*tmp_SA_length)++;
	}


}

#endif

