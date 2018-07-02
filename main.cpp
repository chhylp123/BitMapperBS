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
#include<stdint.h>
#include<ctype.h>
#include "bwt.h"



void convert_2bits_to_chars(long long input, char* output, unsigned int length_of_chars)
{
	unsigned int i = 0, j = 0;
	char ACGT_C[4];

	ACGT_C[0] = 'A';
	ACGT_C[1] = 'C';
	ACGT_C[2] = 'G';
	ACGT_C[3] = 'T';


	output[length_of_chars] = '\0';

	j = length_of_chars - 1;

	for (i = 0; i <length_of_chars; i++)
	{
		output[j] = ACGT_C[input % 4];
		input = input / 4;
		j--;
	}
}



void dec_bit(bwt_string_type input)
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


void output_bit_data(uint64_t input_value)
{
	int bit_values[64];

	uint64_t mode_1 = (uint64_t)1;

	for (int i = 0; i < 64; i++)
	{
		bit_values[i] = input_value&mode_1;
		input_value = input_value >> 1;
	}

	for (int i = 63; i >= 0; i--)
	{
		fprintf(stdout, "%u", bit_values[i]);
	}

	fprintf(stdout, "\n");

}


void Generate_patterns() {

	char filename[100];



	fprintf(stdout, "Please input the text file name (patterns will be randomly generated from file):\n");
	scanf("%s", filename);


	fprintf(stdout, "\n\n\n*********************Warning*********************\n\n");



















	fprintf(stdout,
		"We strongly recommend that input text (genome) should only consist of the characters in {a, c, g, t, A, C, G, T}!\n");
	fprintf(stdout,
		"If not, the generated patterns may include the characters which do not belong to {a, c, g, t, A, C, G, T}. \n");


	fprintf(stdout,
		"You could preprocess input text (genome) using the program named 'preprocess'. \n");

	fprintf(stdout, "\n*********************Warning*********************\n\n\n");







	FILE *_ih_fp = fopen(filename, "r");

	fseek(_ih_fp, 0, SEEK_END);
	unsigned int n = ftell(_ih_fp);

	///fprintf(stdout, "r_n=%u\n", n);

	fseek(_ih_fp, 0, SEEK_SET);


	unsigned char* text = (unsigned char *)malloc((n)*(sizeof(unsigned char)));

	fread(text, sizeof(unsigned char), n, _ih_fp);


	fclose(_ih_fp);

	_ih_fp = fopen("patterns.txt", "w");

	unsigned int query_number=0;
	unsigned int query_length = 0;


	fprintf(stdout, "Please input the required pattern length:\n");
	scanf("%u", &query_length);

	fwrite(&query_length, sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "Please input the required pattern nunmber\n");
	scanf("%u", &query_number);

	fwrite(&query_number, sizeof(unsigned int), 1, _ih_fp);

















	unsigned int occur = 0;

	long long tmp_loc = 0;




	occur = 0;





	srandom((unsigned int)time(NULL));

	while (occur<query_number)
	{



		tmp_loc=random();

		if (tmp_loc + query_length<=n)
		{
			occur++;

			fwrite(text + tmp_loc, sizeof(unsigned char), query_length, _ih_fp);

			///fprintf(_ih_fp, "\n");

		}

	}





}



int build_index()
{
	char filename[100], filenames[100], filenameo[100], filenameb[100];
	char *refer;
	char ch;
	bitmapper_bs_iter i, j;

	FILE *f1;






	fprintf(stdout, "\n\n\n*********************Warning*********************\n\n");

	fprintf(stdout,
		"The input text must only consist of the characters in {a, c, g, t, A, C, G, T}!\n");
	fprintf(stdout,
		"If not, you could preprocess it using the program named 'preprocess'. \n");

	fprintf(stdout, "\n*********************Warning*********************\n\n\n");



	fprintf(stdout,
		"Input file name:\n");
	scanf("%s", filename);

	bitmapper_bs_iter compress_build_sa;

	fprintf(stdout, "Please input the sampling distance D of SA (1 < D < 9):\n");
	scanf("%u", &compress_build_sa);

	if (compress_build_sa>=9 || compress_build_sa<=1)
	{
		fprintf(stdout, "Do not input allowed sampling distance D! FMtree will exit ...\n");
		return 1;
	}


	bitmapper_bs_iter text_length = 0;
	f1 = fopen(filename, "r");


	if (NULL == f1)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filename);
		return 1;
	}


	bitmapper_bs_iter number_of_total_characters = 0;

	///while (!feof(f1))
	while (1)
	{
		fscanf(f1, "%c", &ch);

		if (feof(f1))
		{
			break;
		}

		if ((ch == 'A') || (ch == 'C') || (ch == 'G') || (ch == 'T')||
			(ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't'))
			text_length++;
		else
		{

			fprintf(stdout, "%u-th character does not belong to {a, c, g, t, A, C, G, T}. Its  ASCII Code is %u, and it is %c.\n",
				number_of_total_characters, ch, ch);
		}

		number_of_total_characters++;

	}

	printf("text_length=%llu\n", text_length);

	fclose(f1);


	refer = (char *)malloc(sizeof(char)*(text_length + 1));




	f1 = fopen(filename, "r");
	for (i = 0; i <text_length; i++)
	{
		fscanf(f1, "%c", &ch);

		refer[i] = ch;
	}


	indenpendent_creadte_index(text_length, &refer, compress_build_sa, filename);


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


void search_from_bwt_more_than_3()
{
	long long i, j;
	int length_read;
	char* reads;
	bitmapper_bs_iter top, bot, t;
	bitmapper_bs_iter pre_top, pre_bot;


	FILE* _ih_fp = fopen("patterns.txt", "r");

	unsigned int numocc;
	unsigned int query_length = 0;

	fread(&(query_length), sizeof(unsigned int), 1, _ih_fp);
	fread(&(numocc), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "query_length=%u\n", query_length);
	fprintf(stdout, "numocc=%u\n", numocc);


	char* patterns = (char *)malloc((query_length*numocc)*(sizeof(char)));

	fread(patterns, sizeof(char), (query_length*numocc), _ih_fp);


	reads = patterns;

	reads = patterns;

	length_read = query_length;

	unsigned int num_reads = numocc;

	bitmapper_bs_iter* locates;///unsigned int* locates;

	double start = clock();


	bitmapper_bs_iter number_of_locations = 0;
	bitmapper_bs_iter count_number_of_locations = 0;



	bitmapper_bs_iter tmp_SA_length = 0;





	long long number_of_hits = 0;




	struct  timeval  start_timeval;
	struct  timeval  end_timeval;
	unsigned long timer;
	gettimeofday(&start_timeval, NULL);




	for (i = 1; i <= num_reads; i++)
	{

		if (i % 10000 == 0)
		{
			fprintf(stdout, "#i=%llu\n", i);
		}

		if (length_read < 17)
		{
			fprintf(stdout, "ERROR: the length of read is less than 17!!!!\n");
			break;
		}

		/**
		number_of_hits
			= count(reads, length_read, &top, &bot, &pre_top, &pre_bot);
		**/


		number_of_hits
			= count_hash_table(reads, length_read, &top, &bot, &pre_top, &pre_bot);
		


		///fprintf(stderr, "i=%llu, number_of_hits=%llu\n",i, number_of_hits);




		if (number_of_hits==0)
		{
			reads = reads + length_read;
			continue;
		}




		locates = (bitmapper_bs_iter *)malloc(number_of_hits*sizeof(bitmapper_bs_iter));///locates = (unsigned int *)malloc(number_of_hits*sizeof(unsigned int));

		tmp_SA_length = 0;




		///locate(reads, top, bot, pre_top, pre_bot, locates, length_read, &tmp_SA_length);



		/**
		locate_debug(reads, top, bot,
			pre_top, pre_bot, locates, length_read, &tmp_SA_length);
		**/

		

		/**
		qsort(locates, bot - top, sizeof(bitmapper_bs_iter), bitmapper_bs_iter_compareEntrySize);

		for (t = 0; t<bot - top; t++)
		{

			fprintf(stderr, "i=%llu, site=%llu\n", i, locates[t]);
			fflush(stderr);
		}
		
        **/




		free(locates);

		/**
		if (tmp_SA_length != number_of_hits)
		{
			fprintf(stdout, "i=%llu, tmp_SA_length=%llu, number_of_hits=%llu, top=%llu, bot=%llu\n",
				i, tmp_SA_length, number_of_hits, top, bot);
		}
		**/


		number_of_locations = number_of_locations + tmp_SA_length;
		///number_of_locations = number_of_locations + number_of_hits;
		count_number_of_locations = count_number_of_locations + number_of_hits;






		reads = reads + length_read;

	}

	gettimeofday(&end_timeval, NULL);
	timer = 1000000 * (end_timeval.tv_sec - start_timeval.tv_sec) + end_timeval.tv_usec - start_timeval.tv_usec;


	double finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;


	fprintf(stdout, "\n\n\n*********************Result*********************\n");


	fprintf(stdout, "searching Time: %ld microsecond\n", timer);

	fprintf(stdout, "searching Time: %5.3f seconds\n", duration);

	fprintf(stdout, "number of matched locations=%ld\n", number_of_locations);

	fprintf(stdout, "count: number of matched locations=%ld\n", count_number_of_locations);

	


	fprintf(stdout, "************************************************\n");

	fprintf(stdout, "bitmapper_index_params.c_0_times=%llu\n", bitmapper_index_params.c_0_times);
	fprintf(stdout, "bitmapper_index_params.c_1_times=%llu\n", bitmapper_index_params.c_1_times);
	fprintf(stdout, "bitmapper_index_params.c_2_times=%llu\n", bitmapper_index_params.c_2_times);

	///debug_information();


}




/**



void test_hash()
{
	long long i, j;
	int length_read;
	char* reads;
	unsigned int top, bot, t;
	unsigned int pre_top, pre_bot;


	FILE* _ih_fp = fopen("patterns.txt", "r");

	unsigned int numocc;
	unsigned int query_length = 0;

	fread(&(query_length), sizeof(unsigned int), 1, _ih_fp);
	fread(&(numocc), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "query_length=%u\n", query_length);
	fprintf(stdout, "numocc=%u\n", numocc);


	char* patterns = (char *)malloc((query_length*numocc)*(sizeof(char)));

	fread(patterns, sizeof(char), (query_length*numocc), _ih_fp);


	reads = patterns;

	reads = patterns;

	length_read = query_length;

	unsigned int num_reads = numocc;

	unsigned int* locates;

	double start = clock();


	long long number_of_locations = 0;



	unsigned int tmp_SA_length = 0;





	long long number_of_hits = 0;




	struct  timeval  start_timeval;
	struct  timeval  end_timeval;
	unsigned long timer;
	gettimeofday(&start_timeval, NULL);


	int WINDOW_SIZE = 11;
	unsigned int pre_hashTableMaxSize = pow(4, WINDOW_SIZE);
	char* seeds = (char*)malloc(sizeof(char)*(WINDOW_SIZE + 1));
	unsigned int counting, last_bot, last_counting;









	char ch;

	FILE* f1;

	char filename[100];

	fprintf(stdout,
		"Input file name:\n");
	scanf("%s", filename);

	unsigned int text_length = 0;
	f1 = fopen(filename, "r");


	if (NULL == f1)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filename);
		return;
	}


	unsigned int number_of_total_characters = 0;

	while (!feof(f1))
	{
		fscanf(f1, "%c", &ch);
		if ((ch == 'A') || (ch == 'C') || (ch == 'G') || (ch == 'T') ||
			(ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't'))
			text_length++;
		else
		{

			fprintf(stdout, "%u-th character does not belong to {a, c, g, t, A, C, G, T}. Its  ASCII Code is %u, and it is %c.\n",
				number_of_total_characters, ch, ch);
		}

		number_of_total_characters++;

	}

	printf("text_length=%llu\n", text_length);

	fclose(f1);


	char* refer = (char *)malloc(sizeof(char)*(text_length + 1));




	f1 = fopen(filename, "r");
	for (i = 0; i <text_length; i++)
	{
		fscanf(f1, "%c", &ch);

		refer[i] = ch;
	}



	last_bot = 0;


	for (i = 0; i < pre_hashTableMaxSize; i++)
	{

		convert_2bits_to_chars(i, seeds, WINDOW_SIZE);



		counting = count(seeds, WINDOW_SIZE, &top, &bot, &pre_top, &pre_bot);



		fprintf(stderr, "**********************************\n");
		fprintf(stderr, "counting=%u\n", counting);
		fprintf(stderr, "seeds=%s\n", seeds);
		fprintf(stderr, "top=%u\n", top);
		fprintf(stderr, "bot=%u\n", bot);





		if (last_bot != top&&counting!=0)
		{





			tmp_SA_length = 0;



			locates = (unsigned int*)malloc(sizeof(unsigned int)*(top - last_bot));

			bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
				locates,
				last_bot,
				top, 0,
				&tmp_SA_length, bitmapper_index_params.compress_sa);

			fprintf(stderr, "###################################\n");
			fprintf(stderr, "tmp_SA_length=%u\n", tmp_SA_length);
			fprintf(stderr, "counting=%u\n", counting);
			fprintf(stderr, "seeds=%s\n", seeds);
			fprintf(stderr, "top=%u\n", top);
			fprintf(stderr, "bot=%u\n", bot);
			fprintf(stderr, "last_bot=%u\n", last_bot);
			unsigned int j;
			for (j = 0; j < tmp_SA_length; j++)
			{
				fprintf(stderr, "SA[%u]=%u\n", last_bot + j, locates[j]);
				fprintf(stderr, "%c", refer[locates[j] + 0]);
				fprintf(stderr, "%c", refer[locates[j] + 1]);
				fprintf(stderr, "%c", refer[locates[j] + 2]);
				fprintf(stderr, "%c", refer[locates[j] + 3]);
				fprintf(stderr, "%c", refer[locates[j] + 4]);
				fprintf(stderr, "%c", refer[locates[j] + 5]);
				fprintf(stderr, "%c", refer[locates[j] + 6]);
				fprintf(stderr, "%c", refer[locates[j] + 7]);
				fprintf(stderr, "%c", refer[locates[j] + 8]);
				fprintf(stderr, "%c", refer[locates[j] + 9]);
				fprintf(stderr, "%c", refer[locates[j] + 10]);
				fprintf(stderr, "%c", refer[locates[j] + 11]);
				fprintf(stderr, "%c", refer[locates[j] + 12]);
				fprintf(stderr, "\n");
			}

			free(locates);









			fprintf(stderr, "ERRORR!\n");
















		}

		if (counting!=0)
		{
			last_bot = bot;
		}




	}



	return;









	for (i = 1; i <= num_reads; i++)
	{

		if (i % 10000 == 0)
		{
			fprintf(stdout, "i=%llu\n", i);
		}



		number_of_hits
			= count(reads, length_read, &top, &bot, &pre_top, &pre_bot);

		if (number_of_hits == 0)
		{
			reads = reads + length_read;
			continue;
		}


		locates = (unsigned int *)malloc(number_of_hits*sizeof(unsigned int));

		tmp_SA_length = 0;


		locate(reads, top, bot,
			pre_top, pre_bot, locates, length_read, &tmp_SA_length);






		free(locates);



		number_of_locations = number_of_locations + tmp_SA_length;




		reads = reads + length_read;

	}

	gettimeofday(&end_timeval, NULL);
	timer = 1000000 * (end_timeval.tv_sec - start_timeval.tv_sec) + end_timeval.tv_usec - start_timeval.tv_usec;


	double finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;


	fprintf(stdout, "\n\n\n*********************Result*********************\n");


	fprintf(stdout, "searching Time: %ld microsecond\n", timer);

	fprintf(stdout, "searching Time: %5.3f seconds\n", duration);

	fprintf(stdout, "number of matched locations=%ld\n", number_of_locations);


	fprintf(stdout, "************************************************\n");

}

**/



void search_from_bwt_less_than_4()
{
	long long i, j;
	int length_read;
	char* reads;
	unsigned int top, bot, t;
	unsigned int pre_top, pre_bot;


	FILE* _ih_fp = fopen("patterns.txt", "r");

	unsigned int numocc;
	unsigned int query_length = 0;

	fread(&(query_length), sizeof(unsigned int), 1, _ih_fp);
	fread(&(numocc), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "query_length=%u\n", query_length);
	fprintf(stdout, "numocc=%u\n", numocc);


	char* patterns = (char *)malloc((query_length*numocc)*(sizeof(char)));

	fread(patterns, sizeof(char), (query_length*numocc), _ih_fp);


	reads = patterns;

	reads = patterns;

	length_read = query_length;

	unsigned int num_reads = numocc;

	unsigned int* locates;

	double start = clock();


	long long number_of_locations = 0;



	unsigned int tmp_SA_length = 0;





	long long number_of_hits = 0;




	struct  timeval  start_timeval;
	struct  timeval  end_timeval;
	unsigned long timer;
	gettimeofday(&start_timeval, NULL);




	for (i = 1; i <= num_reads; i++)
	{

		if (i % 10000 == 0)
		{
			fprintf(stdout, "i=%llu\n", i);
		}


		/**
		number_of_hits
			= count(reads, length_read, &top, &bot, &pre_top, &pre_bot);
		**/



		number_of_hits
			= count_less_than_4(reads, length_read, &top, &bot);


		if (number_of_hits == 0)
		{
			reads = reads + length_read;
			continue;
		}




		locates = (unsigned int *)malloc(number_of_hits*sizeof(unsigned int));

		tmp_SA_length = 0;



        /**
		locate_less_than_4(reads, top, bot,
			pre_top, pre_bot, locates, length_read, &tmp_SA_length);
        **/


		/**
		qsort(locates, bot - top, sizeof(unsigned int), unsigned_int_compareEntrySize);

		for (t = 0; t<bot - top; t++)
		{

		fprintf(stderr, "i=%llu, site=%u\n", i, locates[t]);
		fflush(stderr);
		}
		**/




		free(locates);






		number_of_locations = number_of_locations + tmp_SA_length;




		reads = reads + length_read;

	}

	gettimeofday(&end_timeval, NULL);
	timer = 1000000 * (end_timeval.tv_sec - start_timeval.tv_sec) + end_timeval.tv_usec - start_timeval.tv_usec;


	double finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;


	fprintf(stdout, "\n\n\n*********************Result*********************\n");


	fprintf(stdout, "searching Time: %ld microsecond\n", timer);

	fprintf(stdout, "searching Time: %5.3f seconds\n", duration);

	fprintf(stdout, "number of matched locations=%ld\n", number_of_locations);


	fprintf(stdout, "************************************************\n");

}







int search_bwt()
{

	char filename[100];
	fprintf(stdout, "Please input the prefix of index name:\n");
	scanf("%s", filename);

	load_index(filename);

	///return 0;

	///if (bitmapper_index_params.compress_sa >= 4)
	if (bitmapper_index_params.compress_sa >= 4)
	{

		///test_hash();

		search_from_bwt_more_than_3();
	}
	else
	{
		///fprintf(stdout, "compress_sa < 4! \n");
		search_from_bwt_less_than_4();
	}



	return 0;
}









int main()
{
	int mode;
	fprintf(stdout, "input mode:\n1 for index\n2 for search\n3 make patterns:");
	scanf("%d", &mode);


	if (mode == 1) build_index();
	else if (mode == 2) search_bwt();
	else if (mode == 3) Generate_patterns();
	else
	{
		fprintf(stdout, "Do not input required options! FMtree will exit ...\n");
	}

	return 0;
}
