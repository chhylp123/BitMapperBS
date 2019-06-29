

///g++ -O3 bam_prase.cpp -lhts -o bam_prase

#include "Auxiliary.h"
#include"bam_prase.h"
#include "Process_sam_out.h"
#include <pthread.h>
#include <htslib/sam.h>


bam_operations bitmapperBS_bam;
pthread_mutex_t queueMutex_bam;





void get_bam_header(bam_hdr_t** bam_header)
{
	samFile* in;
	in = sam_open(bitmapperBS_bam.output_name, "r");

	///读入sam的首部
	if (((*bam_header) = sam_hdr_read(in)) == 0) {
		fprintf(stderr, "fail to read %s\n", bitmapperBS_bam.output_name);
		return;
	}

	///一定要关闭，不关的话写不到文件里去
	sam_close(in);
}

void init_bam_header_back(char* tmp_sam_file_name, char* bam_file_name)
{

	samFile* in;
	in = sam_open(tmp_sam_file_name, "r");
	bitmapperBS_bam.bam_file = sam_open(bam_file_name, "wb");
	bitmapperBS_bam.bam_header = NULL;


	bitmapperBS_bam.output_name = tmp_sam_file_name;

	///读入sam的首部
	if ((bitmapperBS_bam.bam_header = sam_hdr_read(in)) == 0) {
		fprintf(stderr, "fail to read %s\n", tmp_sam_file_name);
		return;
	}

	///写入sam的首部
	if (sam_hdr_write(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header) != 0) {
		fprintf(stderr, "fail to write %s\n", bam_file_name);
		return;
	}

	///一定要关闭，不关的话写不到文件里去
	sam_close(in);

}

void init_bam_header(char* bam_file_name, _rg_name_l  *_ih_refGenName, int refChromeCont, int argc, char *argv[])
{

	if (bam_file_name[0] == '\0')
	{
		bitmapperBS_bam.bam_file = sam_open_format("-", "wb", NULL);
	}
	else
	{
		bitmapperBS_bam.bam_file = sam_open_format(bam_file_name, "wb", NULL);
	}
	

	bitmapperBS_bam.bam_header = NULL;


	bitmapperBS_bam.output_name = bam_file_name;


	char** chrome_name = (char**)malloc(sizeof(char*)*refChromeCont);
	bitmapper_bs_iter* chrome_length = (bitmapper_bs_iter*)malloc(sizeof(bitmapper_bs_iter)*refChromeCont);

	int i = 0;

	for (i = 0; i < refChromeCont; i++)
	{
		chrome_name[i]=_ih_refGenName[i]._rg_chrome_name;
		chrome_length[i] = _ih_refGenName[i]._rg_chrome_length;
	}


	///读入sam的首部
	//if ((bitmapperBS_bam.bam_header = sam_hdr_read(in)) == 0) {
	if ((bitmapperBS_bam.bam_header = chhy_sam_hdr_read(chrome_name, chrome_length, refChromeCont,
	argc, argv, versionN)) == 0) {
		fprintf(stderr, "fail to phrase %s\n", bam_file_name);
		return;
	}



	///写入sam的首部
	if (sam_hdr_write(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header) != 0) {
		fprintf(stderr, "fail to write %s\n", bam_file_name);
		return;
	}

	free(chrome_name);
	free(chrome_length);

}



///void init_bam_file_from_sam(char* file_name, char* outputFileName)
void init_bam_file_from_sam(char* file_name, char* outputFileName,
_rg_name_l  *_ih_refGenName, int refChromeCont, int argc, char *argv[])
{
	finalizeOutput();


	///init_bam_header(file_name, outputFileName);
	///init_bam_header(file_name, outputFileName, _ih_refGenName, refChromeCont, argc, argv);


}




void init_bam_output_cell(bam_output_cell* cell)
{
	cell->b = bam_init1();  ///初始化吧,分配内存
}


void init_bam_output_group(bam_output_group* group, long long total_size)
{
	group->size = total_size;
	group->number = 0;
	group->b = (bam1_t **)malloc(sizeof(bam1_t *)*group->size);

	int i;
	for (i = 0; i < group->size; i++)
	{
		group->b[i] = bam_init1();
	}
}



BGZF* init_buff_bgzf(BGZF* exist_bgzf)
{
	BGZF *fp;
	fp = (BGZF*)calloc(1, sizeof(BGZF));
	if (fp == NULL)
	{
		return NULL;
	}


	fp->is_write = 1;
	fp->is_compressed = 1;

	fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
	if (fp->uncompressed_block == NULL)
	{
		return NULL;
	}


	fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;

	fp->compress_level = exist_bgzf->compress_level;

	fp->is_be = exist_bgzf->is_be;

	fp->fp = exist_bgzf->fp;

	fp->block_offset = exist_bgzf->block_offset;


	///fprintf(stderr, "exist_bgzf->block_offset: %d\n", exist_bgzf->block_offset);

	return fp;

}


void init_bam_buffer(bam_phrase* buf)
{
	buf->bgzf = init_buff_bgzf(bitmapperBS_bam.bam_file->fp.bgzf);
	buf->b = bam_init1();  ///初始化吧,分配内存

	buf->result.uncomp_data = (unsigned char*)malloc(2 * BGZF_MAX_BLOCK_SIZE);
	buf->result.comp_data = (unsigned char*)buf->result.uncomp_data + BGZF_MAX_BLOCK_SIZE;
}


void write_alignment_directly(char* alignment, long long alignment_length, bam_output_cell* cell)
{


	cell->line.l = alignment_length;
	if (cell->line.l == 0)
	{
		return;

	}

	cell->line.m = cell->line.l;
	cell->line.s = alignment;

	free(cell->b->data);
	memset(cell->b, 0, sizeof(bam1_t));
	int ret = sam_parse1(&cell->line, bitmapperBS_bam.bam_header, cell->b);
	sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, cell->b);

}


void phrase_alignment_only(char* alignment, long long alignment_length, bam_output_group* bam_groups)
{


	bam_groups->line.l = alignment_length;
	if (bam_groups->line.l == 0)
	{
		return;

	}

	bam_groups->line.m = bam_groups->line.l;
	bam_groups->line.s = alignment;


	free(bam_groups->b[bam_groups->number]->data);
	memset(bam_groups->b[bam_groups->number], 0, sizeof(bam1_t));
	int ret = sam_parse1(&bam_groups->line, bam_groups->bam_header, bam_groups->b[bam_groups->number]);
	///sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, cell->b);
	bam_groups->number++;

}



void convert_string_to_bam(char* alignment, long long alignment_length, bam_phrase* bam_groups)
{


	bam_groups->line.l = alignment_length;
	if (bam_groups->line.l == 0)
	{
		return;

	}

	bam_groups->line.m = bam_groups->line.l;
	bam_groups->line.s = alignment;


	free(bam_groups->b->data);
	memset(bam_groups->b, 0, sizeof(bam1_t));
	int ret = sam_parse1(&bam_groups->line, bitmapperBS_bam.bam_header, bam_groups->b);
	///sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, cell->b);
	///至少当cigar不超过0xffff时，一行alignment在bam里不会被分开
	chhy_bam_write1_pure(bitmapperBS_bam.bam_file->fp.bgzf, bam_groups->bgzf, bam_groups->b, &(bam_groups->result));
}

void flush_bam_buffer(bam_phrase* bam_groups)
{
	chhy_lazy_flush_pure(bitmapperBS_bam.bam_file->fp.bgzf, bam_groups->bgzf,&(bam_groups->result));
}




void write_alignment_only(bam_output_group* bam_groups)
{
	long long i = 0;

	pthread_mutex_lock(&queueMutex_bam);

	for (i = 0; i < bam_groups->number; i++)
	{
		sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, bam_groups->b[i]);
	}

	bam_groups->number = 0;

	pthread_mutex_unlock(&queueMutex_bam);
}




void close_bam_file()
{
	if (bitmapperBS_bam.bam_header) 
		bam_hdr_destroy(bitmapperBS_bam.bam_header);

	///一定要关闭，不关的话写不到文件里去
	sam_close(bitmapperBS_bam.bam_file);
}


void close_bam_file_rename(char* old_name, char* new_name)
{
	if (remove(new_name) != 0)
	{
		fprintf(stderr, "ERROE\n");
		exit(0);
	}


	if (bitmapperBS_bam.bam_header)
		bam_hdr_destroy(bitmapperBS_bam.bam_header);

	///一定要关闭，不关的话写不到文件里去
	sam_close(bitmapperBS_bam.bam_file);



	if (rename(old_name, new_name)!=0)
	{
		fprintf(stderr, "ERROE\n");
		exit(0);
	}
}



void write_alignment(const char* alignment)
{


	kstring_t line;
	int r;

	line.l = strlen(alignment);
	if (line.l == 0)
	{
		return;

	}

	line.m = line.l;
	line.s = (char*)malloc(line.l + 1);
	memcpy(line.s, alignment, line.l + 1);


	bam1_t *b = bam_init1();  ///初始化吧,分配内存
	int ret = sam_parse1(&line, bitmapperBS_bam.bam_header, b);

	line.l = 0;
	sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, b);
	free(line.s);

	bam_destroy1(b);   ///释放内存
}





///以回车分割
void write_alignment_group(char* alignment, const long long length)
{


	char* tab;

	long long current_length;
	
	kstring_t line;

	bam1_t *b = bam_init1();  ///初始化吧,分配内存
	

	while ((tab = strchr(alignment, '\n')))
	{

		current_length = tab - alignment;

		*tab = '\0';
		
		line.l = current_length;
		if (line.l == 0)
		{
			return;

		}
		
		line.m = line.l;
		line.s = alignment;

		memset(b, 0, sizeof(bam1_t));

		int ret = sam_parse1(&line, bitmapperBS_bam.bam_header, b);

		sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, b);

		alignment = tab + 1;

	}


	bam_destroy1(b);   ///释放内存

}




///以回车分割
void write_alignment_group_pre_allocate(char* alignment, const long long length, const long long number)
{


	char* tab;

	long long current_length;

	kstring_t line;

	///bam1_t *b = bam_init1();  ///初始化吧,分配内存
	bam1_t** b_array = (bam1_t**)malloc(sizeof(bam1_t*)*number);

	long long i;

	for (i = 0; i < number; i++)
	{
		b_array[i] = bam_init1();
	}

	i = 0;

	while ((tab = strchr(alignment, '\n')))
	{
		if (i >= number)
		{
			fprintf(stderr, "ERROR\n");
			exit(0);
		}
		current_length = tab - alignment;

		*tab = '\0';

		line.l = current_length;
		if (line.l == 0)
		{
			return;

		}

		line.m = line.l;
		line.s = alignment;

		memset(b_array[i], 0, sizeof(bam1_t));

		int ret = sam_parse1(&line, bitmapperBS_bam.bam_header, b_array[i]);

		alignment = tab + 1;

		i++;

	}


	for (i = 0; i < number; i++)
	{
		sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, b_array[i]);
	}



	for (i = 0; i < number; i++)
	{
		bam_destroy1(b_array[i]);   ///释放内存
	}

}





void write_alignment_muti_thread(const char* alignment)
{


	kstring_t line;
	int r;

	line.l = strlen(alignment);
	if (line.l == 0)
	{
		return;

	}

	line.m = line.l;
	line.s = (char*)malloc(line.l + 1);
	memcpy(line.s, alignment, line.l + 1);

	bam1_t *b = bam_init1();  ///初始化吧,分配内存
	int ret = sam_parse1(&line, bitmapperBS_bam.bam_header, b);
	line.l = 0;


	pthread_mutex_lock(&queueMutex_bam);
	sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, b);
	pthread_mutex_unlock(&queueMutex_bam);


	
	free(line.s);

	bam_destroy1(b);   ///释放内存
}


/**
void write_alignment(const char* alignment)
{
	
	kstring_t line;
	int r;
	///bitmapperBS_bam.bam_file->line.l = strlen(alignment);
	line.l = strlen(alignment);
	if (bitmapperBS_bam.bam_file->line.l == 0)
	{
		return;
	}
	bitmapperBS_bam.bam_file->line.m = bitmapperBS_bam.bam_file->line.l;
	bitmapperBS_bam.bam_file->line.s = (char*)malloc(bitmapperBS_bam.bam_file->line.l + 1);
	memcpy(bitmapperBS_bam.bam_file->line.s, alignment, bitmapperBS_bam.bam_file->line.l + 1);
	bam1_t *b = bam_init1();  ///初始化吧,分配内存
	int ret = sam_parse1(&bitmapperBS_bam.bam_file->line, bitmapperBS_bam.bam_header, b);
	bitmapperBS_bam.bam_file->line.l = 0;
	sam_write1(bitmapperBS_bam.bam_file, bitmapperBS_bam.bam_header, b);
	free(bitmapperBS_bam.bam_file->line.s);
	bam_destroy1(b);   ///释放内存
}
**/


/**
int main(int argc, char **argv)
{
	///min_example_string_to_bam();
	init_bam_header("output.sam", "output.bam");
	write_alignment("SRR1532535.52049\t16\t4\t61322565\t255\t90M\t*\t0\t0\tAACAATCCCCCCACTTTAACTTCCCATAATACTAAAATAACAAACATAAATCACCACTCCTAATCCCCAATATTTCTTTTATATTAATTA\tCDCC><DDFHHGHEHIHIHHGFGJJIIJGGFIJJJJJIJJJJJIIJIJJJJIGIGEJJJJJJJIJJJJJJJJJJJJIHHHHHFFFFFCCC\tNM:i:0\tZS:Z:-+\0");
	close_bam_file();
}
**/