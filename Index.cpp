/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Auxiliary.h"
#include "Ref_Genome.h"
#include "Index.h"
#include "bwt.h"

/**********************************************/

FILE		*_ih_fp			= NULL;
FILE		*_ih_fp_bs_ref = NULL;
extern unsigned int	*_ih_hashTable		= NULL;
extern bitmapper_bs_iter total_SA_length;
unsigned int*  _ih_hashTable_key=NULL;
int 		_ih_maxHashTableSize	= 0;
char		*_ih_refGen		= NULL;

///基因组总长
bitmapper_bs_iter refGenLength;
///基因组2-bit形式总厂
bitmapper_bs_iter refGenLength_2_bit;

//_rg_name_l  *_ih_refGenName=NULL;



unsigned int hashVal(char *seq)
{
  unsigned int i=0;
  unsigned int val=0, numericVal=0;
  while(i<WINDOW_SIZE)
    {
      switch (seq[i])
	{
	case 'A':
	  numericVal = 0;
	  break;
	case 'C':
	  numericVal = 1;
	  break;
	case 'G' :
	  numericVal = 2;
	  break;
	case 'T':
	  numericVal = 3;
	  break;
	default:
	  return (unsigned int)-1;
	  break;
	}
      val = (val << 2)|numericVal;
      i++;
    }
  return val;
}

int  hashVal_EN(char *seq,unsigned int *N_site,unsigned int *hash_key)
{
  unsigned int i=0;
  unsigned int val=0, numericVal=0;
  while(i<WINDOW_SIZE)
    {
      switch (seq[i])
	{
	case 'A':
	  numericVal = 0;
	  break;
	case 'C':
	  numericVal = 1;
	  break;
	case 'G' :
	  numericVal = 2;
	  break;
	case 'T':
	  numericVal = 3;
	  break;
	default:
      *N_site=i;
	  (*hash_key)=(unsigned int)-1;
	  return 1;
	  break;
	}
      val = (val << 2)|numericVal;
      i++;
    }
    *hash_key=val;
}

void add_hashVal(char *seq,unsigned int *N_site,unsigned int *hash_key)
{
    unsigned int mask=(unsigned int)-1;
    mask=mask>>(32-WINDOW_SIZE*2);
    switch (seq[WINDOW_SIZE-1])
	{
	case 'A':
	  (*hash_key)=(*hash_key)<<2;
	  (*hash_key)=(*hash_key)&mask;
	  break;
	case 'C':
      (*hash_key)=(*hash_key)<<2;
      (*hash_key)++;
	  (*hash_key)=(*hash_key)&mask;
	  break;
	case 'G' :
      (*hash_key)=(*hash_key)<<2;
      (*hash_key)=(*hash_key)+2;
	  (*hash_key)=(*hash_key)&mask;
	  break;
	case 'T':
	    (*hash_key)=(*hash_key)<<2;
	    (*hash_key)=(*hash_key)+3;
	  (*hash_key)=(*hash_key)&mask;
	  break;
	default:
      *N_site=WINDOW_SIZE-1;
	  (*hash_key)=(unsigned int )-1;
	  break;
	}
}

void initSavingIHashTable(char *fileName, bitmapper_bs_iter reference_length, _rg_name_l **refChromeName1,
	bitmapper_bs_iter refChromeCont)
{
	  bitmapper_bs_iter tmp;

	  _ih_fp = fopen (fileName, "w");
	  tmp = fwrite(&refChromeCont, sizeof(refChromeCont), 1, _ih_fp);

	  bitmapper_bs_iter len = 0;
	  bitmapper_bs_iter i = 0;


	  for (i = 0; i < refChromeCont; i++)
	  {
		  len = strlen((*refChromeName1)[i]._rg_chrome_name);
		  tmp = fwrite(&len, sizeof(len), 1, _ih_fp);
		  tmp = fwrite((*refChromeName1)[i]._rg_chrome_name, sizeof(char), len, _ih_fp);
		  tmp = fwrite(&((*refChromeName1)[i]._rg_chrome_length), 
		  sizeof(((*refChromeName1)[i]._rg_chrome_length)), 1, _ih_fp);
	   }


	   //写入参考基因组与参考基因组长度
	   tmp = fwrite(&reference_length, sizeof(reference_length), 1, _ih_fp);

}



FILE* get_index_file()
{
	return _ih_fp;
}

int Load_Methy_Index(int errThreshould, _rg_name_l  **msf_ih_refGenName, bitmapper_bs_iter* msf_refChromeCont, char* indexName)
{


	_ih_fp = fopen(indexName, "r");

	if (_ih_fp == NULL)
	{
		fprintf(stdout, "Cannot open %s!\n", indexName);
		return 0;
	}
	else
	{
		fprintf(stdout, "Open %s sucessfully!\n", indexName);
	}


	bitmapper_bs_iter len;

	bitmapper_bs_iter tmpSize;
	bitmapper_bs_iter refChromeCont;

	bitmapper_bs_iter tmp;//tmp各种打酱油，各种做临时暂存变量
	bitmapper_bs_iter i = 0, j = 0;



	//读入染色体数目
	tmp = fread(&refChromeCont, sizeof(refChromeCont), 1, _ih_fp);

	///接下来是读入各染色体信息
	_rg_name_l* _ih_refGenName = (_rg_name_l*)malloc(sizeof(_rg_name_l)*refChromeCont);

	i = 0;
	len = 0;
	bitmapper_bs_iter tmp_start = 0;
	bitmapper_bs_iter tmp_end = 0;
	for (i = 0; i < refChromeCont; i++)
	{
		tmp = fread(&len, sizeof(len), 1, _ih_fp);

		char* tmp_refGenName = (char*)malloc(sizeof(char)*(len + 1));
		tmp = fread(tmp_refGenName, sizeof(char), len, _ih_fp);
		tmp_refGenName[len] = '\0';
		strcpy(_ih_refGenName[i]._rg_chrome_name, tmp_refGenName);

		tmp = fread(&(_ih_refGenName[i]._rg_chrome_length), sizeof (_ih_refGenName[i]._rg_chrome_length), 1, _ih_fp);
		_ih_refGenName[i].start_location = tmp_start;

		tmp_end = tmp_start + _ih_refGenName[i]._rg_chrome_length - 1;
		_ih_refGenName[i].end_location = tmp_end;
		tmp_start = tmp_end + 1;
		free(tmp_refGenName);
	}


	///染色体的各种信息
	*msf_ih_refGenName = _ih_refGenName;
	///有多少条染色体
	*msf_refChromeCont = refChromeCont;



	///这个是基因组总长
	tmp = fread(&refGenLength, sizeof(refGenLength), 1, _ih_fp);

	fprintf(stdout, "refGenLength = %llu \n", refGenLength);


	return 1;
}


void SavingMethyIndex(char *fileName, bitmapper_bs_iter reference_length, _rg_name_l **refChromeName1,
	bitmapper_bs_iter refChromeCont, char* refGen)
{
	bitmapper_bs_iter tmp;

	FILE* methy_ih_fp;

	char methy_fileName[SEQ_MAX_LENGTH];

	sprintf(methy_fileName, "%s.methy", fileName);

	methy_ih_fp = fopen(methy_fileName, "w");
	tmp = fwrite(&refChromeCont, sizeof(refChromeCont), 1, methy_ih_fp);

	bitmapper_bs_iter len = 0;
	bitmapper_bs_iter i = 0;


	for (i = 0; i < refChromeCont; i++)
	{
		len = strlen((*refChromeName1)[i]._rg_chrome_name);
		tmp = fwrite(&len, sizeof(len), 1, methy_ih_fp);
		tmp = fwrite((*refChromeName1)[i]._rg_chrome_name, sizeof(char), len, methy_ih_fp);
		tmp = fwrite(&((*refChromeName1)[i]._rg_chrome_length),
			sizeof(((*refChromeName1)[i]._rg_chrome_length)), 1, methy_ih_fp);
	}


	//写入参考基因组与参考基因组长度
	tmp = fwrite(&reference_length, sizeof(reference_length), 1, methy_ih_fp);

	char tmp_char;

	for (i = 0; i < reference_length; i++)
	{
		refGen[i] = toupper(refGen[i]);

		switch (refGen[i])
		{
		case 'A':
			tmp_char = 0;
			break;
		case 'C':
			tmp_char = 1;
			break;
		case 'G':
			tmp_char = 2;
			break;
		case 'T':
			tmp_char = 3;
			break;
		default:
			tmp_char = 4;
			break;
		}
		tmp = fwrite(&tmp_char, sizeof(tmp_char), 1, methy_ih_fp);
	}

	fseek(methy_ih_fp, 0, SEEK_END);
	bitmapper_bs_iter file_length = ftell(methy_ih_fp);
	fseek(methy_ih_fp, 0, SEEK_SET);
	/**
	if (file_length != reference_length)
	{
		fprintf(stdout, "ERRORE: file_length!=reference_length\n");
		exit(0);
	}
	**/

	fclose(methy_ih_fp);

}



void saveIHashTable(unsigned int  whole_length, unsigned int *hashTable, unsigned int * hashTable_order,  unsigned int maxSize , _rg_name_l **refChromeName1,
int refChromeCont, unsigned int reference_length, char** refGen)
{

  int tmp;
  //这里写入的是索引大小
  //写字符0代表没有额外信息
  unsigned char extraInfo = 0;
  tmp = fwrite (&extraInfo, sizeof(extraInfo), 1, _ih_fp);
  tmp = fwrite(&refChromeCont, sizeof(refChromeCont), 1, _ih_fp);
  short len = 0;
  unsigned int i = 0;
  for (i = 0; i < refChromeCont; i++)
  {
	  len = strlen((*refChromeName1)[i]._rg_chrome_name);
	  tmp = fwrite(&len, sizeof(len), 1, _ih_fp);
	  tmp = fwrite((*refChromeName1)[i]._rg_chrome_name,sizeof(char), len, _ih_fp);
	  tmp = fwrite(&((*refChromeName1)[i]._rg_chrome_length), sizeof(unsigned int), 1, _ih_fp);
  }


    //写入参考基因组与参考基因组长度
  tmp = fwrite(&reference_length, sizeof(reference_length), 1, _ih_fp);
  //这里把参考基因组本身和索引写一块了，搜索的时候就不再需要
  tmp = fwrite(*refGen, sizeof(char), reference_length, _ih_fp);


    //写入Hash表内容项的长度
  tmp = fwrite(&whole_length, sizeof(whole_length), 1, _ih_fp);
  //写入Hash表内容项
  tmp = fwrite(hashTable, sizeof(unsigned int), whole_length, _ih_fp);


//写入Hash表头的长度(k-mer)
  tmp = fwrite(&maxSize, sizeof(maxSize), 1, _ih_fp);
//写入Hash表头的内容(k-mer)
  tmp = fwrite(hashTable_order, sizeof(unsigned int), maxSize, _ih_fp);

  if (tmp == 0)
    fprintf(stderr, "Write error while saving hash table.\n");
}


void build_inde_from_disk(char* file_name)
{
	FILE *f1;

	f1 = fopen(file_name, "r");



	if (NULL == f1)
	{
		fprintf(stdout, "Failed to open %s! Bitmapper_BS will exit ...\n", file_name);
		exit(1);
	}
	else
	{
		fprintf(stdout, "Open %s sucessfully! \n", file_name);
	}


	bitmapper_bs_iter number_of_total_characters = 0;

	char ch;
	bitmapper_bs_iter text_length = 0;

	while (1)
	{
		fscanf(f1, "%c", &ch);

		if (feof(f1))
		{
			break;
		}

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



	bitmapper_bs_iter i = 0;
	f1 = fopen(file_name, "r");
	for (i = 0; i <text_length; i++)
	{
		fscanf(f1, "%c", &ch);

		refer[i] = ch;
	}

	bitmapper_bs_iter compress_build_sa = 8;

	indenpendent_creadte_index(text_length, &refer, compress_build_sa, file_name);



}


///到这个函数的时候, refGen里面已经没有N了才对
void generate_directional_BS_genome_to_disk(char* refGen, bitmapper_bs_iter length, char* file_name)
{
	char convert_complement[256];

	long long i = 0;

	FILE *f1;

	f1 = fopen(file_name, "w");


	for (i = 0; i < length; ++i)
	{
		if (refGen[i] != 'A'&&
			refGen[i] != 'C'&&
			refGen[i] != 'G'&&
			refGen[i] != 'T'&&
			refGen[i] != 'a'&&
			refGen[i] != 'c'&&
			refGen[i] != 'g'&&
			refGen[i] != 't')
		{
			fprintf(stdout, "(generate_directional_BS_genome) Text has non-A-C-G-T characters!\n");
			fprintf(stdout, "(generate_directional_BS_genome) i = %c\n", i);
			fprintf(stdout, "(generate_directional_BS_genome) character = %c\n", refGen[i]);
			fprintf(stdout, "(generate_directional_BS_genome) ASCII code = %u\n", (unsigned int)refGen[i]);
			exit(1);
		}
	}




	for (i = 0; i < 256; i++)
	{
		convert_complement[i] = 'N';
	}

	convert_complement['A'] = 'T';
	convert_complement['a'] = 'T';

	convert_complement['C'] = 'G';
	convert_complement['c'] = 'G';

	convert_complement['G'] = 'C';
	convert_complement['g'] = 'C';

	convert_complement['T'] = 'A';
	convert_complement['t'] = 'A';


	

	/********************************************Output complement strand*****************************************/
	for (i = 0; i < length; ++i)
	{

		char buffer = convert_complement[refGen[i]];

		if (buffer == 'C'
			||
			buffer == 'c')
		{
			buffer = 'T';
		}

		
		fprintf(f1, "%c", buffer);
	}

	/********************************************Output complement strand******************************************/



	/********************************************Output reverse forward strand****************************************************/
	for (i = length - 1; i >= 0; --i)
	{



		char buffer = refGen[i];

		if (buffer == 'C'
			||
			buffer == 'c')
		{
			buffer = 'T';
		}

		fprintf(f1, "%c", buffer);

	}

	/********************************************Output reverse forward strand****************************************************/



	fflush(f1);

	fclose(f1);

}


///把refGen里的N随机转换
void replace_N(char* refGen, bitmapper_bs_iter length)
{
	long long i = 0;

	char buffer;


	srand((unsigned int)time(0));
	char randome_char[4];
	randome_char[0] = 'A';
	randome_char[1] = 'C';
	randome_char[2] = 'G';
	randome_char[3] = 'T';


	for (i = 0; i < length; ++i)
	{
		if (refGen[i] != 'A'&&
			refGen[i] != 'C'&&
			refGen[i] != 'G'&&
			refGen[i] != 'T'&&
			refGen[i] != 'a'&&
			refGen[i] != 'c'&&
			refGen[i] != 'g'&&
			refGen[i] != 't')
		{
			buffer = randome_char[rand() % 4];

			refGen[i] = buffer;
		}
	}


}



///把refGen转成2-bit编码
void convert_to_2bit(char* refGen, bitmapper_bs_iter length, char* file_name)
{
	bitmapper_bs_iter i = 0;
	bitmapper_bs_iter i_2_bit = 0;

	unsigned char convert_2bit[256];

	convert_2bit['A'] = 0;
	convert_2bit['a'] = 0;

	convert_2bit['C'] = 1;
	convert_2bit['c'] = 1;

	convert_2bit['G'] = 2;
	convert_2bit['g'] = 2;

	convert_2bit['T'] = 3;
	convert_2bit['t'] = 3;


	unsigned char tmp_char;


	while (i + 4 <= length)
	{
		tmp_char = 0;


		tmp_char = tmp_char | convert_2bit[refGen[i]];
		tmp_char = tmp_char << 2;


		tmp_char = tmp_char | convert_2bit[refGen[i + 1]];
		tmp_char = tmp_char << 2;


		tmp_char = tmp_char | convert_2bit[refGen[i + 2]];
		tmp_char = tmp_char << 2;


		tmp_char = tmp_char | convert_2bit[refGen[i + 3]];
		///tmp_char = tmp_char << 2;



		refGen[i_2_bit] = tmp_char;

		i = i + 4;

		i_2_bit++;
	}

	///这里不等于，则代表还有尾巴没处理完
	if (i != length)
	{
		long long inner_i = 0;
		long long last_length = length - i;

		tmp_char = 0;

		///last_length最大就是3吧
		for (inner_i = 0; inner_i < last_length; inner_i++)
		{
			tmp_char = tmp_char << 2;

			tmp_char = tmp_char | convert_2bit[refGen[i + inner_i]];

		}

		tmp_char = tmp_char << ((4 - last_length) * 2);

		refGen[i_2_bit] = tmp_char;

		i_2_bit++;

	}


	/**
	fprintf(stdout, "i_2_bit = %llu\n", i_2_bit);

	for (i = 0; i < i_2_bit; i++)
	{
		fprintf(stdout, "refGen[%llu] = %llu\n", i, refGen[i]);
	}
	**/

	char file_name_g[1000];

	strcpy(file_name_g, file_name);
	strcpy(file_name_g + strlen(file_name), ".pac");

	FILE* f_g = fopen(file_name_g, "w");

	if (f_g == NULL)
	{
		fprintf(stdout, "Cannot open %s \n", file_name_g);
		exit(1);
	}

	fwrite(&i_2_bit, sizeof(i_2_bit), 1, f_g);
	fwrite(refGen, sizeof(char), i_2_bit, f_g);
	fflush(f_g);
	fclose(f_g);
}


void createIndex(char *fileName, char *indexName)
{
  //开始计时
  double startTime = Get_T();
  //这个一开始充当的是参考基因组本身
  char *refGen;
  //染色体名称数组
  _rg_name_l *refChromeName = NULL;
  //染色体条数
  bitmapper_bs_iter refChromeCont = 0;

  bitmapper_bs_iter i;


  //获得参考基因组文件长度
  ///这里的结果比较粗略，实际上是文件大小,这个会比基因组略大一点
  bitmapper_bs_iter reference_length = initLoadingRefGenome(fileName);

  fprintf(stdout, "Generating Index from %s \n", fileName);
	

  if (!reference_length)
  {
    fprintf(stdout, "The error of Reading Index:\n");
    return;
  }


  


  ///这个函数结束之后，reference_length是真正的基因组长度
  loadRefGenome(&refGen, &refChromeName, &refChromeCont, &reference_length);



  for(i=0;i<refChromeCont;++i)
  {
      fprintf(stdout, "Chrome Name =%s ********\n", (refChromeName)[i]._rg_chrome_name);
      fflush(stdout);
      fprintf(stdout, "Chrome length =%d ********\n", (refChromeName)[i]._rg_chrome_length);
      fflush(stdout);
  }

  fprintf(stdout, "The length of Genome is %u:\n", reference_length);
  fflush(stdout);

  //这就是把那个文件打开而已..，写进了一些属性信息
  ///主要是写基因组长度
  ///写染色体信息
  initSavingIHashTable(indexName, reference_length, &refChromeName, refChromeCont);

  SavingMethyIndex(indexName, reference_length, &refChromeName, refChromeCont, refGen);

  ///char file_name[30] = "tmp_bs_ref";

  char file_name[1000];

  strcpy(file_name, indexName);
  strcpy(file_name + strlen(indexName), ".bs");

  char command_string[1000] = "rm ";

  strcpy(command_string + strlen(command_string), file_name);


  replace_N(refGen, reference_length);


  generate_directional_BS_genome_to_disk(refGen, reference_length, file_name);


  ///这中间还得写一段，就是把基因组转成二进制形式写文件吧
  ///基因组写硬盘
  convert_to_2bit(refGen, reference_length, file_name);

  free(refGen);

  ///索引写硬盘
  build_inde_from_disk(file_name);


  

  fprintf(stdout, "command_string: %s\n", command_string);

  int error = system(command_string);
  error = system("rm *.sa5");

  



//释放内存，然后把当前可用内存使用量更新
  //free(refGen);
  //free(refChromeName);
  fclose(_ih_fp);






  fprintf(stdout, "\nDONE in %0.2fs!\n", (Get_T()-startTime));
}

int  Load_Index(int errThreshould, _rg_name_l  **msf_ih_refGenName, bitmapper_bs_iter* msf_refChromeCont, char* indexName)
{

	bitmapper_bs_iter len;

	bitmapper_bs_iter tmpSize;
	bitmapper_bs_iter refChromeCont;

	bitmapper_bs_iter tmp;//tmp各种打酱油，各种做临时暂存变量
	bitmapper_bs_iter i = 0, j = 0;



   //读入染色体数目
   tmp = fread(&refChromeCont, sizeof(refChromeCont), 1, _ih_fp);

   ///接下来是读入各染色体信息
   _rg_name_l* _ih_refGenName = (_rg_name_l*)malloc(sizeof(_rg_name_l)*refChromeCont);

   i = 0;
   len = 0;
   bitmapper_bs_iter tmp_start = 0;
   bitmapper_bs_iter tmp_end = 0;
	for (i = 0; i < refChromeCont; i++)
	{
		tmp = fread(&len, sizeof(len), 1, _ih_fp);

		char* tmp_refGenName = (char*)malloc(sizeof(char)*(len + 1));
		tmp = fread(tmp_refGenName, sizeof(char), len, _ih_fp);
		tmp_refGenName[len] ='\0';
		strcpy(_ih_refGenName[i]._rg_chrome_name,tmp_refGenName);

		tmp = fread(&(_ih_refGenName[i]._rg_chrome_length), sizeof (_ih_refGenName[i]._rg_chrome_length), 1, _ih_fp);
		_ih_refGenName[i].start_location=tmp_start;

		tmp_end=tmp_start+_ih_refGenName[i]._rg_chrome_length-1;
		_ih_refGenName[i].end_location=tmp_end;
		tmp_start=tmp_end+1;
		free(tmp_refGenName);
	}


  ///染色体的各种信息
  *msf_ih_refGenName=_ih_refGenName;
  ///有多少条染色体
  *msf_refChromeCont = refChromeCont;



  ///这个是基因组总长
  tmp = fread(&refGenLength, sizeof(refGenLength), 1, _ih_fp);

  fprintf(stdout, "refGenLength = %llu \n", refGenLength);


  ///这一块是索引的名字
  char file_name[1000];
  strcpy(file_name, indexName);
  strcpy(file_name + strlen(indexName), ".bs");


  ///这一块是参考组的名字
  char file_name_g[1000];
  strcpy(file_name_g, file_name);
  strcpy(file_name_g + strlen(file_name), ".pac");


  ///打开参考组，读入参考组
  _ih_fp_bs_ref = fopen(file_name_g, "r");
  if (_ih_fp_bs_ref == NULL)
  {
	  fprintf(stdout, "Cannot open %s \n", file_name_g);
	  exit(1);
  }
  else
  {
	  fprintf(stdout, "Open %s sucessfully!\n", file_name_g);
  }
  ///这个是基因组2bit形式总长
  fread(&refGenLength_2_bit, sizeof(refGenLength_2_bit), 1, _ih_fp_bs_ref);
  //加1000是为了保持全敏感性
  _ih_refGen = (char*)malloc(sizeof(char)*(refGenLength_2_bit + 1 + 1000));
  fread(_ih_refGen, sizeof(char), refGenLength_2_bit, _ih_fp_bs_ref);


  ///读入FM-index
  load_index(file_name);


  bitmapper_bs_iter site_0, length_0;

  locate_one_position(&site_0, 0, &length_0);

  total_SA_length = site_0;

  ///fprintf(stderr, "total_SA_length: %llu, refGenLength: %llu \n", total_SA_length, refGenLength);

  ///fprintf(stderr, "site_0: %llu, length_0: %llu \n", site_0, length_0);

  return 1;
}



int Start_Load_Index(char *fileName)
{


  _ih_fp = fopen (fileName, "r");

  if (_ih_fp == NULL)
  {
	  fprintf(stdout, "Cannot open %s!\n", fileName);
	  return 0;
  }
  else
  {
	  fprintf(stdout, "Open %s sucessfully!\n", fileName);
  }

  return 1;
}



unsigned char *getRefGenome()
{
	return (unsigned char*)_ih_refGen;

}

bitmapper_bs_iter getRefGenomeLength()
{
  return refGenLength;

}

bitmapper_bs_iter getRefGenomeLength_2bit()
{
	return refGenLength_2_bit;

}


Commpress_HashTable	getHashTable()
{

    Commpress_HashTable chhy_HashTable;
    chhy_HashTable.locs=_ih_hashTable;
    chhy_HashTable.hashTable_order=_ih_hashTable_key;

  return chhy_HashTable;
}
