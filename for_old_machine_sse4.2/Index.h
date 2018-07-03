/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/


#ifndef __HASH_TABLE__
#define __HASH_TABLE__
#include "Schema.h"
typedef struct HashTable
{
	long long hv;
	int *locs;
} HashTable;

typedef struct
{
	unsigned int *locs;
	unsigned int* hashTable_order;
} Commpress_HashTable;

unsigned int	hashVal(char *seq);
unsigned char			*getRefGenome();
int				Start_Load_Indexz(char *fileName);
Commpress_HashTable	getHashTable();

//下面是一堆函数指针啊
void 			createIndex(char *fileName, char *indexName);
int				Load_Index(int errThreshould, _rg_name_l  **_ih_refGenName, bitmapper_bs_iter* msf_refChromeCont, char* indexName);
///void			(*finalizeLoadingHashTable)();
///unsigned int	*(*getCandidates)(int hv);
int Start_Load_Index(char *fileName);
bitmapper_bs_iter getRefGenomeLength();
bitmapper_bs_iter getRefGenomeLength_2bit();
int  hashVal_EN(char *seq, unsigned int *N_site, unsigned int *hash_key);
void add_hashVal(char *seq, unsigned int *N_site, unsigned int *hash_key);

#endif
