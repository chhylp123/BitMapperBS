/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#ifndef _REF_GENOME_
#define _REF_GENOME_
#include "Schema.h"
#include "bwt.h"


bitmapper_bs_iter initLoadingRefGenome(char *fileName);
void loadRefGenome(char **refGen, struct _rg_name_l **refGenName, bitmapper_bs_iter* refChromeCont, bitmapper_bs_iter* refGenOff);
#endif
