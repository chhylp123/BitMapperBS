/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Auxiliary.h"
#include "Ref_Genome.h"
#include <unistd.h>

FILE *_rg_fp;


bitmapper_bs_iter initLoadingRefGenome(char *fileName)
{
  _rg_fp = fopen (fileName, "r");

  if (_rg_fp == NULL)
  {
	  return 0;
  }
  bitmapper_bs_iter file_length = 0;
  fseek(_rg_fp, 0, SEEK_END);
  file_length = ftell(_rg_fp);
  fseek(_rg_fp, 0, SEEK_SET);
  return file_length ;
}


void loadRefGenome(char **refGen, _rg_name_l **refGenName, bitmapper_bs_iter* refChromeCont, bitmapper_bs_iter* refGenOff)
{
  char ch;
  //计数参考基因组大小
  bitmapper_bs_iter _rg_contGen = 0;
  //计数当前染色体大小
  bitmapper_bs_iter _rg_contChromeLength = 0;
  //计数染色体条数
  bitmapper_bs_iter _rg_contChrome = 0;



  (*refGen) = (char*)malloc(sizeof(char)*(*refGenOff));
  (*refGenName) = (_rg_name_l*)malloc(sizeof(_rg_name_l) * 1);


  while( fscanf(_rg_fp, "%c", &ch) > 0 )
    {

      if (ch == '>')
	  {
        fprintf(stderr, "%c\n", ch);
        _rg_name_l*  tmp_rg_name = (struct _rg_name_l*)malloc(sizeof(struct _rg_name_l)* (++_rg_contChrome));

	  int i = 0;

	  for (i = 0; i < (_rg_contChrome-1); i++)
	  {
		  strcpy(tmp_rg_name[i]._rg_chrome_name, (*refGenName)[i]._rg_chrome_name);
		  tmp_rg_name[i]._rg_chrome_length = (*refGenName)[i]._rg_chrome_length;
	  }

	  char *tmp;
	  tmp = fgets(tmp_rg_name[_rg_contChrome-1]._rg_chrome_name, SEQ_MAX_LENGTH, _rg_fp);
      if((_rg_contChrome-2)>=0)
      {
        tmp_rg_name[_rg_contChrome-2]._rg_chrome_length=_rg_contChromeLength;
      }
	  if (tmp == NULL)
		  fprintf(stderr, "Error reading the reference.\n");
	  free((*refGenName));
	  (*refGenName) = tmp_rg_name;
	  //下面这个就是确定字符串_rg_name的末尾罢了
	  int k;
	  for (k = 0; k<strlen((*refGenName)[_rg_contChrome - 1]._rg_chrome_name); k++)
	  {
		  if ((*refGenName)[_rg_contChrome-1]._rg_chrome_name[k] == ' '||(*refGenName)[_rg_contChrome-1]._rg_chrome_name[k] == '\n')
		  {
			  //到字符串末尾了
			  (*refGenName)[_rg_contChrome-1]._rg_chrome_name[k] = '\0';
			  break;
		  }
	  }
	  fprintf(stderr, "Chrome Name =%s ********\n", (*refGenName)[_rg_contChrome - 1]._rg_chrome_name);
	  _rg_contChromeLength = 0;
	}
      //判断ch是否为空
     else if (!isspace(ch)&&(ch!='\n'))
	{

	  ch = toupper(ch);
	  (*refGen)[_rg_contGen++] = ch;
	  _rg_contChromeLength++;
	}
    }




    //处理最后一个参考基因组
    (*refGenName)[_rg_contChrome-1]._rg_chrome_length=_rg_contChromeLength;

  ///不要这个东西, 要了还麻烦
  ///(*refGen)[_rg_contGen] = '\0';
  *refChromeCont = _rg_contChrome;
  *refGenOff = _rg_contGen;
  fclose(_rg_fp);
}
