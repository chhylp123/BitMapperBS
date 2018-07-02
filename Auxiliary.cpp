/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <zlib.h>
#include <string.h>
#include "Auxiliary.h"

//SEQ_LENGTH记录的是read的长度
unsigned short 			SEQ_LENGTH = 0;


/**********************************************/
double Get_T(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec+t.tv_usec/1000000.0;
}

/**********************************************/
char reverseComplementChar(char c)
{
  char ret;
  switch (c)
    {
    case 'A':
      ret = 'T';
      break;
    case 'T':
      ret = 'A';
      break;
    case 'C':
      ret = 'G';
      break;
    case 'G':
      ret = 'C';
      break;
    default:
      ret = 'N';
      break;
    }
  return ret;
}
/**********************************************/
void reverseComplement (char *seq, char *rcSeq , int length)
{
  int i;
  for (i=0; i<length; i++)
  {
	rcSeq[i]=reverseComplementChar (seq[length-1-i]) ;
  }

  rcSeq[length] = 0;
}




/**********************************************/
void reverse (char *seq, char *rcSeq , int length)
{
  int i;
  int l = length;
  if(l != strlen(seq))
    l = strlen(seq);
  for (i=0; i<l; i++)
    {
      rcSeq[i]=seq[l-1-i] ;
    }
  rcSeq[l] = '\0';

}
/**********************************************/
void stripPath(char *full, char **path, char **fileName)
{
  int i;
  int pos = -1;

  for (i=strlen(full)-1; i>=0; i--)
    {
      if (full[i]=='/')
	{
	  pos = i;
	  break;
	}

    }

  if (pos != -1)
    {
      sprintf(*fileName, "%s%c", (full+pos+1), '\0');
      full[pos+1]='\0';
      sprintf(*path,"%s%c", full, '\0');
    }
  else
    {
      sprintf(*fileName, "%s%c", full, '\0');
      sprintf(*path,"%c", '\0');
    }
}
