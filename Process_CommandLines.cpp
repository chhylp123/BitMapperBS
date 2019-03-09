/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#include <stdio.h>
#include<unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include "Auxiliary.h"
#include "Process_CommandLines.h"
#include "bwt.h"


int				is_index = 0;
int				is_search = 0;
int             is_methy = 0;
int				is_pairedEnd;
int				cropSize = 0;
int				minDistance_pair=0;
int				maxDistance_pair=500;
int				over_all_seed_length = 30;
unsigned int			THREAD_COUNT = 1;
char				*Read_File1;
char				*Read_File2;
char				*Mapped_File = "output";
char				*Mapped_FilePath = "";
char                *folder_path = NULL;
char				*Un_Mapped_File = "unmapped";
char				fileName[2][NAME_LENGTH];
unsigned char			thread_e=255;
double thread_e_f = 0.08;
unsigned char			maxHits=0;
unsigned char			WINDOW_SIZE = 11;  //the default value is 11 instead of 12.32位机上这个不能超过15,因为该软件居然丧心病狂的用的int来存hash值...
unsigned char           default_ws=1;
unsigned char			MERGE_STEP= 8;
double bs_score_threshold = -1;
double bs_edit_distance_threshold = -1;
int is_local = 1;
int bs_available_seed_length = -1;
int pbat = 0;
int bam_output = 0;
///int methylation_size = 1000000;
int methylation_size = 100000;
///double methylation_buffer_times = 2.5;
double methylation_buffer_times = 2.5;
int output_methy = 0;
bitmapper_bs_iter genome_cuts = 64;
bitmapper_bs_iter cut_length;
double maxVariantFrac = 0.0;
int minVariantDepth = 65536;
int CpG = 1;
int CHG = 0;
int CHH = 0;
int unmapped_out = 0;
int ambiguous_out = 0;
int mapstats = 0;


void Print_H();

int CommandLine_process (int argc, char *argv[])
{

  int o;
  int index;
  char *fastaFile = NULL;
  folder_path = NULL;

  static struct option longOptions[] =
    {
      {"pe",		no_argument,  	    &is_pairedEnd,		1},
	  { "sensitive", no_argument,      &is_local,  0},
	  { "fast", no_argument, &is_local, 1},
	  { "pbat", no_argument, &pbat, 1},
	  { "bam", no_argument, &bam_output, 1 },
	  { "sam", no_argument, &bam_output, 0 },
	  { "methy_out", no_argument, &output_methy, 1},
	  { "nomethy_out", no_argument, &output_methy, 0},
	  { "CpG", no_argument, &CpG, 1 },
	  { "CHG", no_argument, &CHG, 1 },
	  { "CHH", no_argument, &CHH, 1 },
	  { "mapstats", no_argument, &mapstats, 1},
	  { "unmapped_out", no_argument, &unmapped_out, 1 },
	  { "ambiguous_out", no_argument, &ambiguous_out, 1 },
	  { "methy_extract", required_argument, 0, 'f' },
      {"index",		required_argument,  0, 			'i'},
      {"search",	required_argument,  0,			'g'},
      {"help",		no_argument,	    0,			'h'},
      {"version",	no_argument,	    0,			'v'},
      {"seq",		required_argument,  0,			'x'},
      {"seq1",		required_argument,  0,			'x'},
      {"seq2",		required_argument,  0,			'y'},
      {"ws",		required_argument,  0,			'w'},
      {"min",		required_argument,  0,			'l'},
      {"max",		required_argument,  0,			'm'},
      {"crop",		required_argument,  0,			'c'},
      {"threads",	required_argument,  0,			't'},
	  { "seed", required_argument, 0, 's' },
	  { "query", required_argument, 0, 'q' },
	  { "maxVariantFrac", required_argument, 0, 'a' },
	  { "minVariantDepth", required_argument, 0, 'b' },
	  { "index_folder", required_argument, 0, 'd' },
      {0,  0,  0, 0},
    };

  if (argc == 1)
  {
    Print_H();
    return 0;
  }
  
  while ( (o = getopt_long ( argc, argv, "hvn:e:o:u:i:s:x:y:w:l:m:c:a:b:d:g:p:r:s:t:q:d:", longOptions, &index)) != -1 )
    {
      switch (o)
	  {
	case 'd':
		folder_path = (char*)malloc(NAME_LENGTH);
		strcpy(folder_path, optarg);
		if (folder_path[strlen(folder_path)-1] == '/')
		{
			folder_path[strlen(folder_path) - 1] = 0;
		}
		break;
	case 'i':
	  is_index = 1;
	  fastaFile = (char*)malloc(NAME_LENGTH);
	  strcpy(fastaFile, optarg);
	  ///fastaFile = optarg;
	  break;
	case 'f':
		is_methy = 1;
		fastaFile = (char*)malloc(NAME_LENGTH);
		strcpy(fastaFile, optarg);
		///fastaFile = optarg;
		break;
	case 'g':
	  is_search = 1;
	  fastaFile = (char*)malloc(NAME_LENGTH);
	  strcpy(fastaFile, optarg);
	  ///fastaFile = optarg;
	  break;
	case 'c':
	  cropSize = atoi(optarg);
	  break;
	case 'w':
	  WINDOW_SIZE = atoi(optarg);
	  default_ws=0;
	  break;
	case 'x':
	  Read_File1 = optarg;
	  break;
	case 'q':
	  Read_File1 = optarg;
		break;
	case 'y':
	  Read_File2 = optarg;
	  break;
	case 'u':
	  Un_Mapped_File = optarg;
	  break;
	case 'o':
	  Mapped_File = (char*)malloc(NAME_LENGTH);
	  Mapped_FilePath = (char*)malloc(NAME_LENGTH);
	  stripPath (optarg, &Mapped_FilePath, &Mapped_File);
	  break;
	case 'n':
	  maxHits = atoi(optarg);
	  break;
	case 'e':
	  thread_e = 1;
	  thread_e_f = atof(optarg);
	  break;
	case 'a':
	  ///bs_available_seed_length = atoi(optarg);
	  maxVariantFrac = atof(optarg);
	  break;
	case 'b':
	  minVariantDepth = atoi(optarg);
	  break;
	case 'l':
	  minDistance_pair = atoi(optarg);
	  break;
	case 'm':
	  maxDistance_pair = atoi(optarg);
	  break;
	case 's':
	  over_all_seed_length = atoi(optarg);
	  break;
	case 'h':
	  Print_H();
	  return 0;
	  break;
	case 'v':
	  fprintf(stdout, "BitMapper %s\n", versionN);
	  return 0;
	  break;
    case 't':
    THREAD_COUNT = atoi(optarg);
    if (THREAD_COUNT == 0 || THREAD_COUNT > sysconf( _SC_NPROCESSORS_ONLN ))
        THREAD_COUNT = sysconf( _SC_NPROCESSORS_ONLN );
    break;
	case '?':
	  fprintf(stderr, "Unavailable parameter: %s\n", longOptions[index].name);
	  abort();
	  return 0;
	  break;
	}

    }

  if (is_index + is_search + is_methy!= 1)
    {
      fprintf(stdout, "Please select indexing or searching mode!\n");
      return 0;
    }

  if (WINDOW_SIZE > 16 || WINDOW_SIZE < 11)
    {
      ///fprintf(stdout, "Please set window size in [11..16]\n");
      ///return 0;
	  fprintf(stdout, "Warning: we strongly recommend users to set window size in [11..16]!\n");
    }


  if ( is_index )
    {
        if (fastaFile == NULL)
        {
          fprintf(stdout, "Please indicate the reference file for indexing!\n");
          return 0;
        }
    }
  else if (is_search)
  {


	  if (is_local == -1)
	  {
		  fprintf(stdout, "Please indicate the alignment mode. 1. local mode: --local. 2. end-to-end mode: --end-to-end.\n");
		  return 0;
	  }
	  


    if (fastaFile == NULL)
	{
	  fprintf(stdout, "Please indicate the reference file for searching!\n");
	  return 0;
	}

    if (Read_File1 == NULL && Read_File2 == NULL)
	{
	  fprintf(stdout, "Please indicate the read files for searching!\n");
	  return 0;
	}


      if (!is_pairedEnd && Read_File2 != NULL)
	{
	  fprintf(stdout, "Please not indicate the read2 files for single-end mode!\n");
	  return 0;
	}

      if (is_pairedEnd && (minDistance_pair <0 || maxDistance_pair < 0 || minDistance_pair > maxDistance_pair))
	{
	  fprintf(stdout, "Please set a valid range for paired-end mode.\n");
	  return 0;
	}

      if (is_pairedEnd && Read_File1 == NULL)
	{
	  fprintf(stdout, "Please indicate the read1 file for single-end mode.\n");
	  return 0;
	}


    }
  else
  {
	  if (fastaFile == NULL)
	  {
		  fprintf(stdout, "Please indicate the reference file for extracting!\n");
		  return 0;
	  }


	  if (Read_File1 == NULL)
	  {
		  fprintf(stdout, "Please indicate the alignment files for extracting!\n");
		  return 0;
	  }
  }

  if (unmapped_out == 1 || ambiguous_out == 1)
  {
	  if (output_methy == 1)
	  {
		  fprintf(stdout, "Only unique mapped reads can be processed by methylation extraction!\n");
		  fprintf(stdout, "Please do not use '--unmapped_out' or '--ambiguous_out'!\n");
		  return 0;
	  }
  }

  if (fastaFile[strlen(fastaFile)-1] == '/')
  {
	  fastaFile[strlen(fastaFile) - 1] = '\0';
  }
  if (is_index)
  {

	  if (folder_path != NULL)
	  {
		  char command[NAME_LENGTH];
		  sprintf(command, "mkdir %s", folder_path);
		  system(command);
		  sprintf(fileName[0], "%s", fastaFile);
		  ///sprintf(fileName[1], "%s/%s.index", folder_path, folder_path);
		  sprintf(fileName[1], "%s/genome.index", folder_path);
		 
	  }
	  else
	  {
		  sprintf(fileName[0], "%s", fastaFile);
		  sprintf(fileName[1], "%s.index", fileName[0]);
	  }

  }
  else
  {
	  sprintf(fileName[0], "%s", fastaFile);
	  sprintf(fileName[1], "%s.index", fileName[0]);
  }

 
  
  fprintf(stderr, "fileName[0]: %s\n", fileName[0]);
  fprintf(stderr, "fileName[1]: %s\n", fileName[1]);
  
  

  return 1;
}


void Print_H()
{
  char *errorType;


  fprintf(stdout,"BitMapperBS: a fast and accurate read aligner for whole-genome bisulte sequencing.\n\n");
  fprintf(stdout,"Usage: bitmapperBS [options]\n\n");
  errorType="edit distance";


  fprintf(stdout,"General Options:\n");
  fprintf(stdout," -v|--version\t\tCurrent Version.\n");
  fprintf(stdout," -h\t\t\tShow the help file.\n");
  fprintf(stdout,"\n\n");

  fprintf(stdout,"Options of indexing step:\n");
  ///fprintf(stdout, "Usage: bitmapperBS --index <ref.fa>\n\n");
  fprintf(stdout," --index [file]\t\tGenerate an index from the specified fasta file. \n");
  fprintf(stdout," --index_folder [folder]Set the folder that stores the genome indexes. If this option is not set, \n\t\t\tthe indexes would be stores in the same folder of genome (input fasta file). \n");
  fprintf(stdout,"\n\n");

  fprintf(stdout,"Options of read mapping step:\n");
  ///fprintf(stdout, "Usage: bitmapperBS --search <ref.fa> --seq1 --seq1\n\n");
  fprintf(stdout," --search [file/folder]\tSearch in the specified genome. If the indexes of this genome are built without \"--index_folder\", \n\t\t\tplease provide the path to the fasta file of genome (index files should be in the same folder). \n\t\t\tOtherwise please provide the path to the index folder (set by \"--index_folder\" during indexing).\n");
  fprintf(stdout," --fast \t\tSet bitmapperBS in fast mode (default). This option is only available in paired-end mode.\n");
  fprintf(stdout," --sensitive \t\tSet bitmapperBS in sensitive mode. This option is only available in paired-end mode.\n");
  fprintf(stdout," --pe \t\t\tSearch will be done in paired-end mode.\n");
  fprintf(stdout," --seq [file]\t\tInput sequences in fastq/fastq.gz format [file]. This option is used  \n\t\t\tfor single-end reads.\n");
  fprintf(stdout," --seq1 [file]\t\tInput sequences in fastq/fastq.gz format [file] (First file). \n\t\t\tUse this option to indicate the first file of \n\t\t\tpaired-end reads. \n");
  fprintf(stdout," --seq2 [file]\t\tInput sequences in fastq/fastq.gz format [file] (Second file). \n\t\t\tUse this option to indicate the second file of \n\t\t\tpaired-end reads.  \n");
  fprintf(stdout," -o [file]\t\tOutput of the mapped sequences in SAM or BAM format. The default is \"output\" in SAM format.\n");
  fprintf(stdout, " --sam \t\t\tOutput mapping results in SAM format (default).\n");
  fprintf(stdout, " --bam \t\t\tOutput mapping results in BAM format.\n");
  fprintf(stdout, " --methy_out \t\tOutput the intermediate methylation result files, instead of SAM or BAM files.\n");
  ///fprintf(stdout," -u [file]\t\tSave unmapped sequences in fasta/fastq format.\n");
  ///fprintf(stdout," --seqcomp \t\tIndicates that the input sequences are compressed (gz).\n");
  ///fprintf(stdout," --outcomp \t\tIndicates that output file should be compressed (gz).\n");
  fprintf(stdout," -e [float]\t\tSet the edit distance rate of read length. This value is between 0 and 1 (default: 0.08 = 8%% of read length).\n");
  //fprintf(stdout, " -s [double]\t\tbs_score_threshold (default 0.3 of the read length).\n");
  //fprintf(stdout, " -d [double]\t\tbs_edit_distance_threshold (default 0.1 of the read length).\n");
  ///fprintf(stdout," --min [int]\t\tMin distance allowed between a pair of end sequences (default: 0).\n");
  ///fprintf(stdout," --max [int]\t\tMax distance allowed between a pair of end sequences (default: 500).\n");
  fprintf(stdout, " --min [int]\t\tMin observed template length between a pair of end sequences (default: 0).\n");
  fprintf(stdout, " --max [int]\t\tMax observed template length between a pair of end sequences (default: 500).\n");
  ///fprintf(stdout," --crop [int]\t\tTrim the reads to the given length.\n");
  fprintf(stdout," --threads, -t [int]\tSet the number of CPU threads (default: 1).\n");
  ///fprintf(stdout, " --seed, -s [int]\tSet the length of seed.\n");
  fprintf(stdout, " --pbat \t\tMapping the BS-seq from pbat protocol.\n");
  fprintf(stdout, " --unmapped_out \tReport unmapped reads.\n");
  fprintf(stdout, " --ambiguous_out \tRandom report one of hit of each ambiguous mapped read.\n");
  fprintf(stdout, " --mapstats \t\tOutput the statistical information of read alignment into file named \"OUTPUT_FILE.mapstats\", \n\t\t\twhere \"OUTPUT_FILE\" is the name of output SAM or BAM file (defined by the option \"-o\").\n");
  
  fprintf(stdout,"\n\n");


  /**
  fprintf(stdout, "Options of methylation extracting step:\n");
  fprintf(stdout, " --methy_extract [file]\tExtract in the specified genome. Provide the path to the fasta file. \n\t\t\tIndex file should be in the same directory.\n");
  fprintf(stdout, " --seq [file]\t\tInput the intermediate methylation result files. Theses files must be generated by option '--methy_out' in read mapping step.\n\t\t\tFor example, if you use the option '-o output_methy' in read mapping step, here you must use the option '--seq output_methy'.\n");
  fprintf(stdout, " --CpG \t\t\tOutput CpG context methylation metrics (default).\n");
  fprintf(stdout, " --CHG \t\t\tOutput CHG context methylation metrics.\n");
  fprintf(stdout, " --CHH \t\t\tOutput CHH context methylation metrics.\n");
  fprintf(stdout, "\n\n");
  **/
  
}
