/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/







#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Process_CommandLines.h"
#include "Auxiliary.h"
#include "Process_Reads.h"
#include "Process_sam_out.h"
#include "Index.h"
#include "Schema.h"
#include "bam_prase.h"
#include "Levenshtein_Cal.h"



int main(int argc, char *argv[])
{

	_rg_name_l  *chhy_ih_refGenName;//存储各条染色体的信息
	bitmapper_bs_iter refChromeCont;//染色体条数
	pthread_t *_r_threads;
	int read_format;

  if (!CommandLine_process(argc, argv))
    return 1;
  

  if(is_index)
    {
      createIndex(fileName[0], fileName[1]);
    }
  else if (is_search)
    {


	  fprintf(stderr,"Read qualities are encoded by Phred+%d...\n", Q_base);
	  fprintf(stderr,"GapOpenPenalty: %d, GapExtensionPenalty: %d, MistMatchPenaltyMax: %d\n", 
	  GapOpenPenalty, GapExtensionPenalty, MistMatchPenaltyMax);
	  fprintf(stderr,"MistMatchPenaltyMin: %d, N_Penalty: %d\n", 
	  MistMatchPenaltyMin, N_Penalty);


      double totalLoadingTime = 0;
      double totalMappingTime = 0;
      double startTime;
      double loadingTime;
      double mappingTime;
      char outputFileName[NAME_LENGTH];
	  



      startTime = Get_T();



	  ///这就是打开read文件
	  if (!initiReadAllReads(Read_File1, Read_File2, is_pairedEnd, &read_format))
      {
		 fprintf(stderr, "Cannot open read files. \n");
         return 1;
      }

	  if(Mapped_File)
	  {
		  sprintf(outputFileName, "%s%s",Mapped_FilePath , Mapped_File);
	  }
	  else
	  {
		  outputFileName[0] = '\0';
	  }
	  

	  

        if (!is_pairedEnd)
        {

          	if (!Start_Load_Index(fileName[1]))
			{
				fprintf(stderr, "Load hash table index failed!\n");
				return 1;
			}


            mappingTime = 0;
            loadingTime = 0;

			

            fprintf(stderr,"Start load hash table!\n");
	
			Load_Index(thread_e, &chhy_ih_refGenName, &refChromeCont, fileName[1]);


            totalLoadingTime += Get_T()-startTime;
            fprintf(stderr, "Start alignment!\n");

			if (output_methy == 0)
			{
				if (bam_output == 0)
				{
					Output_gene(outputFileName);   
					OutPutSAM_Nounheader(chhy_ih_refGenName, refChromeCont, argc, argv);
				}
				else
				{
					init_bam_header(outputFileName, chhy_ih_refGenName, refChromeCont, argc, argv);
				}
				
			}
            
            

			Prepare_alignment(outputFileName, fileName[0], chhy_ih_refGenName, refChromeCont, read_format, is_pairedEnd);

			startTime = Get_T();


			if (THREAD_COUNT == 1)
			{
				if (pbat == 0)
				{
					Map_Single_Seq(0);
					///Map_Single_Seq_end_to_end_cover(0);
				}
				else
				{
					Map_Single_Seq_pbat(0);
				}
			}
            else
            {
				if (pbat == 0)
				{
					Map_Single_Seq_muti_thread(0);
				}
				else
				{
					Map_Single_Seq_pbat_muti_thread(0);
				}
            }


            totalMappingTime += Get_T()-startTime;

        }
        else
        {
			if (output_methy == 1 && methylation_size % 2 != 0)
			{
				fprintf(stderr, "The variable 'methylation_size' must be an even! Please change it to be even!\n");
				exit(1);
			}

			///_r_fp1
			if (pbat == 1)
			{
				exchange_two_reads();
			}

            minDistance_pair = minDistance_pair - SEQ_LENGTH;
            maxDistance_pair = maxDistance_pair - SEQ_LENGTH;


            if (!Start_Load_Index(fileName[1]))
            {
                fprintf(stderr, "Load hash table index failed!\n");
                return 1;
            }


            mappingTime = 0;
            loadingTime = 0;
             //loadHashTable = &Load_Index;
             fprintf(stderr,"Start load hash table!\n");
             // loadHashTable = &Load_Index;
			 Load_Index(thread_e, &chhy_ih_refGenName, &refChromeCont, fileName[1]);



			 

            totalLoadingTime += Get_T()-startTime;

			if (is_local == 1)
			{
				fprintf(stderr, "Start alignment in default fast mode.\n");
			}
			else
			{
				fprintf(stderr, "Start alignment in sensitive mode.\n");
			}

            

			if (output_methy == 0)
			{
				if (bam_output == 0)
				{
					Output_gene(outputFileName); 
					OutPutSAM_Nounheader(chhy_ih_refGenName, refChromeCont, argc, argv);
				}
				else
				{
					init_bam_header(outputFileName, chhy_ih_refGenName, refChromeCont, argc, argv);
				}
				
			}


            startTime = Get_T();

			Prepare_alignment(outputFileName, fileName[0], chhy_ih_refGenName, refChromeCont, read_format, is_pairedEnd);


             if(THREAD_COUNT==1)
            {

				Map_Pair_Seq(0);
            }
            else
            {
				Map_Pair_Seq_muti_thread(0);

            }



			 if (output_methy == 1)
			 {
				 ///out_paired_distance_statistic();
			 }


            totalMappingTime += Get_T()-startTime;

        }


		if (bam_output == 1)
		{
			close_bam_file();
		}

	  ///这个一打开速度就奇慢，不知道为什么
      ///finalizeOutput();

      fprintf(stderr, "-----------------------------------------------------------------------------------------------------------\n");
      fprintf(stderr, "%19s%16.2f%18.2f\n\n", "Total:",totalLoadingTime, totalMappingTime);

      fprintf(stderr, "%-42s%10.2f\n","Total Time:", totalMappingTime+totalLoadingTime);

	  long long number_of_read;
	  long long number_of_unique_mapped_read;
	  long long number_of_ambiguous_mapped_read;
	  long long number_of_unmapped_read;
	  long long number_of_mapped_bases;
	  long long number_of_mapped_errors;
	  get_mapping_informations(&number_of_read, &number_of_unique_mapped_read, &number_of_ambiguous_mapped_read,
		  &number_of_unmapped_read, &number_of_mapped_bases, &number_of_mapped_errors);

	  fprintf(stderr, "%-48s%lld\n", "No. of Reads:", number_of_read);
	  fprintf(stderr, "%-48s%lld (%0.2f\%)\n", "No. of Unique Mapped Reads:", 
		  number_of_unique_mapped_read, ((double)number_of_unique_mapped_read / (double)number_of_read)*100);
	  fprintf(stderr, "%-48s%lld (%0.2f\%)\n", "No. of Ambiguous Mapped Reads:", 
		  number_of_ambiguous_mapped_read, ((double)number_of_ambiguous_mapped_read / (double)number_of_read) * 100);
	  fprintf(stderr, "%-48s%lld (%0.2f\%)\n", "No. of Unmapped Reads:", 
		  number_of_unmapped_read, ((double)number_of_unmapped_read / (double)number_of_read) * 100);

	  fprintf(stderr, "%-47s %0.2f\%\n", "Mismatch and Indel Rate:",
		  ((double)number_of_mapped_errors / (double)number_of_mapped_bases) * 100);

	  if (mapstats == 1)
	  {
		  sprintf(Mapstats_File_Path + strlen(Mapstats_File_Path), "%s", Mapstats_File);
		  fprintf(stderr, "The statistical information will be written to %s ...\n", Mapstats_File_Path);
		
		  FILE* mapstats_fp = fopen(Mapstats_File_Path, "w");
			 
		  if (mapstats_fp != NULL)
		  {
			  fprintf(mapstats_fp, "%-48s%lld\n", "No. of Reads:", number_of_read);
			  fprintf(mapstats_fp, "%-48s%lld (%0.2f\%)\n", "No. of Unique Mapped Reads:",
				  number_of_unique_mapped_read, ((double)number_of_unique_mapped_read / (double)number_of_read) * 100);
			  fprintf(mapstats_fp, "%-48s%lld (%0.2f\%)\n", "No. of Ambiguous Mapped Reads:",
				  number_of_ambiguous_mapped_read, ((double)number_of_ambiguous_mapped_read / (double)number_of_read) * 100);
			  fprintf(mapstats_fp, "%-48s%lld (%0.2f\%)\n", "No. of Unmapped Reads:",
				  number_of_unmapped_read, ((double)number_of_unmapped_read / (double)number_of_read) * 100);

			  fprintf(mapstats_fp, "%-47s %0.2f\%\n", "Mismatch and Indel Rate:",
				  ((double)number_of_mapped_errors / (double)number_of_mapped_bases) * 100);

			  fclose(mapstats_fp);
		  }
	  }

    }
	else if (is_methy)
	{

		double total_time = 0;
		double startTime = Get_T();

		fprintf(stderr, "minVariantDepth: %d\n", minVariantDepth);
		fprintf(stderr, "maxVariantFrac: %f\n", maxVariantFrac);

		sprintf(fileName[1], "%s.methy", fileName[1]);

		if (!Load_Methy_Index(thread_e, &chhy_ih_refGenName, &refChromeCont, fileName[1]))
		{
			fprintf(stderr, "Load hash table index failed!\n");
			return 1;
		}

		///need_context[0]是CpG
		///need_context[1]是CHG
		///need_context[2]是CHH
		///int need_context[3] = { 1, 0, 0 };
		int need_context[3];
		need_context[0] = CpG;
		need_context[1] = CHG;
		need_context[2] = CHH;

		get_genome_cuts(Read_File1);

		fprintf(stderr, "genome_cuts: %d\n", genome_cuts);
		fprintf(stderr, "PE_distance: %d\n", maxDistance_pair);

		init_output_methy(Read_File1, need_context);

		Prepare_methy(fileName[0], chhy_ih_refGenName, refChromeCont);

		
		

		///int PE_distance = 500;
		
		if (!is_pairedEnd)
		{
			
			fprintf(stderr, "Extract from single-end alignment ...\n");

			if (THREAD_COUNT == 1)
			{
				methy_extract(0, Read_File1, need_context);
			}
			else
			{
				methy_extract_mutiple_thread(0, Read_File1, need_context);
			}
		}
		else
		{
			fprintf(stderr, "Extract from paired-end alignment ...\n");

			if (THREAD_COUNT == 1)
			{
				methy_extract_PE(0, Read_File1, maxDistance_pair, need_context);
			}
			else
			{
				methy_extract_PE_mutiple_thread(0, Read_File1, maxDistance_pair, need_context);
			}
			
		}


		total_time += Get_T() - startTime;

		fprintf(stderr, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stderr, "%19s%16.2f\n\n", "Total:", total_time);
	}


	

  return 0;
}
