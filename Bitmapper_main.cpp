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
#include "Levenshtein_Cal.h"



int main(int argc, char *argv[])
{

	_rg_name_l  *chhy_ih_refGenName;//存储各条染色体的信息
	bitmapper_bs_iter refChromeCont;//染色体条数
	pthread_t *_r_threads;
	int read_format;

  if (!CommandLine_process(argc, argv))
    return 1;



	  ///fprintf(stderr, "output_methy:%d\n", output_methy);
  

  if(is_index)
    {
      createIndex(fileName[0], fileName[1]);
    }
  else if (is_search)
    {




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
		 fprintf(stdout, "Cannot open read files. \n");
         return 1;
      }


	   

       sprintf(outputFileName, "%s%s",Mapped_FilePath , Mapped_File);

	   if (output_methy == 0)
	   {
		   Output_gene(outputFileName);
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
             //loadHashTable = &Load_Index;
            fprintf(stdout,"Start load hash table!\n");

			fprintf(stdout, "%d\n", READS_QUENUE_MAX_LENGTH);

             // loadHashTable = &Load_Index;
			Load_Index(thread_e, &chhy_ih_refGenName, &refChromeCont, fileName[1]);


            totalLoadingTime += Get_T()-startTime;
            fprintf(stdout, "Start alignment!\n");

			if (output_methy == 0)
			{
				//这个Prepare_alignment()就是将这个文件里的变量各种配置到新文件中
				OutPutSAM_Nounheader(chhy_ih_refGenName, refChromeCont, argc, argv);
			}
            
            

			Prepare_alignment(outputFileName, fileName[0], chhy_ih_refGenName, refChromeCont, read_format);

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
            fprintf(stdout, "sucess!\n");

        }
        else
        {
			if (output_methy == 1 && methylation_size % 2 != 0)
			{
				fprintf(stdout, "The variable 'methylation_size' must be an even! Please change it to be even!\n");
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
             fprintf(stdout,"Start load hash table!\n");
             // loadHashTable = &Load_Index;
			 Load_Index(thread_e, &chhy_ih_refGenName, &refChromeCont, fileName[1]);



			 

            totalLoadingTime += Get_T()-startTime;

			if (is_local == 1)
			{
				fprintf(stdout, "Start alignment in default fast mode.\n");
			}
			else
			{
				fprintf(stdout, "Start alignment in sensitive mode.\n");
			}


			///fprintf(stderr, "is_local: %d\n", is_local);
            

			if (output_methy == 0)
			{

				OutPutSAM_Nounheader(chhy_ih_refGenName, refChromeCont, argc, argv);
			}


            startTime = Get_T();

			Prepare_alignment(outputFileName, fileName[0], chhy_ih_refGenName, refChromeCont, read_format);


             if(THREAD_COUNT==1)
            {

				 Map_Pair_Seq(0);
            }
            else
            {
				Map_Pair_Seq_muti_thread(0);

            }



            totalMappingTime += Get_T()-startTime;
            fprintf(stdout, "sucess!\n");

        }

	  ///这个一打开速度就奇慢，不知道为什么
      ///finalizeOutput();

      fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
      fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:",totalLoadingTime, totalMappingTime);

      fprintf(stdout, "%-42s%10.2f\n","Total Time:", totalMappingTime+totalLoadingTime);

	  long long number_of_read;
	  long long number_of_unique_mapped_read;
	  long long number_of_ambiguous_mapped_read;
	  long long number_of_unmapped_read;
	  get_mapping_informations(&number_of_read, &number_of_unique_mapped_read, &number_of_ambiguous_mapped_read,
		  &number_of_unmapped_read);

	  fprintf(stdout, "%-48s%lld\n", "No. of Reads:", number_of_read);
	  fprintf(stdout, "%-48s%lld (%0.2f\%)\n", "No. of Unique Mapped Reads:", 
		  number_of_unique_mapped_read, ((double)number_of_unique_mapped_read / (double)number_of_read)*100);
	  fprintf(stdout, "%-48s%lld (%0.2f\%)\n", "No. of Ambiguous Mapped Reads:", 
		  number_of_ambiguous_mapped_read, ((double)number_of_ambiguous_mapped_read / (double)number_of_read) * 100);
	  fprintf(stdout, "%-48s%lld (%0.2f\%)\n", "No. of Unmapped Reads:", 
		  number_of_unmapped_read, ((double)number_of_unmapped_read / (double)number_of_read) * 100);

    }
	else if (is_methy)
	{

		sprintf(fileName[1], "%s.methy", fileName[1]);

		if (!Load_Methy_Index(thread_e, &chhy_ih_refGenName, &refChromeCont, fileName[1]))
		{
			fprintf(stderr, "Load hash table index failed!\n");
			return 1;
		}


		Prepare_methy(fileName[0], chhy_ih_refGenName, refChromeCont);

		methy_extract(0, Read_File1);
	}


  return 0;
}
