/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<stdint.h>
#include <xmmintrin.h>
#include <mmintrin.h>
#include "Levenshtein_Cal.h"


unsigned short tmp_Peq[MAX_Thread][128][8];
unsigned int tmp_Peq_32[MAX_Thread][128][4];
int each_length;

int Calcu_Cigar_MD(int Peq_i,char *pattern,int p_length,char *text,int t_length,
                             unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int err,int match_site,char* cigar,char* MD_Z,int thread_id)
{

    int mis_sites[SEQ_MAX_LENGTH];

    int i = 0;

    ///这里要添加

    int start_location=match_site-t_length+1;

    int tmp_err=0;

    for(i=0;i<t_length;i++)
    {
        if(text[i]!=pattern[i+start_location])
        {
            mis_sites[tmp_err]=i;
            tmp_err++;
        }
    }

    if(tmp_err==err)
    {
        i=0;

        int pre_i=-1;


        for(i=0;i<tmp_err;i++)
        {
            if(mis_sites[i]-pre_i>1)
            {
                sprintf(MD_Z+strlen(MD_Z),"%d%c",mis_sites[i]-pre_i-1,pattern[mis_sites[i]+start_location]);
            }
            else
            {
                sprintf(MD_Z+strlen(MD_Z),"%c",pattern[mis_sites[i]+start_location]);
            }

            pre_i=mis_sites[i];
        }
        if(pre_i+1<t_length)
        {
            sprintf(MD_Z+strlen(MD_Z),"%d",t_length-pre_i-1);
        }

        sprintf(cigar+strlen(cigar),"%d%c",t_length,'M');
        return start_location;
    }

    start_location=match_site-t_length+1;
    Word D0_arry_64[2000];
    Word HP_arry_64[2000];
    int Route_Size_Whole[2000];
    char Route_Char_Whole[2000];
    char err_match[2000];

    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Word tmp_Peq_1=(Word)1;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }

    Peq['A']=tmp_Peq[thread_id]['A'][Peq_i];
    Peq['C']=tmp_Peq[thread_id]['C'][Peq_i];
    Peq['T']=tmp_Peq[thread_id]['T'][Peq_i];
    Peq['G']=tmp_Peq[thread_id]['G'][Peq_i];
    Word Mask_Pre=(Word)1<<(band_length-2);
    Word Mask=(Word)1<<(band_length-1);
    Word VP=0;
    Word VN=0;
    Word X=0;
    Word HN=0;
    int s=0;
    int j=0;
    i=0;
    int bound=band_length-2-band_down;
    Word err_mask=(Word)1;
    int s1=band_length-2;
    int i_bd=i+band_down;
    int last_high=band_length-t_length+p_length-band_down-1;

    Word bit_Mask=((Word)-1)>>(64-band_length);
    Word tmp_D0,tmp_HP;
    tmp_D0=0;
    tmp_HP=0;
    int t_length_1=t_length-1;
    while(i<t_length_1)
    {
        X=Peq[text[i]]|VN;
        tmp_D0=((VP+(X&VP))^VP)|X;
        HN=VP&tmp_D0;
        tmp_HP=VN|~(VP|tmp_D0);
        X=tmp_D0>>1;
        VN=X&tmp_HP;
        VP=HN|~(X|tmp_HP);
        D0_arry_64[i]=tmp_D0;
        HP_arry_64[i]=tmp_HP;

        Peq['A']=Peq['A']>>1;
        Peq['C']=Peq['C']>>1;
        Peq['G']=Peq['G']>>1;
        Peq['T']=Peq['T']>>1;


        ++i;
        ++i_bd;
        if((i_bd)<p_length)
            Peq[pattern[i_bd]]=Peq[pattern[i_bd]]|Mask;
    }




        X=Peq[text[i]]|VN;
        tmp_D0=((VP+(X&VP))^VP)|X;
        HN=VP&tmp_D0;
        tmp_HP=VN|~(VP|tmp_D0);
        X=tmp_D0>>1;
        VN=X&tmp_HP;
        VP=HN|~(X|tmp_HP);
        D0_arry_64[i]=tmp_D0;
        HP_arry_64[i]=tmp_HP;







    int site=p_length-last_high-1;
    ///对角方向和竖向就是这个；要是横向还要减1
    int search_site=match_site-site;
    ///要是竖向还要减1
    int pre_size=1;
    char pre_char='N';
    Word Mask_1=(Word)1;
    i=t_length-1;
    int sum_err=0;

    ///Route_size[0]存储到是Route_size这个数组到长度
    j=1;
    int err_char_i=0;


    ///这个放最前面是因为这些都是匹配的片段，所以绝大部分时间是匹配到
    if(((D0_arry_64[i]>>search_site)&Mask_1)&&(pattern[match_site]==text[i]))   ///对角方向(匹配)
    {
        i--;
        pre_size=1;
        pre_char='M';
        match_site--;

        ///j--;
    }
    else if(!((D0_arry_64[i]>>search_site)&Mask_1)&&(pattern[match_site]!=text[i]))                         ///对角方向(替换)
    {
        err_match[err_char_i]=pattern[match_site];
        err_char_i++;

        i--;
        pre_size=1;
        pre_char='S';
        match_site--;
        sum_err++;




        ///j--;

    }
    else if((HP_arry_64[i]>>search_site)&Mask_1)    ///水平方向过来,这个表示基因组上缺了一块  I代表参考基因组比read上少，即参考基因组比read缺了一块
    {

        err_match[err_char_i]=pattern[match_site];
        err_char_i++;



        i--;
        search_site++;
        pre_size=1;

        ///'I'代表基因组上缺了一块，但是这里是最后一个位置
        ///在这个地方，'S'与'I'是一样的
        pre_char='S';
        sum_err++;
        start_location++;


    }
    else  ///上方过来到    ///这个属于read上缺了一块
    {
        search_site--;
        pre_size=1;
        pre_char='D';

        err_match[err_char_i]=pattern[match_site];;
        err_char_i++;

        match_site--;
        sum_err++;
        start_location--;

    }

    while(i>=0)
    {

        if(sum_err==err)
        {
            break;
        }

        ///这个放最前面是因为这些都是匹配的片段，所以绝大部分时间是匹配到
        if(((D0_arry_64[i]>>search_site)&Mask_1)&&(pattern[match_site]==text[i]))   ///对角方向(匹配)
        {
            i--;
            match_site--;

            if(pre_char!='M')
            {
                Route_Size_Whole[j]=pre_size;
                Route_Char_Whole[j++]=pre_char;
                pre_size=1;
                pre_char='M';

            }
            else
            {
                pre_size++;
            }

        }
        else if(!((D0_arry_64[i]>>search_site)&Mask_1)&&(pattern[match_site]!=text[i]))                        ///对角方向(替换)
        {
            err_match[err_char_i]=pattern[match_site];
            err_char_i++;

            i--;
            match_site--;
            sum_err++;


            if(pre_char!='S')
            {

                Route_Size_Whole[j]=pre_size;
                Route_Char_Whole[j++]=pre_char;
                pre_size=1;
                pre_char='S';

            }
            else
            {
                pre_size++;
            }

        }
        else if((HP_arry_64[i]>>search_site)&Mask_1)  ///水平方向过来,这个表示基因组上缺了一块  I代表参考基因组比read上少，即参考基因组比read缺了一块
        {
            i--;
            search_site++;
            sum_err++;
            if(pre_char!='I')
            {
                Route_Size_Whole[j]=pre_size;
                Route_Char_Whole[j++]=pre_char;
                pre_size=1;
                pre_char='I';

            }
            else
            {
                pre_size++;
            }
            start_location++;
        }
        else  ///上方过来到
        {
            err_match[err_char_i]=pattern[match_site];;
            err_char_i++;
            search_site--;
            match_site--;
            sum_err++;
            if(pre_char!='D')
            {
                Route_Size_Whole[j]=pre_size;
                Route_Char_Whole[j++]=pre_char;
                pre_size=1;
                pre_char='D';

            }
            else
            {
                pre_size++;
            }
            start_location--;
        }




    }

    Route_Size_Whole[j]=pre_size;
    Route_Char_Whole[j++]=pre_char;
    ///如果是break出来的，也就是后面还有M，则i一定是-1；i只要大于等于0，就代表提前终止了
    if(i>=0)
    {
        Route_Size_Whole[j]=i+1;
        Route_Char_Whole[j++]='M';
    }

    Route_Size_Whole[0]=j-1;
    char pre_cigar=0;
    ///这里要修改
    err_char_i--;




    int ijk=0;
    int is_D=0;
    int pre_length=0;
    for(j--;j>=1;--j)
    {
        if(Route_Char_Whole[j]=='M')
        {
            pre_length=pre_length+Route_Size_Whole[j];
            ///sprintf(MD_Z+strlen(MD_Z),"%d",Route_Size_Whole[j]);
            is_D=0;
        }
        else if(Route_Char_Whole[j]=='S')
        {
            if(pre_length!=0)
                sprintf(MD_Z+strlen(MD_Z),"%d",pre_length);
            pre_length=0;
            if(is_D)
            {
                sprintf(MD_Z+strlen(MD_Z),"%d",0);
            }
            for(ijk=0;ijk<Route_Size_Whole[j];++ijk)
            {
                sprintf(MD_Z+strlen(MD_Z),"%c",err_match[err_char_i]);
                err_char_i--;      ///这里要改
            }
            is_D=0;
        }
        else if(Route_Char_Whole[j]=='D')   ///这里要改
        {
             if(pre_length!=0)
                sprintf(MD_Z+strlen(MD_Z),"%d",pre_length);
            pre_length=0;

            sprintf(MD_Z+strlen(MD_Z),"%c",'^');
            for(ijk=0;ijk<Route_Size_Whole[j];++ijk)
            {
                sprintf(MD_Z+strlen(MD_Z),"%c",err_match[err_char_i]);
                err_char_i--;       ///这里要改
            }
            is_D=1;
        }
        else
        {
            is_D=0;
        }
    }

    if(pre_length!=0)
        sprintf(MD_Z+strlen(MD_Z),"%d",pre_length);


    int  size_SM=0;
    j=Route_Size_Whole[0];
    for(;j>=1;--j)
    {
        if(Route_Char_Whole[j]=='M'||Route_Char_Whole[j]=='S')
        {
            size_SM=0;
             while(Route_Char_Whole[j]=='M'||Route_Char_Whole[j]=='S')
            {
                size_SM=size_SM+Route_Size_Whole[j];
                j--;
            }
            j++;
            sprintf(cigar+strlen(cigar),"%d%c",size_SM,'M');

        }
        else
        {
            sprintf(cigar+strlen(cigar),"%d%c",Route_Size_Whole[j],Route_Char_Whole[j]);
        }
    }









    return start_location;

}




int Start_location_Calcu_Cigar_MD(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int pre_err, char* cigar ,char* MD_Z,int thread_id,int return_site)
{
    ///这个与反向不同，已经知道start_location了
    int mis_sites[SEQ_MAX_LENGTH];
    int i = 0;

    ///基因组传过来的就已经是起始位置-errthold
    int tmp_err=0;
    for(i=0;i<t_length;i++)
    {
        if(text[i]!=pattern[thread_e+i])
        {
            mis_sites[tmp_err]=thread_e+i;
            tmp_err++;
        }
    }

    if(tmp_err==pre_err)
    {
        i=0;

        int pre_i=thread_e-1;


        for(i=0;i<tmp_err;i++)
        {
            if(mis_sites[i]-pre_i>1)
            {
                sprintf(MD_Z+strlen(MD_Z),"%d%c",mis_sites[i]-pre_i-1,pattern[mis_sites[i]]);
            }
            else
            {
                sprintf(MD_Z+strlen(MD_Z),"%c",pattern[mis_sites[i]]);
            }

            pre_i=mis_sites[i];
        }
        if(pre_i+1+thread_e<p_length)
        {
            sprintf(MD_Z+strlen(MD_Z),"%d",p_length-pre_i-1-thread_e);
        }
        sprintf(cigar+strlen(cigar),"%d%c",t_length,'M');

        return 1;
    }









    Word D0_arry_64[2000];
    Word HP_arry_64[2000];
    int Route_Size_Whole[2000];
    char Route_Char_Whole[2000];
    char err_match[2000];


    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Word tmp_Peq_1=(Word)1;
    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;
    for (r =0; r<band_length; r++)
    {
        Peq[pattern[p_length-r-1]]=Peq[pattern[p_length-r-1]]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;
    Word Mask_Pre=(Word)1<<(band_length-2);
    Word Mask=(Word)1<<(band_length-1);
    Word VP=0;
    Word VN=0;
    Word X=0;
    Word D0=0;

    Word HN=0;
    Word HP=0;
    int s=0;
    i = 0;
    int j=0;
    int bound=band_length-2-band_down;
    int err=0;
    Word err_mask=(Word)1;
    int s1=band_length-2;
    i=t_length-1;
    int i_bd=p_length-band_length;
    int last_high=band_length-t_length+p_length-band_down-1;
    while(i>0)
    {
        X=Peq[text[i]]|VN;
        D0=((VP+(X&VP))^VP)|X;

        D0_arry_64[i]=D0;

        HN=VP&D0;
        HP=VN|~(VP|D0);


        HP_arry_64[i]=HP;

        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);


        Peq['A']=Peq['A']>>1;
        Peq['T']=Peq['T']>>1;
        Peq['G']=Peq['G']>>1;
        Peq['C']=Peq['C']>>1;


        --i;
        --i_bd;
        Peq[pattern[i_bd]]=Peq[pattern[i_bd]]|Mask;
    }


        X=Peq[text[i]]|VN;
        D0=((VP+(X&VP))^VP)|X;

        D0_arry_64[i]=D0;

        HN=VP&D0;
        HP=VN|~(VP|D0);


        HP_arry_64[i]=HP;

        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);




    ///int match_site=errthold;
    ///应该左移才对
    int match_site=thread_e;
    ///int site=last_high-1;
    ///对角方向和竖向就是这个；要是横向还要减1
    ///int search_site=errthold;
    ///参考基因组的低位对应D0的高位
    int search_site=2*thread_e-match_site;
    ///要是竖向还要减1
    int pre_size=1;
    char pre_char='N';
    Word Mask_1=(Word)1;
    i=0;


    int sum_err=0;

    int err_char_i=0;
    char err_match_char[1000];
    err_match_char[0]='\0';
    int err_match_length[100];
    err_match_length[0]=0;


    char MD_Z_match[1000];
    int MD_Z_char_i=0;


    while(i<t_length)
    {


        if(sum_err==pre_err)
        {
            break;
        }


        if(((D0_arry_64[i]>>search_site)&Mask_1)&&(pattern[match_site]==text[i]))   ///对角方向(匹配)
        {

            i++;
            match_site++;
            if(err_match_char[err_char_i]!='M')
            {
                err_char_i++;
                err_match_char[err_char_i]='M';
                err_match_length[err_char_i]=1;
                err_match_length[0]++;
            }
            else
            {
                err_match_length[err_char_i]++;
            }
            ///j--;
        }
        else if(!((D0_arry_64[i]>>search_site)&Mask_1)&&(pattern[match_site]!=text[i]))                         ///对角方向(替换)
        {

            MD_Z_match[MD_Z_char_i]=pattern[match_site];
            MD_Z_char_i++;

            i++;
            match_site++;
            sum_err++;



            if(err_match_char[err_char_i]!='S')
            {
                err_char_i++;
                err_match_char[err_char_i]='S';
                err_match_length[err_char_i]=1;
                err_match_length[0]++;
            }
            else
            {
                err_match_length[err_char_i]++;
            }

        }
        else if((HP_arry_64[i]>>search_site)&Mask_1)    ///水平方向过来
        {

            i++;


            search_site++;


            sum_err++;


            if(err_match_char[err_char_i]!='I')
            {
                err_char_i++;
                err_match_char[err_char_i]='I';
                err_match_length[err_char_i]=1;
                err_match_length[0]++;
            }
            else
            {
                err_match_length[err_char_i]++;
            }

        }
        else  ///上方过来到
        {


            MD_Z_match[MD_Z_char_i]=pattern[match_site];
            MD_Z_char_i++;

            search_site--;
            match_site++;
            sum_err++;




            if(err_match_char[err_char_i]!='D')
            {
                err_char_i++;
                err_match_char[err_char_i]='D';
                err_match_length[err_char_i]=1;
                err_match_length[0]++;
            }
            else
            {
                err_match_length[err_char_i]++;
            }

            continue;

        }

    }


     if(i<t_length)
    {
        if(err_match_char[err_char_i]=='M')
        {
            err_match_length[err_char_i]+=t_length-i;
        }
        else
        {
            err_char_i++;
            err_match_char[err_char_i]='M';
            err_match_length[err_char_i]=t_length-i;
            err_match_length[0]++;
        }

    }





    int ijk=0;
    int is_D=0;
    int pre_length=0;
    MD_Z_char_i=0;
    for(i=1;i<=err_char_i;i++)
    {
        if(err_match_char[i]=='M')
        {
            pre_length=pre_length+err_match_length[i];
            ///sprintf(MD_Z+strlen(MD_Z),"%d",Route_Size_Whole[j]);
            is_D=0;
        }
        else if(err_match_char[i]=='S')
        {
            if(pre_length!=0)
                sprintf(MD_Z+strlen(MD_Z),"%d",pre_length);
            pre_length=0;
            if(is_D)
            {
                sprintf(MD_Z+strlen(MD_Z),"%d",0);
            }
            for(ijk=0;ijk<err_match_length[i];++ijk)
            {
                sprintf(MD_Z+strlen(MD_Z),"%c",MD_Z_match[MD_Z_char_i]);
                MD_Z_char_i++;      ///这里要改
            }
            is_D=0;
        }
        else if(err_match_char[i]=='D')   ///这里要改
        {
             if(pre_length!=0)
                sprintf(MD_Z+strlen(MD_Z),"%d",pre_length);
            pre_length=0;

            sprintf(MD_Z+strlen(MD_Z),"%c",'^');
            for(ijk=0;ijk<err_match_length[i];++ijk)
            {
                sprintf(MD_Z+strlen(MD_Z),"%c",MD_Z_match[MD_Z_char_i]);
                MD_Z_char_i++;       ///这里要改
            }
            is_D=1;
        }
        else
        {
            ///pre_length=pre_length+Route_Size_Whole[j];
            is_D=0;
        }

        ///fprintf(_out_fp,"%d%c",map.CIGAR_SIZE[j],map.CIGAR_CHAR[j]);
    }

    if(pre_length!=0)
        sprintf(MD_Z+strlen(MD_Z),"%d",pre_length);
    int  size_SM=0;

    /**
    for(i=1;i<=err_char_i;i++)
    {

            sprintf(cigar+strlen(cigar),"%d%c",err_match_length[i],err_match_char[i]);


    }

    sprintf(cigar+strlen(cigar),"#####");
    **/

    for(i=1;i<=err_char_i;i++)
    {
        if(err_match_char[i]=='M'||err_match_char[i]=='S')
        {
            size_SM=0;
            while(err_match_char[i]=='M'||err_match_char[i]=='S')
            {
                if(i>err_char_i)
                    break;
                size_SM=size_SM+err_match_length[i];

                ///sprintf(cigar+strlen(cigar),"i=%d,length=%d,char=%c",i,err_match_length[i],err_match_char[i]);


                i++;
            }
            i--;
            sprintf(cigar+strlen(cigar),"%d%c",size_SM,'M');

        }
        else
        {
            sprintf(cigar+strlen(cigar),"%d%c",err_match_length[i],err_match_char[i]);
        }


    }

}



/****************************************4条，固定区间长度的算法***************************************************************/
int Brief_4_Banded_BPM_Non_SSE(char *pattern1,char *pattern2,char* pattern3, char* pattern4, int p_length,char *text,int t_length,
                              int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    Word tmp_Peq_1=(Word)1;
    Word tmp_Peq_2=tmp_Peq_1<<(each_length);    ///注意each_length对应的是band_length+1而不是band_length
    Word tmp_Peq_3=tmp_Peq_2<<(each_length);
    Word tmp_Peq_4=tmp_Peq_3<<(each_length);

    int i=0;

    for (r =0; r<band_length; r++)
    {
        Peq[pattern1[r]]=Peq[pattern1[r]]|tmp_Peq_1;
        Peq[pattern2[r]]=Peq[pattern2[r]]|tmp_Peq_2;
        Peq[pattern3[r]]=Peq[pattern3[r]]|tmp_Peq_3;
        Peq[pattern4[r]]=Peq[pattern4[r]]|tmp_Peq_4;


        tmp_Peq_1=tmp_Peq_1<<1;
        tmp_Peq_2=tmp_Peq_2<<1;
        tmp_Peq_3=tmp_Peq_3<<1;
        tmp_Peq_4=tmp_Peq_4<<1;
    }


    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;

    Word All_1=(Word)-1;
    Word single_err_mask1=All_1>>(64-each_length);
    Word single_err_mask2=single_err_mask1<<(each_length);
    Word single_err_mask3=single_err_mask2<<(each_length);
    Word single_err_mask4=single_err_mask3<<(each_length);

    if(each_length<=16)
    {
        tmp_Peq[thread_id]['A'][0]=Peq_A&single_err_mask1;
        tmp_Peq[thread_id]['A'][1]=(Peq_A&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['A'][2]=(Peq_A&single_err_mask3)>>(each_length+each_length);
        tmp_Peq[thread_id]['A'][3]=(Peq_A)>>(each_length+each_length+each_length);

        tmp_Peq[thread_id]['T'][0]=Peq_T&single_err_mask1;
        tmp_Peq[thread_id]['T'][1]=(Peq_T&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['T'][2]=(Peq_T&single_err_mask3)>>(each_length+each_length);
        tmp_Peq[thread_id]['T'][3]=(Peq_T)>>(each_length+each_length+each_length);


        tmp_Peq[thread_id]['G'][0]=Peq_G&single_err_mask1;
        tmp_Peq[thread_id]['G'][1]=(Peq_G&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['G'][2]=(Peq_G&single_err_mask3)>>(each_length+each_length);
        tmp_Peq[thread_id]['G'][3]=(Peq_G)>>(each_length+each_length+each_length);


        tmp_Peq[thread_id]['C'][0]=Peq_C&single_err_mask1;
        tmp_Peq[thread_id]['C'][1]=(Peq_C&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['C'][2]=(Peq_C&single_err_mask3)>>(each_length+each_length);
        tmp_Peq[thread_id]['C'][3]=(Peq_C)>>(each_length+each_length+each_length);
    }
    else
    {

        tmp_Peq_32[thread_id]['A'][0]=Peq_A&single_err_mask1;
        tmp_Peq_32[thread_id]['A'][1]=(Peq_A&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['A'][2]=(Peq_A&single_err_mask3)>>(each_length+each_length);
        tmp_Peq_32[thread_id]['A'][3]=(Peq_A)>>(each_length+each_length+each_length);

        tmp_Peq_32[thread_id]['T'][0]=Peq_T&single_err_mask1;
        tmp_Peq_32[thread_id]['T'][1]=(Peq_T&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['T'][2]=(Peq_T&single_err_mask3)>>(each_length+each_length);
        tmp_Peq_32[thread_id]['T'][3]=(Peq_T)>>(each_length+each_length+each_length);


        tmp_Peq_32[thread_id]['G'][0]=Peq_G&single_err_mask1;
        tmp_Peq_32[thread_id]['G'][1]=(Peq_G&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['G'][2]=(Peq_G&single_err_mask3)>>(each_length+each_length);
        tmp_Peq_32[thread_id]['G'][3]=(Peq_G)>>(each_length+each_length+each_length);


        tmp_Peq_32[thread_id]['C'][0]=Peq_C&single_err_mask1;
        tmp_Peq_32[thread_id]['C'][1]=(Peq_C&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['C'][2]=(Peq_C&single_err_mask3)>>(each_length+each_length);
        tmp_Peq_32[thread_id]['C'][3]=(Peq_C)>>(each_length+each_length+each_length);

    }



    Word Mask1=(Word)1<<(band_length-1);
    Word Mask2=Mask1<<(each_length);    ///each_length对应的是band_length+1
    Word Mask3=Mask2<<(each_length);
    Word Mask4=Mask3<<(each_length);




    /*******************************这个是D0右移用的，好像也是加法也能用用的**************************************/
    ///这个是低band_length为1，高each_length-band_length为0
    Word  tmp_MASK=All_1>>(64-band_length);
    Word  Hppp_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    Word _Hppp_mask=~Hppp_mask;


    ///这个是循环右移的掩码
    tmp_MASK=All_1>>(64-band_length+1);
    Word Hpp_Mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;



    ///下面这两个变量是计算err用的
    ///这个处理方法必须要求需要计算编辑距离的区域的长度不能大于128
    ///计算err用的掩码就是低位全1
    tmp_MASK=(Word)1;
    Word err_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;



    Word VP=(Word)0;
    Word VN=(Word)0;
    Word X=(Word)0;
    Word D0=(Word)0;
    Word HN=(Word)0;
    Word HP=(Word)0;

    Word _TMP;

    Word TMP;


    Word Err=(Word)0;


    int s=0;
    i = 0;
    int j=0;
    ///int err[16];
    int err1;
    int err2;
    int err3;
    int err4;

    int i_bd=i+band_down;


    int last_high=band_length-t_length+p_length-band_down-1;

    int pre_end1=last_high+errthold;
    int pre_end2=pre_end1<<each_length;
    int pre_end3=pre_end2<<each_length;
    int pre_end4=pre_end3<<each_length;


    Word err_arry;


    /**运算**/
    int t_length_1=t_length-1;
    while(i<t_length_1)
    {


        X=Peq[text[i]]|VN;

        TMP=X&VP;


        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;

        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;

        D0=(TMP^VP)|X;

        HN=VP&D0;

        HP=VN|~(VP|D0);

        X=(D0>>1)&Hppp_mask;


        VN=X&HP;

        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;
        err4=Err&single_err_mask4;

        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3)&&(err4>pre_end4))
            return -1;

        Peq['A']=(Peq['A']>>1)&Hpp_Mask;
        Peq['T']=(Peq['T']>>1)&Hpp_Mask;
        Peq['G']=(Peq['G']>>1)&Hpp_Mask;
        Peq['C']=(Peq['C']>>1)&Hpp_Mask;


        ++i;
        ++i_bd;

        Peq[pattern1[i_bd]]=Peq[pattern1[i_bd]]|Mask1;
        Peq[pattern2[i_bd]]=Peq[pattern2[i_bd]]|Mask2;
        Peq[pattern3[i_bd]]=Peq[pattern3[i_bd]]|Mask3;
        Peq[pattern4[i_bd]]=Peq[pattern4[i_bd]]|Mask4;


    }


      X=Peq[text[i]]|VN;

        TMP=X&VP;


        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;

        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;

        D0=(TMP^VP)|X;

        HN=VP&D0;

        HP=VN|~(VP|D0);

        X=(D0>>1)&Hppp_mask;


        VN=X&HP;

        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;
        err4=Err&single_err_mask4;

        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3)&&(err4>pre_end4))
            return -1;



    int site=p_length-last_high-1;

    err1=Err&single_err_mask1;
    err2=(Err&single_err_mask2)>>each_length;
    err3=(Err&single_err_mask3)>>(each_length+each_length);
    err4=(Err&single_err_mask4)>>(each_length+each_length+each_length);


     if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }



    i=0;
    Word tmp_VP;


    while(i<last_high)
    {


        tmp_VP=VP>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err+err_arry;



        tmp_VP=VN>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err-err_arry;


        ++i;


        err1=Err&single_err_mask1;
        err2=(Err&single_err_mask2)>>each_length;
        err3=(Err&single_err_mask3)>>(each_length+each_length);
        err4=(Err&single_err_mask4)>>(each_length+each_length+each_length);


         if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site+i;
                return_sites_error[3]=err4;
        }

    }


    return 1;


}






/****************************************3条，固定区间长度的算法***************************************************************/
int Brief_3_Banded_BPM_Non_SSE(char *pattern1,char *pattern2,char* pattern3, int p_length,char *text,int t_length,
                              int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    ///这两个纯粹是为了计算peq用的
    Word tmp_Peq_1=(Word)1;
    Word tmp_Peq_2=tmp_Peq_1<<(each_length);    ///注意each_length对应的是band_length+1而不是band_length
    Word tmp_Peq_3=tmp_Peq_2<<(each_length);

    int i=0;


    for (r =0; r<band_length; r++)
    {
        Peq[pattern1[r]]=Peq[pattern1[r]]|tmp_Peq_1;
        Peq[pattern2[r]]=Peq[pattern2[r]]|tmp_Peq_2;
        Peq[pattern3[r]]=Peq[pattern3[r]]|tmp_Peq_3;


        tmp_Peq_1=tmp_Peq_1<<1;
        tmp_Peq_2=tmp_Peq_2<<1;
        tmp_Peq_3=tmp_Peq_3<<1;
    }


    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;

    Word All_1=(Word)-1;
    Word single_err_mask1=All_1>>(64-each_length);
    Word single_err_mask2=single_err_mask1<<(each_length);
    Word single_err_mask3=single_err_mask2<<(each_length);

    if(each_length<=16)
    {
        tmp_Peq[thread_id]['A'][0]=Peq_A&single_err_mask1;
        tmp_Peq[thread_id]['A'][1]=(Peq_A&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['A'][2]=(Peq_A)>>(each_length+each_length);

        tmp_Peq[thread_id]['T'][0]=Peq_T&single_err_mask1;
        tmp_Peq[thread_id]['T'][1]=(Peq_T&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['T'][2]=(Peq_T)>>(each_length+each_length);

        tmp_Peq[thread_id]['G'][0]=Peq_G&single_err_mask1;
        tmp_Peq[thread_id]['G'][1]=(Peq_G&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['G'][2]=(Peq_G)>>(each_length+each_length);


        tmp_Peq[thread_id]['C'][0]=Peq_C&single_err_mask1;
        tmp_Peq[thread_id]['C'][1]=(Peq_C&single_err_mask2)>>each_length;
        tmp_Peq[thread_id]['C'][2]=(Peq_C)>>(each_length+each_length);
    }
    else
    {

        tmp_Peq_32[thread_id]['A'][0]=Peq_A&single_err_mask1;
        tmp_Peq_32[thread_id]['A'][1]=(Peq_A&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['A'][2]=(Peq_A)>>(each_length+each_length);

        tmp_Peq_32[thread_id]['T'][0]=Peq_T&single_err_mask1;
        tmp_Peq_32[thread_id]['T'][1]=(Peq_T&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['T'][2]=(Peq_T)>>(each_length+each_length);

        tmp_Peq_32[thread_id]['G'][0]=Peq_G&single_err_mask1;
        tmp_Peq_32[thread_id]['G'][1]=(Peq_G&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['G'][2]=(Peq_G)>>(each_length+each_length);


        tmp_Peq_32[thread_id]['C'][0]=Peq_C&single_err_mask1;
        tmp_Peq_32[thread_id]['C'][1]=(Peq_C&single_err_mask2)>>each_length;
        tmp_Peq_32[thread_id]['C'][2]=(Peq_C)>>(each_length+each_length);

    }



    Word Mask1=(Word)1<<(band_length-1);
    Word Mask2=Mask1<<(each_length);    ///each_length对应的是band_length+1
    Word Mask3=Mask2<<(each_length);




    /*******************************这个是D0右移用的，好像也是加法也能用用的**************************************/
    ///这个是低band_length为1，高each_length-band_length为0
    Word  tmp_MASK=All_1>>(64-band_length);
    Word  Hppp_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    Word _Hppp_mask=~Hppp_mask;


    ///这个是循环右移的掩码
    tmp_MASK=All_1>>(64-band_length+1);
    Word Hpp_Mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;



    ///下面这两个变量是计算err用的
    ///这个处理方法必须要求需要计算编辑距离的区域的长度不能大于128
    ///计算err用的掩码就是低位全1
    tmp_MASK=(Word)1;
    Word err_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;



    Word VP=(Word)0;
    Word VN=(Word)0;
    Word X=(Word)0;
    Word D0=(Word)0;
    Word HN=(Word)0;
    Word HP=(Word)0;

    Word _TMP;

    Word TMP;

    Word Err=(Word)0;


    int s=0;

    i = 0;
    int j=0;
    ///int err[16];
    int err1;
    int err2;
    int err3;

    int i_bd=i+band_down;


    int last_high=band_length-t_length+p_length-band_down-1;



    int pre_end1=last_high+errthold;
    int pre_end2=pre_end1<<each_length;
    int pre_end3=pre_end2<<each_length;


    Word err_arry;


    /**运算**/
    int t_length_1=t_length-1;
    while(i<t_length_1)
    {


        X=Peq[text[i]]|VN;

        TMP=X&VP;

        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;

        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;

        D0=(TMP^VP)|X;

        HN=VP&D0;

        HP=VN|~(VP|D0);

        X=(D0>>1)&Hppp_mask;

        VN=X&HP;

        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;

        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3))
            return -1;

        Peq['A']=(Peq['A']>>1)&Hpp_Mask;
        Peq['T']=(Peq['T']>>1)&Hpp_Mask;
        Peq['G']=(Peq['G']>>1)&Hpp_Mask;
        Peq['C']=(Peq['C']>>1)&Hpp_Mask;





        ++i;
        ++i_bd;



        Peq[pattern1[i_bd]]=Peq[pattern1[i_bd]]|Mask1;
        Peq[pattern2[i_bd]]=Peq[pattern2[i_bd]]|Mask2;
        Peq[pattern3[i_bd]]=Peq[pattern3[i_bd]]|Mask3;


    }

        X=Peq[text[i]]|VN;

        TMP=X&VP;

        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;

        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;

        D0=(TMP^VP)|X;

        HN=VP&D0;

        HP=VN|~(VP|D0);

        X=(D0>>1)&Hppp_mask;

        VN=X&HP;

        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;

        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3))
            return -1;


    int site=p_length-last_high-1;

    err1=Err&single_err_mask1;
    err2=(Err&single_err_mask2)>>each_length;
    err3=(Err&single_err_mask3)>>(each_length+each_length);


     if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }



    i=0;
    Word tmp_VP;


    while(i<last_high)
    {

        tmp_VP=VP>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err+err_arry;



        tmp_VP=VN>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err-err_arry;


        ++i;


        err1=Err&single_err_mask1;
        err2=(Err&single_err_mask2)>>each_length;
        err3=(Err&single_err_mask3)>>(each_length+each_length);


         if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }

    }

    return 1;


}





/****************************************2条，固定区间长度的算法***************************************************************/
int Brief_2_Banded_BPM_Non_SSE(char *pattern1,char *pattern2,int p_length,char *text,int t_length,
                              int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    ///这两个纯粹是为了计算peq用的
    Word tmp_Peq_1=(Word)1;
    ///Word tmp_Peq_2=tmp_Peq_1<<(band_length+1);
    Word tmp_Peq_2=tmp_Peq_1<<(each_length);    ///注意each_length对应的是band_length+1而不是band_length

    int i=0;


    for (r =0; r<band_length; r++)
    {
        Peq[pattern1[r]]=Peq[pattern1[r]]|tmp_Peq_1;
        Peq[pattern2[r]]=Peq[pattern2[r]]|tmp_Peq_2;


        tmp_Peq_1=tmp_Peq_1<<1;
        tmp_Peq_2=tmp_Peq_2<<1;
    }


    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;

    Word All_1=(Word)-1;
    Word single_err_mask1=All_1>>(64-each_length);
    Word single_err_mask2=single_err_mask1<<(each_length);

    if(each_length<=16)
    {
        tmp_Peq[thread_id]['A'][0]=Peq_A&single_err_mask1;
        tmp_Peq[thread_id]['A'][1]=Peq_A>>each_length;

        tmp_Peq[thread_id]['T'][0]=Peq_T&single_err_mask1;
        tmp_Peq[thread_id]['T'][1]=Peq_T>>each_length;

        tmp_Peq[thread_id]['G'][0]=Peq_G&single_err_mask1;
        tmp_Peq[thread_id]['G'][1]=Peq_G>>each_length;

        tmp_Peq[thread_id]['C'][0]=Peq_C&single_err_mask1;
        tmp_Peq[thread_id]['C'][1]=Peq_C>>each_length;
    }
    else
    {

        tmp_Peq_32[thread_id]['A'][0]=Peq_A&single_err_mask1;
        tmp_Peq_32[thread_id]['A'][1]=Peq_A>>each_length;

        tmp_Peq_32[thread_id]['T'][0]=Peq_T&single_err_mask1;
        tmp_Peq_32[thread_id]['T'][1]=Peq_T>>each_length;

        tmp_Peq_32[thread_id]['G'][0]=Peq_G&single_err_mask1;
        tmp_Peq_32[thread_id]['G'][1]=Peq_G>>each_length;

        tmp_Peq_32[thread_id]['C'][0]=Peq_C&single_err_mask1;
        tmp_Peq_32[thread_id]['C'][1]=Peq_C>>each_length;


    }



    Word Mask1=(Word)1<<(band_length-1);
    Word Mask2=Mask1<<(each_length);    ///each_length对应的是band_length+1




    /*******************************这个是D0右移用的，好像也是加法也能用用的**************************************/
    ///这个是低band_length为1，高each_length-band_length为0
    Word  tmp_MASK=All_1>>(64-band_length);
    Word  Hppp_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    Word _Hppp_mask=~Hppp_mask;


    ///这个是循环右移的掩码
    tmp_MASK=All_1>>(64-band_length+1);
    Word Hpp_Mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;



    ///下面这两个变量是计算err用的
    ///这个处理方法必须要求需要计算编辑距离的区域的长度不能大于128
    ///计算err用的掩码就是低位全1
    tmp_MASK=(Word)1;
    Word err_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;



    Word VP=(Word)0;
    Word VN=(Word)0;
    Word X=(Word)0;
    Word D0=(Word)0;
    Word HN=(Word)0;
    Word HP=(Word)0;

    Word _TMP;

    Word TMP;


    Word Err=(Word)0;


    int s=0;

    i = 0;
    int j=0;
    ///int err[16];
    int err1;
    int err2;

    int i_bd=i+band_down;

    int last_high=band_length-t_length+p_length-band_down-1;



    int pre_end1=last_high+errthold;
    int pre_end2=pre_end1<<each_length;

    Word err_arry;


    /**运算**/
    int t_length_1=t_length-1;
    while(i<t_length_1)
    {



        X=Peq[text[i]]|VN;

        TMP=X&VP;


        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;

        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;


        D0=(TMP^VP)|X;



        HN=VP&D0;


        HP=VN|~(VP|D0);


        X=(D0>>1)&Hppp_mask;


        VN=X&HP;


        VP=HN|~(X|HP);


        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;

        if((err1>pre_end1)&&(err2>pre_end2))
            return -1;

        Peq['A']=(Peq['A']>>1)&Hpp_Mask;
        Peq['T']=(Peq['T']>>1)&Hpp_Mask;
        Peq['G']=(Peq['G']>>1)&Hpp_Mask;
        Peq['C']=(Peq['C']>>1)&Hpp_Mask;





        ++i;
        ++i_bd;



        Peq[pattern1[i_bd]]=Peq[pattern1[i_bd]]|Mask1;
        Peq[pattern2[i_bd]]=Peq[pattern2[i_bd]]|Mask2;


    }



        X=Peq[text[i]]|VN;

        TMP=X&VP;


        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;

        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;


        D0=(TMP^VP)|X;



        HN=VP&D0;


        HP=VN|~(VP|D0);


        X=(D0>>1)&Hppp_mask;


        VN=X&HP;


        VP=HN|~(X|HP);


        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;

        if((err1>pre_end1)&&(err2>pre_end2))
            return -1;


    int site=p_length-last_high-1;


    err1=Err&single_err_mask1;
    err2=(Err&single_err_mask2)>>each_length;


     if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }



    i=0;
    Word tmp_VP;

    while(i<last_high)
    {

        tmp_VP=VP>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err+err_arry;



        tmp_VP=VN>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err-err_arry;


        ++i;


        err1=Err&single_err_mask1;
        err2=(Err&single_err_mask2)>>each_length;


         if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }

    }

    return 1;


}





/****************************************4条，固定区间长度的算法, 前向***************************************************************/
int Start_location_Brief_4_Banded_BPM_Non_SSE(char *pattern1,char *pattern2, char* pattern3, char* pattern4,int p_length,char *text,int t_length,
                              int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    Word tmp_Peq_1=(Word)1;
    Word tmp_Peq_2=tmp_Peq_1<<(each_length);    ///注意each_length对应的是band_length+1而不是band_length
    Word tmp_Peq_3=tmp_Peq_2<<(each_length);
    Word tmp_Peq_4=tmp_Peq_3<<(each_length);


    int i=0;

    for (r =0; r<band_length; r++)
    {
        Peq[pattern1[p_length-r-1]]=Peq[pattern1[p_length-r-1]]|tmp_Peq_1;
        Peq[pattern2[p_length-r-1]]=Peq[pattern2[p_length-r-1]]|tmp_Peq_2;
        Peq[pattern3[p_length-r-1]]=Peq[pattern3[p_length-r-1]]|tmp_Peq_3;
        Peq[pattern4[p_length-r-1]]=Peq[pattern4[p_length-r-1]]|tmp_Peq_4;

        tmp_Peq_1=tmp_Peq_1<<1;
        tmp_Peq_2=tmp_Peq_2<<1;
        tmp_Peq_3=tmp_Peq_3<<1;
        tmp_Peq_4=tmp_Peq_4<<1;
    }


    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;

    Word All_1=(Word)-1;
    Word single_err_mask1=All_1>>(64-each_length);
    Word single_err_mask2=single_err_mask1<<(each_length);
    Word single_err_mask3=single_err_mask2<<(each_length);
    Word single_err_mask4=single_err_mask3<<(each_length);


    Word Mask1=(Word)1<<(band_length-1);
    Word Mask2=Mask1<<(each_length);    ///each_length对应的是band_length+1
    Word Mask3=Mask2<<(each_length);
    Word Mask4=Mask3<<(each_length);




    /*******************************这个是D0右移用的，好像也是加法也能用用的**************************************/
    ///这个是低band_length为1，高each_length-band_length为0
    Word  tmp_MASK=All_1>>(64-band_length);
    Word  Hppp_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;


    Word _Hppp_mask=~Hppp_mask;


    ///这个是循环右移的掩码
    tmp_MASK=All_1>>(64-band_length+1);
    Word Hpp_Mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;



    ///下面这两个变量是计算err用的
    ///这个处理方法必须要求需要计算编辑距离的区域的长度不能大于128
    ///计算err用的掩码就是低位全1
    tmp_MASK=(Word)1;
    Word err_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;



    Word VP=(Word)0;
    Word VN=(Word)0;
    Word X=(Word)0;
    Word D0=(Word)0;
    Word HN=(Word)0;
    Word HP=(Word)0;

    Word _TMP;

    Word TMP;

    Word Err=(Word)0;

    i=t_length-1;
    int j=0;
    int err1;
    int err2;
    int err3;
    int err4;

    int i_bd=p_length-band_length;
    int last_high=band_length-t_length+p_length-band_down-1;

    int pre_end1=last_high+errthold;
    int pre_end2=pre_end1<<each_length;
    int pre_end3=pre_end2<<each_length;
    int pre_end4=pre_end3<<each_length;
    Word err_arry;

    while(i>0)
    {
        X=Peq[text[i]]|VN;
        TMP=X&VP;
        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;
        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;
        D0=(TMP^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=(D0>>1)&Hppp_mask;
        VN=X&HP;
        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;
        err4=Err&single_err_mask4;

        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3)&&(err4>pre_end4))
            return -1;


        Peq['A']=(Peq['A']>>1)&Hpp_Mask;
        Peq['T']=(Peq['T']>>1)&Hpp_Mask;
        Peq['G']=(Peq['G']>>1)&Hpp_Mask;
        Peq['C']=(Peq['C']>>1)&Hpp_Mask;

        --i;
        --i_bd;

        Peq[pattern1[i_bd]]=Peq[pattern1[i_bd]]|Mask1;
        Peq[pattern2[i_bd]]=Peq[pattern2[i_bd]]|Mask2;
        Peq[pattern3[i_bd]]=Peq[pattern3[i_bd]]|Mask3;
        Peq[pattern4[i_bd]]=Peq[pattern4[i_bd]]|Mask4;


    }


        X=Peq[text[i]]|VN;
        TMP=X&VP;
        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;
        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;
        D0=(TMP^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=(D0>>1)&Hppp_mask;
        VN=X&HP;
        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;
        err4=Err&single_err_mask4;

        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3)&&(err4>pre_end4))
            return -1;


    int site=last_high;

    err1=Err&single_err_mask1;
    err2=(Err&single_err_mask2)>>each_length;
    err3=(Err&single_err_mask3)>>(each_length+each_length);
    err4=(Err&single_err_mask4)>>(each_length+each_length+each_length);




     if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }




    i=0;
    Word tmp_VP;


    while(i<last_high)
    {
        tmp_VP=VP>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err+err_arry;
        tmp_VP=VN>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err-err_arry;
        ++i;
        err1=Err&single_err_mask1;
        err2=(Err&single_err_mask2)>>each_length;
        err3=(Err&single_err_mask3)>>(each_length+each_length);
        err4=(Err&single_err_mask4)>>(each_length+each_length+each_length);

         if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site-i;
                return_sites_error[3]=err4;
        }

    }

    return 1;


}




/****************************************3条，固定区间长度的算法, 前向***************************************************************/
int Start_location_Brief_3_Banded_BPM_Non_SSE(char *pattern1,char *pattern2, char* pattern3, int p_length,char *text,int t_length,
                              int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    Word tmp_Peq_1=(Word)1;
    Word tmp_Peq_2=tmp_Peq_1<<(each_length);    ///注意each_length对应的是band_length+1而不是band_length
    Word tmp_Peq_3=tmp_Peq_2<<(each_length);

    int i=0;

    for (r =0; r<band_length; r++)
    {
        Peq[pattern1[p_length-r-1]]=Peq[pattern1[p_length-r-1]]|tmp_Peq_1;
        Peq[pattern2[p_length-r-1]]=Peq[pattern2[p_length-r-1]]|tmp_Peq_2;
        Peq[pattern3[p_length-r-1]]=Peq[pattern3[p_length-r-1]]|tmp_Peq_3;


        tmp_Peq_1=tmp_Peq_1<<1;
        tmp_Peq_2=tmp_Peq_2<<1;
        tmp_Peq_3=tmp_Peq_3<<1;
    }


    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;

    Word All_1=(Word)-1;
    Word single_err_mask1=All_1>>(64-each_length);
    Word single_err_mask2=single_err_mask1<<(each_length);
    Word single_err_mask3=single_err_mask2<<(each_length);



    Word Mask1=(Word)1<<(band_length-1);
    Word Mask2=Mask1<<(each_length);
    Word Mask3=Mask2<<(each_length);




    /*******************************这个是D0右移用的，好像也是加法也能用用的**************************************/
    ///这个是低band_length为1，高each_length-band_length为0
    Word  tmp_MASK=All_1>>(64-band_length);
    Word  Hppp_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;


    Word _Hppp_mask=~Hppp_mask;


    ///这个是循环右移的掩码
    tmp_MASK=All_1>>(64-band_length+1);
    Word Hpp_Mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;



    ///下面这两个变量是计算err用的
    ///这个处理方法必须要求需要计算编辑距离的区域的长度不能大于128
    ///计算err用的掩码就是低位全1
    tmp_MASK=(Word)1;
    Word err_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;



    Word VP=(Word)0;
    Word VN=(Word)0;
    Word X=(Word)0;
    Word D0=(Word)0;
    Word HN=(Word)0;
    Word HP=(Word)0;

    Word _TMP;

    Word TMP;
    Word Err=(Word)0;



    i=t_length-1;
    int j=0;
    int err1;
    int err2;
    int err3;
    int i_bd=p_length-band_length;
    int last_high=band_length-t_length+p_length-band_down-1;

    int pre_end1=last_high+errthold;
    int pre_end2=pre_end1<<each_length;
    int pre_end3=pre_end2<<each_length;
    Word err_arry;

    while(i>0)
    {

        X=Peq[text[i]]|VN;
        TMP=X&VP;
        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;
        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;
        D0=(TMP^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=(D0>>1)&Hppp_mask;
        VN=X&HP;
        VP=HN|~(X|HP);
        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;
        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;
        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3))
            return -1;



        Peq['A']=(Peq['A']>>1)&Hpp_Mask;
        Peq['T']=(Peq['T']>>1)&Hpp_Mask;
        Peq['G']=(Peq['G']>>1)&Hpp_Mask;
        Peq['C']=(Peq['C']>>1)&Hpp_Mask;

        --i;
        --i_bd;

        Peq[pattern1[i_bd]]=Peq[pattern1[i_bd]]|Mask1;
        Peq[pattern2[i_bd]]=Peq[pattern2[i_bd]]|Mask2;
        Peq[pattern3[i_bd]]=Peq[pattern3[i_bd]]|Mask3;


    }


        X=Peq[text[i]]|VN;
        TMP=X&VP;
        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;
        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;
        D0=(TMP^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=(D0>>1)&Hppp_mask;
        VN=X&HP;
        VP=HN|~(X|HP);
        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;
        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;
        err3=Err&single_err_mask3;
        if((err1>pre_end1)&&(err2>pre_end2)&&(err3>pre_end3))
            return -1;


    int site=last_high;
    err1=Err&single_err_mask1;
    err2=(Err&single_err_mask2)>>each_length;
    err3=(Err&single_err_mask3)>>(each_length+each_length);



     if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }




    i=0;
    Word tmp_VP;
    while(i<last_high)
    {
        tmp_VP=VP>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err+err_arry;
        tmp_VP=VN>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err-err_arry;

        ++i;
        err1=Err&single_err_mask1;
        err2=(Err&single_err_mask2)>>each_length;
        err3=(Err&single_err_mask3)>>(each_length+each_length);



         if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }

    }

    return 1;


}




/****************************************2条，固定区间长度的算法, 前向***************************************************************/
int Start_location_Brief_2_Banded_BPM_Non_SSE(char *pattern1,char *pattern2,int p_length,char *text,int t_length,
                              int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    Word tmp_Peq_1=(Word)1;
    Word tmp_Peq_2=tmp_Peq_1<<(each_length);    ///注意each_length对应的是band_length+1而不是band_length


    int i=0;


    for (r =0; r<band_length; r++)
    {
        Peq[pattern1[p_length-r-1]]=Peq[pattern1[p_length-r-1]]|tmp_Peq_1;
        Peq[pattern2[p_length-r-1]]=Peq[pattern2[p_length-r-1]]|tmp_Peq_2;


        tmp_Peq_1=tmp_Peq_1<<1;
        tmp_Peq_2=tmp_Peq_2<<1;
    }


    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;

    Word All_1=(Word)-1;
    Word single_err_mask1=All_1>>(64-each_length);
    Word single_err_mask2=single_err_mask1<<(each_length);

    Word Mask1=(Word)1<<(band_length-1);
    Word Mask2=Mask1<<(each_length);    ///each_length对应的是band_length+1




    /*******************************这个是D0右移用的，好像也是加法也能用用的**************************************/
    ///这个是低band_length为1，高each_length-band_length为0
    Word  tmp_MASK=All_1>>(64-band_length);
    Word  Hppp_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hppp_mask=Hppp_mask|tmp_MASK;

    Word _Hppp_mask=~Hppp_mask;


    ///这个是循环右移的掩码
    tmp_MASK=All_1>>(64-band_length+1);
    Word Hpp_Mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<(each_length);
    Hpp_Mask=Hpp_Mask|tmp_MASK;



    ///下面这两个变量是计算err用的
    ///这个处理方法必须要求需要计算编辑距离的区域的长度不能大于128
    ///计算err用的掩码就是低位全1
    tmp_MASK=(Word)1;
    Word err_mask=((Word)0)|tmp_MASK;

    tmp_MASK=tmp_MASK<<each_length;
    err_mask=err_mask|tmp_MASK;



    Word VP=(Word)0;
    Word VN=(Word)0;
    Word X=(Word)0;
    Word D0=(Word)0;
    Word HN=(Word)0;
    Word HP=(Word)0;

    Word _TMP;

    Word TMP;

    Word Err=(Word)0;

    i=t_length-1;
    int j=0;
    int err1;
    int err2;

    int i_bd=p_length-band_length;

    int last_high=band_length-t_length+p_length-band_down-1;

    int pre_end1=last_high+errthold;
    int pre_end2=pre_end1<<each_length;

    Word err_arry;


    while(i>0)
    {

        X=Peq[text[i]]|VN;
        TMP=X&VP;
        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;
        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;
        D0=(TMP^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=(D0>>1)&Hppp_mask;
        VN=X&HP;
        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;

        if((err1>pre_end1)&&(err2>pre_end2))
            return -1;

        Peq['A']=(Peq['A']>>1)&Hpp_Mask;
        Peq['T']=(Peq['T']>>1)&Hpp_Mask;
        Peq['G']=(Peq['G']>>1)&Hpp_Mask;
        Peq['C']=(Peq['C']>>1)&Hpp_Mask;


        --i;
        --i_bd;


        Peq[pattern1[i_bd]]=Peq[pattern1[i_bd]]|Mask1;
        Peq[pattern2[i_bd]]=Peq[pattern2[i_bd]]|Mask2;


    }


        X=Peq[text[i]]|VN;
        TMP=X&VP;
        _TMP=TMP^VP;
        _TMP=_TMP&_Hppp_mask;
        TMP=(TMP&Hppp_mask)+(VP&Hppp_mask);
        TMP=TMP^_TMP;
        D0=(TMP^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=(D0>>1)&Hppp_mask;
        VN=X&HP;
        VP=HN|~(X|HP);

        err_arry=D0&err_mask;
        Err=Err+err_mask;
        Err=Err-err_arry;

        err1=Err&single_err_mask1;
        err2=Err&single_err_mask2;

        if((err1>pre_end1)&&(err2>pre_end2))
            return -1;


    int site=last_high;

    err1=Err&single_err_mask1;
    err2=(Err&single_err_mask2)>>each_length;


     if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }



    i=0;
    Word tmp_VP;


    while(i<last_high)
    {

        tmp_VP=VP>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err+err_arry;

        tmp_VP=VN>>i;
        err_arry=tmp_VP&err_mask;
        Err=Err-err_arry;
        ++i;
        err1=Err&single_err_mask1;
        err2=(Err&single_err_mask2)>>each_length;
         if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }

    }


    return 1;


}




/***********************************************************************************************************************/
/***************************************************精简的8条start比对**************************************************/
int Start_location_Brief_Reserve_Banded_BPM_8_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,char *pattern8,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][7]=(unsigned short)0;




    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][7]=(unsigned short)0;


    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][7]=(unsigned short)0;



    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][7]=(unsigned short)0;





    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[p_length-r-1]][0]=tmp_Peq[thread_id][pattern1[p_length-r-1]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[p_length-r-1]][1]=tmp_Peq[thread_id][pattern2[p_length-r-1]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[p_length-r-1]][2]=tmp_Peq[thread_id][pattern3[p_length-r-1]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[p_length-r-1]][3]=tmp_Peq[thread_id][pattern4[p_length-r-1]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[p_length-r-1]][4]=tmp_Peq[thread_id][pattern5[p_length-r-1]][4]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern6[p_length-r-1]][5]=tmp_Peq[thread_id][pattern6[p_length-r-1]][5]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern7[p_length-r-1]][6]=tmp_Peq[thread_id][pattern7[p_length-r-1]][6]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern8[p_length-r-1]][7]=tmp_Peq[thread_id][pattern8[p_length-r-1]][7]|tmp_Peq_1;

        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 (tmp_Peq[thread_id]['A'][7], tmp_Peq[thread_id]['A'][6], tmp_Peq[thread_id]['A'][5], tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 (tmp_Peq[thread_id]['T'][7], tmp_Peq[thread_id]['T'][6], tmp_Peq[thread_id]['T'][5], tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 (tmp_Peq[thread_id]['G'][7], tmp_Peq[thread_id]['G'][6], tmp_Peq[thread_id]['G'][5], tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 (tmp_Peq[thread_id]['C'][7], tmp_Peq[thread_id]['C'][6], tmp_Peq[thread_id]['C'][5], tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);

    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask6=_mm_set_epi16((unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask7=_mm_set_epi16((unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask8=_mm_set_epi16(mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    int err6;
    int err7;
    int err8;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));


    i=t_length-1;
    int i_bd=p_length-band_length;



    while(i>0)
    {

        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);


        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
        cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
        if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6)&&_mm_extract_epi16(cmp_result,7))
            return 1;


        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);

        --i;
        --i_bd;

        Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
        Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
        Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
        Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
        Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);
        Peq_SSE[pattern6[i_bd]]=_mm_or_si128(Mask6,Peq_SSE[pattern6[i_bd]]);
        Peq_SSE[pattern7[i_bd]]=_mm_or_si128(Mask7,Peq_SSE[pattern7[i_bd]]);
        Peq_SSE[pattern8[i_bd]]=_mm_or_si128(Mask8,Peq_SSE[pattern8[i_bd]]);

    }



        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);


        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
        cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
        if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6)&&_mm_extract_epi16(cmp_result,7))
            return 1;



    int site=last_high;

    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);
    err6=_mm_extract_epi16(Err_8,5);
    err7=_mm_extract_epi16(Err_8,6);
    err8=_mm_extract_epi16(Err_8,7);
    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }
    if((err6<=errthold)&&(err6<return_sites_error[5]))
    {
            return_sites[5]=site;
            return_sites_error[5]=err6;
    }
    if((err7<=errthold)&&(err7<return_sites_error[6]))
    {
            return_sites[6]=site;
            return_sites_error[6]=err7;
    }
    if((err8<=errthold)&&(err8<return_sites_error[7]))
    {
            return_sites[7]=site;
            return_sites_error[7]=err8;
    }


    i=0;
    int pre_end_error=last_high+errthold;
    while(i<last_high)
    {
        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);
        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;
        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);
        err3=_mm_extract_epi16(Err_8,2);
        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);
        err6=_mm_extract_epi16(Err_8,5);
        err7=_mm_extract_epi16(Err_8,6);
        err8=_mm_extract_epi16(Err_8,7);

        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site-i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site-i;
                return_sites_error[4]=err5;
        }
        if((err6<=errthold)&&(err6<return_sites_error[5]))
        {
                return_sites[5]=site-i;
                return_sites_error[5]=err6;
        }
        if((err7<=errthold)&&(err7<return_sites_error[6]))
        {
                return_sites[6]=site-i;
                return_sites_error[6]=err7;
        }
        if((err8<=errthold)&&(err8<return_sites_error[7]))
        {
                return_sites[7]=site-i;
                return_sites_error[7]=err8;
        }

    }

    return 1;


}

/***************************************************************************************/


/***********************************************************************************************************************/
/***************************************************精简的7条start比对**************************************************/
int Start_location_Brief_Reserve_Banded_BPM_7_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][6]=(unsigned short)0;


    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][6]=(unsigned short)0;

    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][6]=(unsigned short)0;

    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][6]=(unsigned short)0;



    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[p_length-r-1]][0]=tmp_Peq[thread_id][pattern1[p_length-r-1]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[p_length-r-1]][1]=tmp_Peq[thread_id][pattern2[p_length-r-1]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[p_length-r-1]][2]=tmp_Peq[thread_id][pattern3[p_length-r-1]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[p_length-r-1]][3]=tmp_Peq[thread_id][pattern4[p_length-r-1]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[p_length-r-1]][4]=tmp_Peq[thread_id][pattern5[p_length-r-1]][4]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern6[p_length-r-1]][5]=tmp_Peq[thread_id][pattern6[p_length-r-1]][5]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern7[p_length-r-1]][6]=tmp_Peq[thread_id][pattern7[p_length-r-1]][6]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 ((unsigned int)0, tmp_Peq[thread_id]['A'][6], tmp_Peq[thread_id]['A'][5], tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 ((unsigned int)0, tmp_Peq[thread_id]['T'][6], tmp_Peq[thread_id]['T'][5], tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 ((unsigned int)0, tmp_Peq[thread_id]['G'][6], tmp_Peq[thread_id]['G'][5], tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 ((unsigned int)0, tmp_Peq[thread_id]['C'][6], tmp_Peq[thread_id]['C'][5], tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);


    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask6=_mm_set_epi16((unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask7=_mm_set_epi16((unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)0,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    int err6;
    int err7;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));


    i=t_length-1;
    int i_bd=p_length-band_length;



    while(i>0)
    {

        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6))
                return 1;




        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);

        --i;
        --i_bd;

        Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
        Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
        Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
        Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
        Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);
        Peq_SSE[pattern6[i_bd]]=_mm_or_si128(Mask6,Peq_SSE[pattern6[i_bd]]);
        Peq_SSE[pattern7[i_bd]]=_mm_or_si128(Mask7,Peq_SSE[pattern7[i_bd]]);

    }



        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6))
                return 1;



    int site=last_high;

    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);
    err6=_mm_extract_epi16(Err_8,5);
    err7=_mm_extract_epi16(Err_8,6);
    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }
    if((err6<=errthold)&&(err6<return_sites_error[5]))
    {
            return_sites[5]=site;
            return_sites_error[5]=err6;
    }
    if((err7<=errthold)&&(err7<return_sites_error[6]))
    {
            return_sites[6]=site;
            return_sites_error[6]=err7;
    }



    i=0;
    int pre_end_error=last_high+errthold;
    while(i<last_high)
    {
        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;
        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);

        err3=_mm_extract_epi16(Err_8,2);

        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);
        err6=_mm_extract_epi16(Err_8,5);
        err7=_mm_extract_epi16(Err_8,6);

        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site-i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site-i;
                return_sites_error[4]=err5;
        }
        if((err6<=errthold)&&(err6<return_sites_error[5]))
        {
                return_sites[5]=site-i;
                return_sites_error[5]=err6;
        }
        if((err7<=errthold)&&(err7<return_sites_error[6]))
        {
                return_sites[6]=site-i;
                return_sites_error[6]=err7;
        }

    }

    return 1;


}

/***************************************************************************************/


/***********************************************************************************************************************/
/***************************************************精简的6条start比对**************************************************/
int Start_location_Brief_Reserve_Banded_BPM_6_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }

    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][5]=(unsigned short)0;

    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][5]=(unsigned short)0;

    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][5]=(unsigned short)0;

    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][5]=(unsigned short)0;


    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[p_length-r-1]][0]=tmp_Peq[thread_id][pattern1[p_length-r-1]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[p_length-r-1]][1]=tmp_Peq[thread_id][pattern2[p_length-r-1]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[p_length-r-1]][2]=tmp_Peq[thread_id][pattern3[p_length-r-1]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[p_length-r-1]][3]=tmp_Peq[thread_id][pattern4[p_length-r-1]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[p_length-r-1]][4]=tmp_Peq[thread_id][pattern5[p_length-r-1]][4]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern6[p_length-r-1]][5]=tmp_Peq[thread_id][pattern6[p_length-r-1]][5]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['A'][5], tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['T'][5], tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['G'][5], tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['C'][5], tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);



    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask6=_mm_set_epi16((unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;
    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    int err6;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));

    i=t_length-1;
    int i_bd=p_length-band_length;


    while(i>0)
    {

        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5))
                return 1;


        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);


        --i;
        --i_bd;


        Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
        Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
        Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
        Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
        Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);
        Peq_SSE[pattern6[i_bd]]=_mm_or_si128(Mask6,Peq_SSE[pattern6[i_bd]]);

    }


        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5))
                return 1;


    int site=last_high;


    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);
    err6=_mm_extract_epi16(Err_8,5);
    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }
    if((err6<=errthold)&&(err6<return_sites_error[5]))
    {
            return_sites[5]=site;
            return_sites_error[5]=err6;
    }



    i=0;
    int pre_end_error=last_high+errthold;
    while(i<last_high)
    {
        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);
        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;
        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);
        err3=_mm_extract_epi16(Err_8,2);
        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);
        err6=_mm_extract_epi16(Err_8,5);
        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site-i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site-i;
                return_sites_error[4]=err5;
        }
        if((err6<=errthold)&&(err6<return_sites_error[5]))
        {
                return_sites[5]=site-i;
                return_sites_error[5]=err6;
        }


    }

    return 1;


}

/***************************************************************************************/


/***********************************************************************************************************************/
/***************************************************精简的5条start比对**************************************************/
int Start_location_Brief_Reserve_Banded_BPM_5_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;


    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;


    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;

    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;


    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[p_length-r-1]][0]=tmp_Peq[thread_id][pattern1[p_length-r-1]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[p_length-r-1]][1]=tmp_Peq[thread_id][pattern2[p_length-r-1]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[p_length-r-1]][2]=tmp_Peq[thread_id][pattern3[p_length-r-1]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[p_length-r-1]][3]=tmp_Peq[thread_id][pattern4[p_length-r-1]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[p_length-r-1]][4]=tmp_Peq[thread_id][pattern5[p_length-r-1]][4]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 ((unsigned int)0, (unsigned int)0, (unsigned int)0, tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);


    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));
    i=t_length-1;
    int i_bd=p_length-band_length;

    while(i>0)
    {

        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);

        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4))
                return 1;

        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);


        --i;
        --i_bd;



        Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
        Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
        Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
        Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
        Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);

    }


        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);

        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4))
                return 1;


    int site=last_high;

    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);

    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }




    i=0;
    int pre_end_error=last_high+errthold;
    while(i<last_high)
    {

        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;
        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);

        err3=_mm_extract_epi16(Err_8,2);
                //fprintf(chhy_w_fp1,"err[2]=%llu\n",_mm_extract_epi16(Err_8,2));

        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);

        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site-i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site-i;
                return_sites_error[4]=err5;
        }



    }

    return 1;


}

/***************************************************************************************/


/***********************************************************************************************************************/
/***************************************************精简的4条start比对(高阈值)**************************************************/
int Start_location_Brief_Reserve_Banded_BPM_4_high_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned int tmp_SSE;
    __m128i Peq_SSE[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned int tmp_Peq_1=(unsigned int)1;

    tmp_Peq_32[thread_id]['A'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][3]=(unsigned int)0;


    tmp_Peq_32[thread_id]['C'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][3]=(unsigned int)0;


    tmp_Peq_32[thread_id]['G'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][3]=(unsigned int)0;

    tmp_Peq_32[thread_id]['T'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][3]=(unsigned int)0;


    tmp_Peq_1=(unsigned int)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq_32[thread_id][pattern1[p_length-r-1]][0]=tmp_Peq_32[thread_id][pattern1[p_length-r-1]][0]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern2[p_length-r-1]][1]=tmp_Peq_32[thread_id][pattern2[p_length-r-1]][1]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern3[p_length-r-1]][2]=tmp_Peq_32[thread_id][pattern3[p_length-r-1]][2]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern4[p_length-r-1]][3]=tmp_Peq_32[thread_id][pattern4[p_length-r-1]][3]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi32 (tmp_Peq_32[thread_id]['A'][3], tmp_Peq_32[thread_id]['A'][2], tmp_Peq_32[thread_id]['A'][1], tmp_Peq_32[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi32 (tmp_Peq_32[thread_id]['T'][3], tmp_Peq_32[thread_id]['T'][2], tmp_Peq_32[thread_id]['T'][1], tmp_Peq_32[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi32 (tmp_Peq_32[thread_id]['G'][3], tmp_Peq_32[thread_id]['G'][2], tmp_Peq_32[thread_id]['G'][1], tmp_Peq_32[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi32 (tmp_Peq_32[thread_id]['C'][3], tmp_Peq_32[thread_id]['C'][2], tmp_Peq_32[thread_id]['C'][1], tmp_Peq_32[thread_id]['C'][0]);

    unsigned int mask_short=(unsigned int)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi32((unsigned int)0,(unsigned int)0,(unsigned int)0,mask_short);
    __m128i Mask2=_mm_set_epi32((unsigned int)0,(unsigned int)0,mask_short,(unsigned int)0);
    __m128i Mask3=_mm_set_epi32((unsigned int)0,mask_short,(unsigned int)0,(unsigned int)0);
    __m128i Mask4=_mm_set_epi32(mask_short,(unsigned int)0,(unsigned int)0,(unsigned int)0);
    __m128i HP_Mask=_mm_set_epi32((unsigned int)1,(unsigned int)1,(unsigned int)1,(unsigned int)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned int a_mask=((unsigned int)-1)>>(32-band_length);
    __m128i add_mask=_mm_set_epi32(a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi32((int)(last_high+errthold),(int)(last_high+errthold),(int)(last_high+errthold),(int)(last_high+errthold));

    i=t_length-1;
    int i_bd=p_length-band_length;

    while(i>0)
    {

        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);

        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
        cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
        if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,6))
                return 1;



        Peq_SSE['A']=_mm_srli_epi32(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi32(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi32(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi32(Peq_SSE['C'],1);


        --i;
        --i_bd;




        Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
        Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
        Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
        Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);

    }


        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);

        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
        cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
        if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,6))
                return 1;


    int site=last_high;

    err1=_mm_extract_epi16(Err_8,1);
    err1=err1<<16;
    err1=err1|_mm_extract_epi16(Err_8,0);

    err2=_mm_extract_epi16(Err_8,3);
    err2=err2<<16;
    err2=err2|_mm_extract_epi16(Err_8,2);

    err3=_mm_extract_epi16(Err_8,5);
    err3=err3<<16;
    err3=err3|_mm_extract_epi16(Err_8,4);

    err4=_mm_extract_epi16(Err_8,7);
    err4=err4<<16;
    err4=err4|_mm_extract_epi16(Err_8,6);

    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }





    i=0;
    int pre_end_error=last_high+errthold;

    while(i<last_high)
    {

        tmp_process=_mm_srli_epi32(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_add_epi32(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi32(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,tmp_process1);
        ++i;


        err1=_mm_extract_epi16(Err_8,1);
        err1=err1<<16;
        err1=err1|_mm_extract_epi16(Err_8,0);

        err2=_mm_extract_epi16(Err_8,3);
        err2=err2<<16;
        err2=err2|_mm_extract_epi16(Err_8,2);

        err3=_mm_extract_epi16(Err_8,5);
        err3=err3<<16;
        err3=err3|_mm_extract_epi16(Err_8,4);

        err4=_mm_extract_epi16(Err_8,7);
        err4=err4<<16;
        err4=err4|_mm_extract_epi16(Err_8,6);

        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site-i;
                return_sites_error[3]=err4;
        }



    }

    return 1;


}

/***************************************************************************************/






/***********************************************************************************************************************/
/***************************************************精简的3条start比对(高阈值)**************************************************/
int Start_location_Brief_Reserve_Banded_BPM_3_high_SSE(char *pattern1,char *pattern2,char *pattern3,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned int tmp_SSE;
    __m128i Peq_SSE[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned int tmp_Peq_1=(unsigned int)1;

    tmp_Peq_32[thread_id]['A'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][2]=(unsigned int)0;


    tmp_Peq_32[thread_id]['C'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][2]=(unsigned int)0;

    tmp_Peq_32[thread_id]['G'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][2]=(unsigned int)0;

    tmp_Peq_32[thread_id]['T'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][2]=(unsigned int)0;


    tmp_Peq_1=(unsigned int)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq_32[thread_id][pattern1[p_length-r-1]][0]=tmp_Peq_32[thread_id][pattern1[p_length-r-1]][0]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern2[p_length-r-1]][1]=tmp_Peq_32[thread_id][pattern2[p_length-r-1]][1]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern3[p_length-r-1]][2]=tmp_Peq_32[thread_id][pattern3[p_length-r-1]][2]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['A'][2], tmp_Peq_32[thread_id]['A'][1], tmp_Peq_32[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['T'][2], tmp_Peq_32[thread_id]['T'][1], tmp_Peq_32[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['G'][2], tmp_Peq_32[thread_id]['G'][1], tmp_Peq_32[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['C'][2], tmp_Peq_32[thread_id]['C'][1], tmp_Peq_32[thread_id]['C'][0]);


    unsigned int mask_short=(unsigned int)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi32((unsigned int)0,(unsigned int)0,(unsigned int)0,mask_short);
    __m128i Mask2=_mm_set_epi32((unsigned int)0,(unsigned int)0,mask_short,(unsigned int)0);
    __m128i Mask3=_mm_set_epi32((unsigned int)0,mask_short,(unsigned int)0,(unsigned int)0);
    __m128i HP_Mask=_mm_set_epi32((unsigned int)0,(unsigned int)1,(unsigned int)1,(unsigned int)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned int a_mask=((unsigned int)-1)>>(32-band_length);
    __m128i add_mask=_mm_set_epi32((unsigned int)0,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi32((unsigned int)0,(int)(last_high+errthold),(int)(last_high+errthold),(int)(last_high+errthold));


    i=t_length-1;
    int i_bd=p_length-band_length;



    while(i>0)
    {

        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);


        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
        cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
        if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4))
                return 1;


        Peq_SSE['A']=_mm_srli_epi32(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi32(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi32(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi32(Peq_SSE['C'],1);

        --i;
        --i_bd;



        Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
        Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
        Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);

    }


        X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);


        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
        cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
        if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4))
                return 1;


    int site=last_high;

    err1=_mm_extract_epi16(Err_8,1);
    err1=err1<<16;
    err1=err1|_mm_extract_epi16(Err_8,0);

    err2=_mm_extract_epi16(Err_8,3);
    err2=err2<<16;
    err2=err2|_mm_extract_epi16(Err_8,2);

    err3=_mm_extract_epi16(Err_8,5);
    err3=err3<<16;
    err3=err3|_mm_extract_epi16(Err_8,4);


    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }






    i=0;
    int pre_end_error=last_high+errthold;
    while(i<last_high)
    {
        tmp_process=_mm_srli_epi32(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_add_epi32(Err_8,tmp_process);

        tmp_process1=_mm_srli_epi32(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,tmp_process1);
        ++i;

        err1=_mm_extract_epi16(Err_8,1);
        err1=err1<<16;
        err1=err1|_mm_extract_epi16(Err_8,0);

        err2=_mm_extract_epi16(Err_8,3);
        err2=err2<<16;
        err2=err2|_mm_extract_epi16(Err_8,2);

        err3=_mm_extract_epi16(Err_8,5);
        err3=err3<<16;
        err3=err3|_mm_extract_epi16(Err_8,4);


        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site-i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site-i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site-i;
                return_sites_error[2]=err3;
        }




    }

    return 1;


}

/***************************************************************************************/


/**************************************************单条start比对***************************************************************************************************/
int Start_location_Reserve_Banded_BPM(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err,int thread_id)
{
    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Word tmp_Peq_1=(Word)1;


    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    for (r =0; r<band_length; r++)
    {
        Peq[pattern[p_length-r-1]]=Peq[pattern[p_length-r-1]]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;
    tmp_Peq[thread_id]['A'][0]=(unsigned short)Peq['A'];
    tmp_Peq[thread_id]['C'][0]=(unsigned short)Peq['C'];
    tmp_Peq[thread_id]['G'][0]=(unsigned short)Peq['G'];
    tmp_Peq[thread_id]['T'][0]=(unsigned short)Peq['T'];


    Word Mask_Pre=(Word)1<<(band_length-2);
    Word Mask=(Word)1<<(band_length-1);
    Word VP=0;
    Word VN=0;
    Word X=0;
    Word D0=0;
    Word HN=0;
    Word HP=0;
    int s=0;
    int i = 0;
    int j=0;

    int bound=band_length-2-band_down;
    int err=0;

    Word err_mask=(Word)1;
    int s1=band_length-2;
    i=t_length-1;
    int i_bd=p_length-band_length;
    int last_high=band_length-t_length+p_length-band_down-1;
    while(i>0)
    {
        X=Peq[text[i]]|VN;
        D0=((VP+(X&VP))^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);

        if(!(D0&err_mask))
         {
            ++err;
            if((err-last_high)>errthold)
                return -1;
        }

        /**
        for (symbol = 0; symbol < 4; symbol++)
        {
            Peq[Peq_index[symbol]]=Peq[Peq_index[symbol]]>>1;
        }
        **/

        Peq['A']=Peq['A']>>1;
        Peq['T']=Peq['T']>>1;
        Peq['G']=Peq['G']>>1;
        Peq['C']=Peq['C']>>1;

        --i;
        --i_bd;
        Peq[pattern[i_bd]]=Peq[pattern[i_bd]]|Mask;
    }


        X=Peq[text[i]]|VN;
        D0=((VP+(X&VP))^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);

        if(!(D0&err_mask))
         {
            ++err;
            if((err-last_high)>errthold)
                return -1;
        }


    int site=last_high;
    int return_site=-1;
    if((err<=errthold)&&(err<*return_err))
    {
            *return_err=err;
            return_site=site;
    }
    int i_last=i;
    i=0;
    while(i<last_high)
    {
        err=err+((VP>>i)&(Word)1);
        err=err-((VN>>i)&(Word)1);
        ++i;

        if((err<=errthold)&&(err<*return_err))
        {
                *return_err=err;
                return_site=site-i;
        }
    }

    return return_site;

}

/*************************************************************************************************************************************************************************/







/***********************************************************************************************************************/
/***************************************************精简的4条比对(适合于高阈值)**************************************************/
int Brief_Reserve_Banded_BPM_4_high_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned int tmp_SSE;
    __m128i Peq_SSE[128];
    //这里只定义为4,即只承认有ATGC是有道理的;
    //若pattern和text里同一个位置都有N，那么他还是不相等，所以N的列必须还为全0
    char Peq_index[4]= {'A','C','G','T'};
    //int p_length=strlen(pattern1);
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned int tmp_Peq_1=(unsigned int)1;

    tmp_Peq_32[thread_id]['A'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][3]=(unsigned int)0;




    tmp_Peq_32[thread_id]['C'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][3]=(unsigned int)0;


    tmp_Peq_32[thread_id]['G'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][3]=(unsigned int)0;



    tmp_Peq_32[thread_id]['T'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][2]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][3]=(unsigned int)0;





    tmp_Peq_1=(unsigned int)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq_32[thread_id][pattern1[r]][0]=tmp_Peq_32[thread_id][pattern1[r]][0]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern2[r]][1]=tmp_Peq_32[thread_id][pattern2[r]][1]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern3[r]][2]=tmp_Peq_32[thread_id][pattern3[r]][2]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern4[r]][3]=tmp_Peq_32[thread_id][pattern4[r]][3]|tmp_Peq_1;

        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi32 (tmp_Peq_32[thread_id]['A'][3], tmp_Peq_32[thread_id]['A'][2], tmp_Peq_32[thread_id]['A'][1], tmp_Peq_32[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi32 (tmp_Peq_32[thread_id]['T'][3], tmp_Peq_32[thread_id]['T'][2], tmp_Peq_32[thread_id]['T'][1], tmp_Peq_32[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi32 (tmp_Peq_32[thread_id]['G'][3], tmp_Peq_32[thread_id]['G'][2], tmp_Peq_32[thread_id]['G'][1], tmp_Peq_32[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi32 (tmp_Peq_32[thread_id]['C'][3], tmp_Peq_32[thread_id]['C'][2], tmp_Peq_32[thread_id]['C'][1], tmp_Peq_32[thread_id]['C'][0]);

    unsigned int mask_short=(unsigned int)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi32((unsigned int)0,(unsigned int)0,(unsigned int)0,mask_short);
    __m128i Mask2=_mm_set_epi32((unsigned int)0,(unsigned int)0,mask_short,(unsigned int)0);
    __m128i Mask3=_mm_set_epi32((unsigned int)0,mask_short,(unsigned int)0,(unsigned int)0);
    __m128i Mask4=_mm_set_epi32(mask_short,(unsigned int)0,(unsigned int)0,(unsigned int)0);
    __m128i HP_Mask=_mm_set_epi32((unsigned int)1,(unsigned int)1,(unsigned int)1,(unsigned int)1);


    //Word VP=0;
    __m128i VP=_mm_setzero_si128();
    //Word VN=0;
    __m128i VN=_mm_setzero_si128();
    //Word X=0;
    __m128i X=_mm_setzero_si128();
    //Word D0=0;
    __m128i D0=_mm_setzero_si128();
    //Word HN=0;
    __m128i HN=_mm_setzero_si128();
    //Word HP=0;
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned int a_mask=((unsigned int)-1)>>(32-band_length);
    __m128i add_mask=_mm_set_epi32(a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    //Word err_mask=(Word)1<<(band_length-1);
    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi32((int)(last_high+errthold),(int)(last_high+errthold),(int)(last_high+errthold),(int)(last_high+errthold));

    //int s1=band_length-2;
    int i_bd=i+band_down;

   int t_length_1=t_length-1;

    while(i<t_length_1)
    {

         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);


        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);

        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,6))
                return 1;



        Peq_SSE['A']=_mm_srli_epi32(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi32(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi32(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi32(Peq_SSE['C'],1);


        ++i;
        ++i_bd;


            Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
            Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
            Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
            Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);

    }


             X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);


        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);

        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,6))
                return 1;



    int site=p_length-last_high-1;

    err1=_mm_extract_epi16(Err_8,1);
    err1=err1<<16;
    err1=err1|_mm_extract_epi16(Err_8,0);

    err2=_mm_extract_epi16(Err_8,3);
    err2=err2<<16;
    err2=err2|_mm_extract_epi16(Err_8,2);

    err3=_mm_extract_epi16(Err_8,5);
    err3=err3<<16;
    err3=err3|_mm_extract_epi16(Err_8,4);

    err4=_mm_extract_epi16(Err_8,7);
    err4=err4<<16;
    err4=err4|_mm_extract_epi16(Err_8,6);



    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }






    i=0;


    while(i<last_high)
    {
        //HP_Mask;


        tmp_process=_mm_srli_epi32(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_add_epi32(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi32(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,tmp_process1);
        ++i;


        err1=_mm_extract_epi16(Err_8,1);
        err1=err1<<16;
        err1=err1|_mm_extract_epi16(Err_8,0);

        err2=_mm_extract_epi16(Err_8,3);
        err2=err2<<16;
        err2=err2|_mm_extract_epi16(Err_8,2);

        err3=_mm_extract_epi16(Err_8,5);
        err3=err3<<16;
        err3=err3|_mm_extract_epi16(Err_8,4);

        err4=_mm_extract_epi16(Err_8,7);
        err4=err4<<16;
        err4=err4|_mm_extract_epi16(Err_8,6);



        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site+i;
                return_sites_error[3]=err4;
        }




    }

    return 1;


}

/***************************************************************************************/

/***********************************************************************************************************************/
/***************************************************精简的3条比对(适合于高阈值)**************************************************/
int Brief_Reserve_Banded_BPM_3_high_SSE(char *pattern1,char *pattern2,char *pattern3,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned int tmp_SSE;
    __m128i Peq_SSE[128];
    //这里只定义为4,即只承认有ATGC是有道理的;
    //若pattern和text里同一个位置都有N，那么他还是不相等，所以N的列必须还为全0
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned int tmp_Peq_1=(unsigned int)1;

    tmp_Peq_32[thread_id]['A'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['A'][2]=(unsigned int)0;




    tmp_Peq_32[thread_id]['C'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['C'][2]=(unsigned int)0;


    tmp_Peq_32[thread_id]['G'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['G'][2]=(unsigned int)0;



    tmp_Peq_32[thread_id]['T'][0]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][1]=(unsigned int)0;
    tmp_Peq_32[thread_id]['T'][2]=(unsigned int)0;





    tmp_Peq_1=(unsigned int)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq_32[thread_id][pattern1[r]][0]=tmp_Peq_32[thread_id][pattern1[r]][0]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern2[r]][1]=tmp_Peq_32[thread_id][pattern2[r]][1]|tmp_Peq_1;
        tmp_Peq_32[thread_id][pattern3[r]][2]=tmp_Peq_32[thread_id][pattern3[r]][2]|tmp_Peq_1;

        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['A'][2], tmp_Peq_32[thread_id]['A'][1], tmp_Peq_32[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['T'][2], tmp_Peq_32[thread_id]['T'][1], tmp_Peq_32[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['G'][2], tmp_Peq_32[thread_id]['G'][1], tmp_Peq_32[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi32 ((unsigned int)0, tmp_Peq_32[thread_id]['C'][2], tmp_Peq_32[thread_id]['C'][1], tmp_Peq_32[thread_id]['C'][0]);


    unsigned int mask_short=(unsigned int)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi32((unsigned int)0,(unsigned int)0,(unsigned int)0,mask_short);
    __m128i Mask2=_mm_set_epi32((unsigned int)0,(unsigned int)0,mask_short,(unsigned int)0);
    __m128i Mask3=_mm_set_epi32((unsigned int)0,mask_short,(unsigned int)0,(unsigned int)0);
    __m128i HP_Mask=_mm_set_epi32((unsigned int)0,(unsigned int)1,(unsigned int)1,(unsigned int)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned int a_mask=((unsigned int)-1)>>(32-band_length);
    __m128i add_mask=_mm_set_epi32((unsigned int)0,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi32((unsigned int)0,(int)(last_high+errthold),(int)(last_high+errthold),(int)(last_high+errthold));

    int i_bd=i+band_down;



   int t_length_1=t_length-1;

    while(i<t_length_1)
    {


         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);

        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);


        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);


        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4))
                return 1;


        Peq_SSE['A']=_mm_srli_epi32(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi32(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi32(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi32(Peq_SSE['C'],1);


        ++i;
        ++i_bd;

                Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
                Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
                Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);

    }




         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);

        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_add_epi32(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);


        X=_mm_srli_epi32(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);


        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_add_epi32(Err_8,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi32(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,4))
                return 1;



    int site=p_length-last_high-1;


    err1=_mm_extract_epi16(Err_8,1);
    err1=err1<<16;
    err1=err1|_mm_extract_epi16(Err_8,0);

    err2=_mm_extract_epi16(Err_8,3);
    err2=err2<<16;
    err2=err2|_mm_extract_epi16(Err_8,2);

    err3=_mm_extract_epi16(Err_8,5);
    err3=err3<<16;
    err3=err3|_mm_extract_epi16(Err_8,4);





    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }







    i=0;


    while(i<last_high)
    {


        tmp_process=_mm_srli_epi32(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_add_epi32(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi32(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_sub_epi32(Err_8,tmp_process1);
        ++i;


        err1=_mm_extract_epi16(Err_8,1);
        err1=err1<<16;
        err1=err1|_mm_extract_epi16(Err_8,0);

        err2=_mm_extract_epi16(Err_8,3);
        err2=err2<<16;
        err2=err2|_mm_extract_epi16(Err_8,2);

        err3=_mm_extract_epi16(Err_8,5);
        err3=err3<<16;
        err3=err3|_mm_extract_epi16(Err_8,4);







        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }





    }

    return 1;


}

/***************************************************************************************/



/***********************************************************************************************************************/
/***************************************************精简的8条比对**************************************************/
int Brief_Reserve_Banded_BPM_8_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,char *pattern8,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    //这里只定义为4,即只承认有ATGC是有道理的;
    //若pattern和text里同一个位置都有N，那么他还是不相等，所以N的列必须还为全0
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][7]=(unsigned short)0;




    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][7]=(unsigned short)0;


    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][7]=(unsigned short)0;



    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][6]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][7]=(unsigned short)0;





    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[r]][0]=tmp_Peq[thread_id][pattern1[r]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[r]][1]=tmp_Peq[thread_id][pattern2[r]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[r]][2]=tmp_Peq[thread_id][pattern3[r]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[r]][3]=tmp_Peq[thread_id][pattern4[r]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[r]][4]=tmp_Peq[thread_id][pattern5[r]][4]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern6[r]][5]=tmp_Peq[thread_id][pattern6[r]][5]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern7[r]][6]=tmp_Peq[thread_id][pattern7[r]][6]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern8[r]][7]=tmp_Peq[thread_id][pattern8[r]][7]|tmp_Peq_1;

        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 (tmp_Peq[thread_id]['A'][7], tmp_Peq[thread_id]['A'][6], tmp_Peq[thread_id]['A'][5], tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 (tmp_Peq[thread_id]['T'][7], tmp_Peq[thread_id]['T'][6], tmp_Peq[thread_id]['T'][5], tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 (tmp_Peq[thread_id]['G'][7], tmp_Peq[thread_id]['G'][6], tmp_Peq[thread_id]['G'][5], tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 (tmp_Peq[thread_id]['C'][7], tmp_Peq[thread_id]['C'][6], tmp_Peq[thread_id]['C'][5], tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);


    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask6=_mm_set_epi16((unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask7=_mm_set_epi16((unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask8=_mm_set_epi16(mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    int err6;
    int err7;
    int err8;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    ///zui gao wei shi 0, fang zhi yi chu
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));
    int i_bd=i+band_down;


    int t_length_1=t_length-1;
    while(i<t_length_1)
    {

         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6)&&_mm_extract_epi16(cmp_result,7))
                return 1;

        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);


        ++i;
        ++i_bd;



            Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
            Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
            Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
            Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
            Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);
            Peq_SSE[pattern6[i_bd]]=_mm_or_si128(Mask6,Peq_SSE[pattern6[i_bd]]);
            Peq_SSE[pattern7[i_bd]]=_mm_or_si128(Mask7,Peq_SSE[pattern7[i_bd]]);
            Peq_SSE[pattern8[i_bd]]=_mm_or_si128(Mask8,Peq_SSE[pattern8[i_bd]]);

    }




         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);

        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6)&&_mm_extract_epi16(cmp_result,7))
                return 1;










    int site=p_length-last_high-1;

    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);
    err6=_mm_extract_epi16(Err_8,5);
    err7=_mm_extract_epi16(Err_8,6);
    err8=_mm_extract_epi16(Err_8,7);
    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }
    if((err6<=errthold)&&(err6<return_sites_error[5]))
    {
            return_sites[5]=site;
            return_sites_error[5]=err6;
    }
    if((err7<=errthold)&&(err7<return_sites_error[6]))
    {
            return_sites[6]=site;
            return_sites_error[6]=err7;
    }
    if((err8<=errthold)&&(err8<return_sites_error[7]))
    {
            return_sites[7]=site;
            return_sites_error[7]=err8;
    }


    i=0;
    int pre_end_error=last_high+errthold;

    while(i<last_high)
    {
        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;

        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);
        err3=_mm_extract_epi16(Err_8,2);
        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);
        err6=_mm_extract_epi16(Err_8,5);
        err7=_mm_extract_epi16(Err_8,6);
        err8=_mm_extract_epi16(Err_8,7);


        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site+i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site+i;
                return_sites_error[4]=err5;
        }
        if((err6<=errthold)&&(err6<return_sites_error[5]))
        {
                return_sites[5]=site+i;
                return_sites_error[5]=err6;
        }
        if((err7<=errthold)&&(err7<return_sites_error[6]))
        {
                return_sites[6]=site+i;
                return_sites_error[6]=err7;
        }
        if((err8<=errthold)&&(err8<return_sites_error[7]))
        {
                return_sites[7]=site+i;
                return_sites_error[7]=err8;
        }

    }

    return 1;


}


/**************************************************精简的7条比对****************************************************/
int Brief_Reserve_Banded_BPM_7_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,char *pattern7,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    //这里只定义为4,即只承认有ATGC是有道理的;
    //若pattern和text里同一个位置都有N，那么他还是不相等，所以N的列必须还为全0
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][6]=(unsigned short)0;




    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][6]=(unsigned short)0;


    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][6]=(unsigned short)0;



    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][5]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][6]=(unsigned short)0;





    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[r]][0]=tmp_Peq[thread_id][pattern1[r]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[r]][1]=tmp_Peq[thread_id][pattern2[r]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[r]][2]=tmp_Peq[thread_id][pattern3[r]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[r]][3]=tmp_Peq[thread_id][pattern4[r]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[r]][4]=tmp_Peq[thread_id][pattern5[r]][4]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern6[r]][5]=tmp_Peq[thread_id][pattern6[r]][5]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern7[r]][6]=tmp_Peq[thread_id][pattern7[r]][6]|tmp_Peq_1;

        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 ((unsigned short)0, tmp_Peq[thread_id]['A'][6], tmp_Peq[thread_id]['A'][5], tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 ((unsigned short)0, tmp_Peq[thread_id]['T'][6], tmp_Peq[thread_id]['T'][5], tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 ((unsigned short)0, tmp_Peq[thread_id]['G'][6], tmp_Peq[thread_id]['G'][5], tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 ((unsigned short)0, tmp_Peq[thread_id]['C'][6], tmp_Peq[thread_id]['C'][5], tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);


    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask6=_mm_set_epi16((unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask7=_mm_set_epi16((unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)0,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    int err6;
    int err7;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));
    int i_bd=i+band_down;


    int t_length_1=t_length-1;
    while(i<t_length_1)
    {

         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);


        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);

         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6))
                return 1;

        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);

        ++i;
        ++i_bd;

        Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
        Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
        Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
        Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
        Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);
        Peq_SSE[pattern6[i_bd]]=_mm_or_si128(Mask6,Peq_SSE[pattern6[i_bd]]);
        Peq_SSE[pattern7[i_bd]]=_mm_or_si128(Mask7,Peq_SSE[pattern7[i_bd]]);

    }



      X=_mm_or_si128 (Peq_SSE[text[i]],VN);
    tmp_process1=_mm_and_si128(X, VP);
    /*******************************errthold!=5启用下面的算法***********************************************************/
    /*****************************************下面这两步是纯粹防止溢出*************************/
    tmp_process1=_mm_and_si128(tmp_process1,add_mask);
    VP=_mm_and_si128(VP,add_mask);
    /*************************************************************************************/
    tmp_process=_mm_adds_epu16(tmp_process1, VP);
    tmp_process=_mm_xor_si128(tmp_process, VP);
    D0=_mm_or_si128 (tmp_process, X);
    HN=_mm_and_si128(D0, VP);
    tmp_process=_mm_or_si128(D0, VP);
    tmp_process=_mm_andnot_si128(tmp_process,for_not);
    HP=_mm_or_si128(tmp_process, VN);
    X=_mm_srli_epi16(D0,1);
    VN=_mm_and_si128(X,HP);
    tmp_process=_mm_or_si128(X,HP);
    tmp_process=_mm_andnot_si128(tmp_process,for_not);
    VP=_mm_or_si128(HN,tmp_process);


    err_arry=_mm_and_si128(D0,err_mask);
    Err_8=_mm_adds_epi16(Err_8,HP_Mask);
    Err_8=_mm_subs_epi16(Err_8,err_arry);

     cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
     if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5)&&_mm_extract_epi16(cmp_result,6))
            return 1;



    int site=p_length-last_high-1;
    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);
    err6=_mm_extract_epi16(Err_8,5);
    err7=_mm_extract_epi16(Err_8,6);
    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }
    if((err6<=errthold)&&(err6<return_sites_error[5]))
    {
            return_sites[5]=site;
            return_sites_error[5]=err6;
    }
    if((err7<=errthold)&&(err7<return_sites_error[6]))
    {
            return_sites[6]=site;
            return_sites_error[6]=err7;
    }



    i=0;


    while(i<last_high)
    {
        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;


        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);

        err3=_mm_extract_epi16(Err_8,2);

        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);
        err6=_mm_extract_epi16(Err_8,5);
        err7=_mm_extract_epi16(Err_8,6);


        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site+i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site+i;
                return_sites_error[4]=err5;
        }
        if((err6<=errthold)&&(err6<return_sites_error[5]))
        {
                return_sites[5]=site+i;
                return_sites_error[5]=err6;
        }
        if((err7<=errthold)&&(err7<return_sites_error[6]))
        {
                return_sites[6]=site+i;
                return_sites_error[6]=err7;
        }

    }

    return 1;


}


/***************************************************************************************/


/**************************************************精简的6条比对****************************************************/
int Brief_Reserve_Banded_BPM_6_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,char *pattern6,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    //这里只定义为4,即只承认有ATGC是有道理的;
    //若pattern和text里同一个位置都有N，那么他还是不相等，所以N的列必须还为全0
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][5]=(unsigned short)0;




    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][5]=(unsigned short)0;


    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][5]=(unsigned short)0;



    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][5]=(unsigned short)0;





    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[r]][0]=tmp_Peq[thread_id][pattern1[r]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[r]][1]=tmp_Peq[thread_id][pattern2[r]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[r]][2]=tmp_Peq[thread_id][pattern3[r]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[r]][3]=tmp_Peq[thread_id][pattern4[r]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[r]][4]=tmp_Peq[thread_id][pattern5[r]][4]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern6[r]][5]=tmp_Peq[thread_id][pattern6[r]][5]|tmp_Peq_1;

        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['A'][5], tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['T'][5], tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['G'][5], tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['C'][5], tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);


    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask6=_mm_set_epi16((unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    int err6;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));

    int i_bd=i+band_down;


    int t_length_1=t_length-1;

    while(i<t_length_1)
    {

         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5))
                return 1;


        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);


        ++i;
        ++i_bd;

            Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
            Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
            Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
            Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
            Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);
            Peq_SSE[pattern6[i_bd]]=_mm_or_si128(Mask6,Peq_SSE[pattern6[i_bd]]);

    }


         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4)&&_mm_extract_epi16(cmp_result,5))
                return 1;



    int site=p_length-last_high-1;
    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);
    err6=_mm_extract_epi16(Err_8,5);
    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }
    if((err6<=errthold)&&(err6<return_sites_error[5]))
    {
            return_sites[5]=site;
            return_sites_error[5]=err6;
    }


    i=0;
    while(i<last_high)
    {
        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;


        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);

        err3=_mm_extract_epi16(Err_8,2);

        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);
        err6=_mm_extract_epi16(Err_8,5);




        /*************************************这里需要看是先抽出来再比较好还是先一起比较再抽出来好*****************************************************************/

        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site+i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site+i;
                return_sites_error[4]=err5;
        }
        if((err6<=errthold)&&(err6<return_sites_error[5]))
        {
                return_sites[5]=site+i;
                return_sites_error[5]=err6;
        }


    }

    return 1;


}


/***************************************************************************************************************************************************/

/**************************************************精简的5条比对****************************************************/
int Brief_Reserve_Banded_BPM_5_SSE(char *pattern1,char *pattern2,char *pattern3,char *pattern4,char *pattern5,int p_length,char *text,int t_length,
                int* return_sites,int* return_sites_error,unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length, int thread_id)
{
    Word Peq[128];
    unsigned short tmp_SSE;
    __m128i Peq_SSE[128];
    //这里只定义为4,即只承认有ATGC是有道理的;
    //若pattern和text里同一个位置都有N，那么他还是不相等，所以N的列必须还为全0
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq_SSE[symbol]=_mm_setzero_si128();
    }


    unsigned short tmp_Peq_1=(unsigned short)1;

    tmp_Peq[thread_id]['A'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['A'][4]=(unsigned short)0;




    tmp_Peq[thread_id]['C'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['C'][4]=(unsigned short)0;


    tmp_Peq[thread_id]['G'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['G'][4]=(unsigned short)0;



    tmp_Peq[thread_id]['T'][0]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][1]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][2]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][3]=(unsigned short)0;
    tmp_Peq[thread_id]['T'][4]=(unsigned short)0;





    tmp_Peq_1=(unsigned short)1;
    for (r =0; r<band_length; r++)
    {
        tmp_Peq[thread_id][pattern1[r]][0]=tmp_Peq[thread_id][pattern1[r]][0]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern2[r]][1]=tmp_Peq[thread_id][pattern2[r]][1]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern3[r]][2]=tmp_Peq[thread_id][pattern3[r]][2]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern4[r]][3]=tmp_Peq[thread_id][pattern4[r]][3]|tmp_Peq_1;
        tmp_Peq[thread_id][pattern5[r]][4]=tmp_Peq[thread_id][pattern5[r]][4]|tmp_Peq_1;

        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_SSE['A']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['A'][4], tmp_Peq[thread_id]['A'][3], tmp_Peq[thread_id]['A'][2], tmp_Peq[thread_id]['A'][1], tmp_Peq[thread_id]['A'][0]);
    Peq_SSE['T']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['T'][4], tmp_Peq[thread_id]['T'][3], tmp_Peq[thread_id]['T'][2], tmp_Peq[thread_id]['T'][1], tmp_Peq[thread_id]['T'][0]);
    Peq_SSE['G']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['G'][4], tmp_Peq[thread_id]['G'][3], tmp_Peq[thread_id]['G'][2], tmp_Peq[thread_id]['G'][1], tmp_Peq[thread_id]['G'][0]);
    Peq_SSE['C']=_mm_set_epi16 ((unsigned short)0, (unsigned short)0, (unsigned short)0, tmp_Peq[thread_id]['C'][4], tmp_Peq[thread_id]['C'][3], tmp_Peq[thread_id]['C'][2], tmp_Peq[thread_id]['C'][1], tmp_Peq[thread_id]['C'][0]);

    unsigned short mask_short=(unsigned short)1<<(band_length-1);
    __m128i Mask1=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short);
    __m128i Mask2=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0);
    __m128i Mask3=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0);
    __m128i Mask4=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i Mask5=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,mask_short,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);
    __m128i HP_Mask=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1,(unsigned short)1);


    __m128i VP=_mm_setzero_si128();
    __m128i VN=_mm_setzero_si128();
    __m128i X=_mm_setzero_si128();
    __m128i D0=_mm_setzero_si128();
    __m128i HN=_mm_setzero_si128();
    __m128i HP=_mm_setzero_si128();
    __m128i tmp_process;
    __m128i tmp_process1;
    __m128i for_not=_mm_set1_epi32(-1);

    __m128i Err_8=_mm_set_epi16((unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0,(unsigned short)0);

    int s=0;

    int i = 0;
    int j=0;


    int err1;
    int err2;
    int err3;
    int err4;
    int err5;
    /******************************************************************************************/
    /**********************************这个是为errthold！=5专门新加的****************************/
    unsigned short a_mask=((unsigned short)-1)>>(16-band_length);
    __m128i add_mask=_mm_set_epi16(a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask,a_mask);
    /******************************************************************************************/

    __m128i err_mask=HP_Mask;
    __m128i err_arry;
    __m128i cmp_result;
    int last_high=band_length-t_length+p_length-band_down-1;
    __m128i pre_end=_mm_set_epi16((short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold),(short)(last_high+errthold));

    int i_bd=i+band_down;

    int t_length_1=t_length-1;

    while(i<t_length_1)
    {

         X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4))
                return 1;


        Peq_SSE['A']=_mm_srli_epi16(Peq_SSE['A'],1);
        Peq_SSE['T']=_mm_srli_epi16(Peq_SSE['T'],1);
        Peq_SSE['G']=_mm_srli_epi16(Peq_SSE['G'],1);
        Peq_SSE['C']=_mm_srli_epi16(Peq_SSE['C'],1);


        ++i;
        ++i_bd;

                Peq_SSE[pattern1[i_bd]]=_mm_or_si128(Mask1,Peq_SSE[pattern1[i_bd]]);
                Peq_SSE[pattern2[i_bd]]=_mm_or_si128(Mask2,Peq_SSE[pattern2[i_bd]]);
                Peq_SSE[pattern3[i_bd]]=_mm_or_si128(Mask3,Peq_SSE[pattern3[i_bd]]);
                Peq_SSE[pattern4[i_bd]]=_mm_or_si128(Mask4,Peq_SSE[pattern4[i_bd]]);
                Peq_SSE[pattern5[i_bd]]=_mm_or_si128(Mask5,Peq_SSE[pattern5[i_bd]]);

    }




      X=_mm_or_si128 (Peq_SSE[text[i]],VN);
        tmp_process1=_mm_and_si128(X, VP);
        /*******************************errthold!=5启用下面的算法***********************************************************/
        /*****************************************下面这两步是纯粹防止溢出*************************/
        tmp_process1=_mm_and_si128(tmp_process1,add_mask);
        VP=_mm_and_si128(VP,add_mask);
        /*************************************************************************************/
        tmp_process=_mm_adds_epu16(tmp_process1, VP);
        tmp_process=_mm_xor_si128(tmp_process, VP);
        D0=_mm_or_si128 (tmp_process, X);
        HN=_mm_and_si128(D0, VP);
        tmp_process=_mm_or_si128(D0, VP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        HP=_mm_or_si128(tmp_process, VN);
        X=_mm_srli_epi16(D0,1);
        VN=_mm_and_si128(X,HP);
        tmp_process=_mm_or_si128(X,HP);
        tmp_process=_mm_andnot_si128(tmp_process,for_not);
        VP=_mm_or_si128(HN,tmp_process);
        err_arry=_mm_and_si128(D0,err_mask);
        Err_8=_mm_adds_epi16(Err_8,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,err_arry);
         cmp_result=_mm_cmpgt_epi16(Err_8,pre_end);
         if(_mm_extract_epi16(cmp_result,0)&&_mm_extract_epi16(cmp_result,1)&&_mm_extract_epi16(cmp_result,2)&&_mm_extract_epi16(cmp_result,3)&&_mm_extract_epi16(cmp_result,4))
                return 1;


    int site=p_length-last_high-1;

    err1=_mm_extract_epi16(Err_8,0);
    err2=_mm_extract_epi16(Err_8,1);
    err3=_mm_extract_epi16(Err_8,2);
    err4=_mm_extract_epi16(Err_8,3);
    err5=_mm_extract_epi16(Err_8,4);
    if((err1<=errthold)&&(err1<return_sites_error[0]))
    {
            return_sites[0]=site;
            return_sites_error[0]=err1;
    }
    if((err2<=errthold)&&(err2<return_sites_error[1]))
    {
            return_sites[1]=site;
            return_sites_error[1]=err2;
    }
    if((err3<=errthold)&&(err3<return_sites_error[2]))
    {
            return_sites[2]=site;
            return_sites_error[2]=err3;
    }
    if((err4<=errthold)&&(err4<return_sites_error[3]))
    {
            return_sites[3]=site;
            return_sites_error[3]=err4;
    }
    if((err5<=errthold)&&(err5<return_sites_error[4]))
    {
            return_sites[4]=site;
            return_sites_error[4]=err5;
    }





    i=0;


    while(i<last_high)
    {

        tmp_process=_mm_srli_epi16(VP,i);
        tmp_process=_mm_and_si128(tmp_process,HP_Mask);
        Err_8=_mm_adds_epi16(Err_8,tmp_process);


        tmp_process1=_mm_srli_epi16(VN,i);
        tmp_process1=_mm_and_si128(tmp_process1,HP_Mask);
        Err_8=_mm_subs_epi16(Err_8,tmp_process1);
        ++i;


        err1=_mm_extract_epi16(Err_8,0);
        err2=_mm_extract_epi16(Err_8,1);

        err3=_mm_extract_epi16(Err_8,2);

        err4=_mm_extract_epi16(Err_8,3);
        err5=_mm_extract_epi16(Err_8,4);





        if((err1<=errthold)&&(err1<return_sites_error[0]))
        {
                return_sites[0]=site+i;
                return_sites_error[0]=err1;
        }
        if((err2<=errthold)&&(err2<return_sites_error[1]))
        {
                return_sites[1]=site+i;
                return_sites_error[1]=err2;
        }
        if((err3<=errthold)&&(err3<return_sites_error[2]))
        {
                return_sites[2]=site+i;
                return_sites_error[2]=err3;
        }
        if((err4<=errthold)&&(err4<return_sites_error[3]))
        {
                return_sites[3]=site+i;
                return_sites_error[3]=err4;
        }
        if((err5<=errthold)&&(err5<return_sites_error[4]))
        {
                return_sites[4]=site+i;
                return_sites_error[4]=err5;
        }



    }

    return 1;


}


/*************************************************************************************************************************************************/

/**************************************************单条比对***************************************************************************************************/
int Reserve_Banded_BPM(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err, int thread_id)
{

    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Word tmp_Peq_1=(Word)1;


    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    for (r =0; r<band_length; r++)
    {
            Peq[pattern[r]]=Peq[pattern[r]]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
        //jump_c[symbol]=128;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;
    tmp_Peq[thread_id]['A'][0]=(unsigned short)Peq['A'];
    tmp_Peq[thread_id]['C'][0]=(unsigned short)Peq['C'];
    tmp_Peq[thread_id]['G'][0]=(unsigned short)Peq['G'];
    tmp_Peq[thread_id]['T'][0]=(unsigned short)Peq['T'];

    Word Mask_Pre=(Word)1<<(band_length-2);
    Word Mask=(Word)1<<(band_length-1);
    Word VP=0;
    Word VN=0;
    Word X=0;
    Word D0=0;
    Word HN=0;
    Word HP=0;


    int s=0;
    int i = 0;
    int j=0;

    int bound=band_length-2-band_down;
    int err=0;

    Word err_mask=(Word)1;
    int s1=band_length-2;
    int i_bd=i+band_down;
    int last_high=band_length-t_length+p_length-band_down-1;
   int t_length_1=t_length-1;
    //while(i<t_length)
    while(i<t_length_1)
    {

        X=Peq[text[i]]|VN;

        D0=((VP+(X&VP))^VP)|X;

        HN=VP&D0;
        HP=VN|~(VP|D0);

        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);
        if(!(D0&err_mask))
         {
            ++err;
            if((err-last_high)>errthold)
                return -1;
        }

        Peq['A']=Peq['A']>>1;
        Peq['C']=Peq['C']>>1;
        Peq['G']=Peq['G']>>1;
        Peq['T']=Peq['T']>>1;


        ++i;
        ++i_bd;
        Peq[pattern[i_bd]]=Peq[pattern[i_bd]]|Mask;
    }


    ///这个循环拿出来是为了防止内存泄露

       X=Peq[text[i]]|VN;
        D0=((VP+(X&VP))^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);
        if(!(D0&err_mask))
         {
            ++err;
            if((err-last_high)>errthold)
                return -1;
        }


    int site=p_length-last_high-1;
    int return_site=-1;
    if((err<=errthold)&&(err<*return_err))
    {
            *return_err=err;
            return_site=site;
    }
    int i_last=i;
    i=0;
    while(i<last_high)
    {
        err=err+((VP>>i)&(Word)1);
        err=err-((VN>>i)&(Word)1);
        ++i;

        if((err<=errthold)&&(err<*return_err))
        {
                *return_err=err;
                return_site=site+i;
        }


    }
    return return_site;

}

/*************************************************************************************************************************************************************************/



/**************************************************单条比对***************************************************************************************************/
int Start_location_Reserve_Banded_BPM_high(char *pattern,int p_length,char *text,int t_length,
                unsigned short errthold,unsigned short band_down,unsigned short band_below,unsigned short band_length,int* return_err, int thread_id)
{

    Word Peq[128];
    char Peq_index[4]= {'A','C','G','T'};
    int symbol = 0;
    int r;
    Word tmp_Peq_1=(Word)1;


    Peq['A']=(Word)0;
    Peq['T']=(Word)0;
    Peq['G']=(Word)0;
    Peq['C']=(Word)0;
    Word Peq_A;
    Word Peq_T;
    Word Peq_C;
    Word Peq_G;

    for (r =0; r<band_length; r++)
    {
        Peq[pattern[p_length-r-1]]=Peq[pattern[p_length-r-1]]|tmp_Peq_1;
        tmp_Peq_1=tmp_Peq_1<<1;
    }
    Peq_A=Peq['A'];
    Peq_C=Peq['C'];
    Peq_T=Peq['T'];
    Peq_G=Peq['G'];

    for(symbol = 0; symbol < 128; symbol++)
    {
        Peq[symbol]=(Word)0;
    }
    Peq['A']=Peq_A;
    Peq['C']=Peq_C;
    Peq['T']=Peq_T;
    Peq['G']=Peq_G;
    tmp_Peq_32[thread_id]['A'][0]=(unsigned int)Peq['A'];
    tmp_Peq_32[thread_id]['C'][0]=(unsigned int)Peq['C'];
    tmp_Peq_32[thread_id]['G'][0]=(unsigned int)Peq['G'];
    tmp_Peq_32[thread_id]['T'][0]=(unsigned int)Peq['T'];



    Word Mask_Pre=(Word)1<<(band_length-2);
    Word Mask=(Word)1<<(band_length-1);
    Word VP=0;
    Word VN=0;
    Word X=0;
    Word D0=0;
    Word HN=0;
    Word HP=0;
    int s=0;
    int i = 0;
    int j=0;

    int bound=band_length-2-band_down;
    int err=0;

    Word err_mask=(Word)1;
    int s1=band_length-2;
    i=t_length-1;
    int i_bd=p_length-band_length;
    int last_high=band_length-t_length+p_length-band_down-1;
    while(i>0)
    {
        X=Peq[text[i]]|VN;
        D0=((VP+(X&VP))^VP)|X;
        HN=VP&D0;
        HP=VN|~(VP|D0);
        X=D0>>1;
        VN=X&HP;
        VP=HN|~(X|HP);
        if(!(D0&err_mask))
         {
            ++err;
            if((err-last_high)>errthold)
                return -1;
        }

        Peq['A']=Peq['A']>>1;
        Peq['T']=Peq['T']>>1;
        Peq['G']=Peq['G']>>1;
        Peq['C']=Peq['C']>>1;



        --i;
        --i_bd;
        Peq[pattern[i_bd]]=Peq[pattern[i_bd]]|Mask;
    }

    X=Peq[text[i]]|VN;
    D0=((VP+(X&VP))^VP)|X;
    HN=VP&D0;
    HP=VN|~(VP|D0);
    X=D0>>1;
    VN=X&HP;
    VP=HN|~(X|HP);
    if(!(D0&err_mask))
     {
        ++err;
        if((err-last_high)>errthold)
            return -1;
    }

    int site=last_high;
    int return_site=-1;
    if((err<=errthold)&&(err<*return_err))
    {
            *return_err=err;
            return_site=site;
    }
    int i_last=i;
    i=0;
    while(i<last_high)
    {
        err=err+((VP>>i)&(Word)1);
        err=err-((VN>>i)&(Word)1);
        ++i;

        if((err<=errthold)&&(err<*return_err))
        {
            *return_err=err;
            return_site=site-i;
        }
    }

    return return_site;

}

/*************************************************************************************************************************************************************************/


