#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include<iostream>
typedef uint64_t Word;



#define SEQ_MAX_LENGTH		1000			// Seq Max Length








/**************************************************单条比对***************************************************************************************************/
inline int Calculate_Cigar(
	char *pattern, int p_length, 
	char *text, int t_length, 
	unsigned short errthold, 
	int* return_err, 
        char* cigar,
		uint16_t* matrix,
        int end_site,
		int* return_start_site,
        char* path
        )
{
	///(*return_err) = 999999;

	///这个是那个需要预处理的向量
	Word Peq[256];

	int band_length = (errthold << 1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;

	


	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///这是把pattern的前2k + 1个字符预处理
	///pattern[0]对应Peq[0]
	///pattern[2k]对应Peq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

	Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);
	memset(matrix, 0, sizeof(uint16_t)*band_length);

	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = 1 << (errthold << 1);
	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	Word inner_i;
	Word column_start;

	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)
	while (i<t_length_1)
	{
		///pattern[0]在Peq[2k], 而pattern[2k]在Peq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);


		///pattern[0]在Peq[2k], 而pattern[2k]在Peq[0]
		//右移实际上是把pattern[0]移掉了
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///这是把新的pattern[2k]加进来, 这貌似是加到Peq[2k]上了
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;

		Peq['T'] = Peq['T'] | Peq['C'];


		column_start = i*band_length;

		for (inner_i = 0; inner_i < band_length; inner_i++)
		{
			matrix[column_start + inner_i] = matrix[column_start - band_length + inner_i] + (~(D0 >> inner_i)&err_mask);
		}
		
	}


	///这个循环拿出来是为了防止内存泄露
	///其实也就是循环里的最后一行语句吧
	///完全可以把pattern增大一位
	///不过这样也好，可以减少计算开销
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);

	///注意这里是(i+1)
	column_start = (i+1)*band_length;

	for (inner_i = 0; inner_i < band_length; inner_i++)
	{
		matrix[column_start + inner_i] = matrix[column_start - band_length + inner_i] + (~(D0 >> inner_i)&err_mask);
	}



        
	std::cout << "Matrix: "<< std::endl;

    int back_track_site = band_length - (p_length - end_site);

    ///char cigar[1000];

    Word v_value, h_value, delta_value, min_value;
    Word direction, is_mismatch; ///0 is match, 1 is mismatch, 2 is xiang shang, 3 is shuiping

    i = t_length;
    int path_length=0;
    int start_site = end_site;



    while(i>0)
    {

        column_start = i * band_length;

        h_value = matrix[column_start - band_length + back_track_site + 1];
        delta_value = matrix[column_start - band_length + back_track_site];

        ///Only this situation has xiangshang keneng
        if(back_track_site!=0)
        {
             

            v_value = matrix[column_start + back_track_site - 1];

            min_value = delta_value;
            direction = 0;

            if(v_value < min_value)
            {
            min_value = v_value;
            direction = 2;
            }


            if(h_value < min_value)
            {
            min_value = h_value;
            direction = 3;
            }

        }
        else
        {
            min_value = delta_value;
            direction = 0;

            if(h_value < min_value)
            {
            min_value = h_value;
            direction = 3;
            }
        }

           
        if(direction==0)
        {

            if(delta_value != matrix[column_start + back_track_site])
            { 
            direction = 1;
            }

            i--;

            start_site--;
   
        }
        if(direction==2)///ru guo xiang shang yi dong, bing bu huan lie
        {
            back_track_site--;
            start_site--;
        }
        else if(direction==3)///ru guo xiang zuo yi dong
        {
            i--;
            back_track_site++;
        }

        path[path_length++] = direction;

    }

    if(direction!=3)
    {
        start_site++;
    }
	(*return_start_site) = start_site;



    std::cout<<"start_site: "<<start_site<<std::endl;
    std::cout<<"end_site  : "<<end_site<<std::endl;

	std::cout << "Direction: ";

    for(i=path_length-1; i>=0; i--)
    {
        std::cout<<(int)path[i];
    }

	std::cout << std::endl;

	
}













/**************************************************单条比对***************************************************************************************************/
inline int BS_Reserve_Banded_BPM(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, int* return_err)
{
	(*return_err) = 999999;

	///这个是那个需要预处理的向量
	Word Peq[256];

	int band_length = (errthold<<1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;




	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///这是把pattern的前2k + 1个字符预处理
	///pattern[0]对应Peq[0]
	///pattern[2k]对应Peq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

    Peq['T'] = Peq['T'] | Peq['C'];

	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = 1 << (errthold << 1);
	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)
	while (i<t_length_1)
	{
		///pattern[0]在Peq[2k], 而pattern[2k]在Peq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);
		///如果斜对角线方向匹配则D0是1
		///如果不匹配则D0是0
		///这个意思是如果最上面那条对角线上的斜对角线方向发生误配,则执行内部程序
		///
		if (!(D0&err_mask))
		{
			++err;

			///即使全部递减，也就减2k
			if ((err - last_high)>errthold)
				return -1;
		}

		///pattern[0]在Peq[2k], 而pattern[2k]在Peq[0]
		//右移实际上是把pattern[0]移掉了
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///这是把新的pattern[2k]加进来, 这貌似是加到Peq[2k]上了
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;

        Peq['T'] = Peq['T'] | Peq['C'];
	}


	///这个循环拿出来是为了防止内存泄露
	///其实也就是循环里的最后一行语句吧
	///完全可以把pattern增大一位
	///不过这样也好，可以减少计算开销
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);
	if (!(D0&err_mask))
	{
		++err;
		if ((err - last_high)>errthold)
			return -1;
	}

	/// last_high = 2k
	/// site = (SEQ_LENGTH + 2k) - 2k -1
	/// site = SEQ_LENGTH - 1
	///此时这个site貌似是最上面那条对角线的位置
	///int site = p_length - last_high - 1;
	int site = t_length - 1;
	int return_site = -1;
	if ((err <= errthold) && (err<*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	int i_last = i;
	i = 0;
	while (i<last_high)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err<*return_err))
		{
			*return_err = err;
			return_site = site + i;
		}


	}
	return return_site;

}







inline int Reserve_Banded_BPM(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, int* return_err)
{
	///这个是那个需要预处理的向量
	Word Peq[256];

	int band_length = (errthold<<1) + 1;
	int i = 0;
	Word tmp_Peq_1 = (Word)1;




	Peq['A'] = (Word)0;
	Peq['T'] = (Word)0;
	Peq['G'] = (Word)0;
	Peq['C'] = (Word)0;


	Word Peq_A;
	Word Peq_T;
	Word Peq_C;
	Word Peq_G;

	///band_length = 2k + 1
	///这是把pattern的前2k + 1个字符预处理
	///pattern[0]对应Peq[0]
	///pattern[2k]对应Peq[2k]
	for (i = 0; i<band_length; i++)
	{
		Peq[pattern[i]] = Peq[pattern[i]] | tmp_Peq_1;
		tmp_Peq_1 = tmp_Peq_1 << 1;
	}

        


	Peq_A = Peq['A'];
	Peq_C = Peq['C'];
	Peq_T = Peq['T'];
	Peq_G = Peq['G'];


	memset(Peq, 0, sizeof(Word)* 256);


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = 1 << (errthold << 1);
	Word VP = 0;
	Word VN = 0;
	Word X = 0;
	Word D0 = 0;
	Word HN = 0;
	Word HP = 0;


	i = 0;

	int err = 0;

	Word err_mask = (Word)1;


	///band_down = 2k
	///i_bd = 2k
	///int i_bd = i + band_down;
	int i_bd = (errthold << 1);


	int last_high = (errthold << 1);


	/// t_length_1 = SEQ_LENGTH - 1
	int t_length_1 = t_length - 1;
	//while(i<t_length)
	while (i<t_length_1)
	{
		///pattern[0]在Peq[2k], 而pattern[2k]在Peq[0]
		X = Peq[text[i]] | VN;

		D0 = ((VP + (X&VP)) ^ VP) | X;

		HN = VP&D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X&HP;
		VP = HN | ~(X | HP);
		///如果斜对角线方向匹配则D0是1
		///如果不匹配则D0是0
		///这个意思是如果最上面那条对角线上的斜对角线方向发生误配,则执行内部程序
		///
		if (!(D0&err_mask))
		{
			++err;

			///即使全部递减，也就减2k
			if ((err - last_high)>errthold)
				return -1;
		}

		///pattern[0]在Peq[2k], 而pattern[2k]在Peq[0]
		//右移实际上是把pattern[0]移掉了
		Peq['A'] = Peq['A'] >> 1;
		Peq['C'] = Peq['C'] >> 1;
		Peq['G'] = Peq['G'] >> 1;
		Peq['T'] = Peq['T'] >> 1;


		++i;
		++i_bd;
		///这是把新的pattern[2k]加进来, 这貌似是加到Peq[2k]上了
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;


	}


	///这个循环拿出来是为了防止内存泄露
	///其实也就是循环里的最后一行语句吧
	///完全可以把pattern增大一位
	///不过这样也好，可以减少计算开销
	X = Peq[text[i]] | VN;
	D0 = ((VP + (X&VP)) ^ VP) | X;
	HN = VP&D0;
	HP = VN | ~(VP | D0);
	X = D0 >> 1;
	VN = X&HP;
	VP = HN | ~(X | HP);
	if (!(D0&err_mask))
	{
		++err;
		if ((err - last_high)>errthold)
			return -1;
	}

	/// last_high = 2k
	/// site = (SEQ_LENGTH + 2k) - 2k -1
	/// site = SEQ_LENGTH - 1
	///此时这个site貌似是最上面那条对角线的位置
	///int site = p_length - last_high - 1;
	int site = t_length - 1;
	int return_site = -1;
	if ((err <= errthold) && (err<*return_err))
	{
		*return_err = err;
		return_site = site;
	}
	int i_last = i;
	i = 0;
	while (i<last_high)
	{
		err = err + ((VP >> i)&(Word)1);
		err = err - ((VN >> i)&(Word)1);
		++i;

		if ((err <= errthold) && (err<*return_err))
		{
			*return_err = err;
			return_site = site + i;
		}


	}
	return return_site;

}








/*************************************************************************************************************************************************************************/


int main(int argc, char *argv[])
{

  ///char* pattern = (char*)malloc(sizeof(char)*1000);
  char* pattern = argv[1];
  int p_length; 
  ///char* text = (char*)malloc(sizeof(char)*1000);
  char* text = argv[2];
  int t_length;
  int errthold;
  int return_err;
  int site;


  
  ///pattern should longer than text
  /**
  std::cout<<"Please input pattern length:"<<std::endl;
  std::cin>>p_length;
  **/
  /**
  std::cout<<"Please input pattern:"<<std::endl;
  std::cin>>pattern;
  **/
  p_length = strlen(pattern);
  t_length = strlen(text);

  /**
  std::cout<<"Please input text length:"<<std::endl;
  std::cin>>t_length;
  **/
  /**
  std::cout<<"Please input text:"<<std::endl;
  std::cin>>text;
  **/

  std::cout<<"Pattern:"<<pattern<<std::endl;
  std::cout<<"Pattern length:"<<p_length<<std::endl;
  std::cout<<"Text   :"<<text<<std::endl;
  std::cout<<"Text length:"<<t_length<<std::endl;

  std::cout<<"Please errthold:"<<std::endl;
  std::cin>>errthold;

  ///p_length = t_length + 2*errthold
  if(p_length != t_length + 2*errthold)
  {
    std::cout<<"ERROR: p_length != t_length + 2*errthold"<<std::endl;
    exit(1);
  }

  std::cout<<"Pattern: "<<pattern<<std::endl;
  std::cout<<"Text   : "<<text<<std::endl;

  ///return_err = 9999;
  site = BS_Reserve_Banded_BPM(pattern, p_length, text, t_length, errthold, &return_err);
  
  std::cout<<"Distance: "<<return_err<<std::endl;
  std::cout<<"site    : "<<site<<std::endl;


  uint16_t* matrix = (uint16_t*)malloc(sizeof(uint16_t)*(errthold*2+1)*(1000));

  char* cigar= (char*)malloc(sizeof(char)*1000);

  char* path = (char*)malloc(sizeof(char)*1000);

  int start_site;

  Calculate_Cigar(pattern, p_length, text, t_length, errthold, &return_err, cigar, matrix, site, &start_site, path);





}
