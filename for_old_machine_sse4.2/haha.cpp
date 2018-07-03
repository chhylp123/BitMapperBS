///这个程序实际就假设pattern里面有可能有N
inline int BS_Reserve_Banded_BPM
(char *pattern, int p_length, char *text, int t_length, unsigned short errthold, int* return_err)
{
	(*return_err) = 999999;

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


	Peq['A'] = Peq_A;
	Peq['C'] = Peq_C;
	Peq['T'] = Peq_T;
	Peq['G'] = Peq_G;


	Word Mask = ((Word)1 << (errthold << 1));

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


	///fprintf(stderr, "sucess(1)\n");


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


	////fprintf(stderr, "sucess(2)\n");

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



