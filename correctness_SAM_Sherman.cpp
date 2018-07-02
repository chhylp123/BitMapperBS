#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include <boost/algorithm/string.hpp>

#include<algorithm>
#include<ctime>

using namespace boost::algorithm;
using namespace std;


typedef unsigned long int _uli;
typedef unsigned int _ui;
typedef unsigned short int _usi;

class Read{
public:
Read() : mism(0), indel_len(0), chr("0"), strand('+'), pos(0) {};
Read(string ch, char st, long int p, int m) : indel_len(0), mism(m), chr(ch), strand(st), pos(p) {};
string chr;
char strand;
long int pos;
int mism;
int indel_len;
};

void make_mask(vector<_uli> &mask){
	_uli SIZE = 64;
        _uli one = 1;
        for(int i = 0; i <= SIZE; i++){
                _uli res = (one << i);
                mask.push_back(res);
        }//for
}//make mask

void count_mism(const string & cigar, const string &md, int & mism, int &indel_len){
	mism = 0; 
	indel_len = 0;
	//break cigar into chunks ending with character
	vector<int> indices;
	int asize = cigar.length();
	for(int i = 0; i < asize; i++){
		if(!isdigit(cigar[i]))
			indices.push_back(i);	
	}//for	
	int ind_size = indices.size();
	//check if there is D, deletion 
	bool isD = false;
	vector<int> deletions;

	for(int i = 0; i < ind_size; i++){
		if(cigar[indices[i]] == 'D' || cigar[indices[i]] == 'I'){
			int cur_index = indices[i];
			//find deletion length
			int prev_ind = 0;
			if(i > 0){
				prev_ind = indices[i - 1] + 1;//index of indel len
			}
				string sub = cigar.substr(prev_ind, cur_index - prev_ind);
				istringstream istr(sub);
				int alen ;
				istr >> alen;
				indel_len += alen;
				if(cigar[indices[i]] == 'D')
					deletions.push_back(alen);
		}//if
	}//for

	int j = 0, del_ind = 0;
	asize = md.length();
	while(j < asize){
		
		if(md[j] == '^'){
			int alen = deletions[del_ind];
			del_ind++;
			j += alen;	
		}//if deletion
		else{
			if(!isdigit(md[j]))
				mism++;
		}
		j++;
	}//scan md
}//count_mism

void adjust_pos(long int &pos, const string &cigar){
	//find clip_left
	//find the first H character
	int i = 0;
	int asize = cigar.length();
	while(i < asize){
		if(!isdigit(cigar[i]))
			break;
		i++;
	}//while

	if(i < asize){
		//check what is the first character
		if(cigar[i] == 'H' || cigar[i] == 'S'){
			string sub = cigar.substr(0, i);
			istringstream istr(sub);
			int clip_left ;
			istr >> clip_left;
			pos -= clip_left;//in BS-seeker
		}//if
	}
}//adjust_pos

void parse_options(int argc, char* argv[], string &inp_pos, string &outp, string &totals, 
	string &pos_thresh)
{
	int res = 0;

	long int i, j;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-'){
			cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; 
			exit(0); 
			return ;
		}
		switch(argv[i][1]){
		case 'r': outp = argv[++i]; break;//output of brat-bw
		case 'p': inp_pos = argv[++i]; break;//original positions
		case 't': totals = argv[++i]; break;//total reads attempted to map
		case 'e': pos_thresh = argv[++i]; break;//threshold for position to be considered correct
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; exit(0); break;
		}
	}//for i

	return ;
}//parse_opt()

int main(int argc, char* argv[]){
	string inp ;//file with correct positions: read_id, read, chr, strand, pos, mism, indel_len, indel_pos
	string inp2 ;//output of brat-bw
	string totals ;//total reads attempted to ma
	string thresh_pos ;

	parse_options(argc, argv, inp, inp2, totals, thresh_pos);
	long int thresh = 50;
	if(!thresh_pos.empty()){
		istringstream istr_pos(thresh_pos);
		istr_pos >> thresh ;
	}


	
	long int total = 100;
	/**
	if(totals.empty()){
		cout << "ERROR: please provide the total number of reads attempted to map (option t)" << endl;
		exit(0);
	}
	else{
		istringstream istr(totals);
		istr >> total;
	}
	**/


	ifstream in;
	in.open(inp.c_str(), ios::in);
	
	if(!in){
		cout << "ERROR: please provide file with original positions (option p)" << endl;
		exit(0);
	}


	string read, chr;
	int strand;
	int indel_len, indel_pos, mism;
	long int id, st = 0, en = 0, pos, last, error ;
	string id_string, tmp_string;
	vector<Read> reads;//max len 76
	reads.resize(total+2);

	///vector<string> fields;

	long int first_i, second_i, third_i;

	string o_chr;
	string o_position;
	string o_err;

	long int o_position_p;

	long long total_good_chhy=0;
	long long total_equal_chhy=0;
	long long total_less_chhy=0;

	string path;

	string line;

	string read_sequence;

	string read_qual;

	string tmp;

	string err;

	int maq_socre;

	int err_int;
	int o_err_int;

	in >> id_string;	
	while(!in.eof()){



		if(id_string[0] != '@')
		{
			///cout<<id_string<<endl;
			for (int i = 0; i < id_string.length(); ++i)
			{
				if(id_string[i]=='_')
				{
					first_i = i;
				}

				if(id_string[i]==':')
				{
					second_i = i;
				}

				if(id_string[i]=='-')
				{
					third_i = i;

					break;
				}

			}

			o_chr = id_string.substr(first_i + 1, second_i - first_i - 1);
			o_position = id_string.substr(second_i + 1, third_i - second_i - 1);
			o_err = id_string.substr(third_i + 1);

			istringstream istr100(o_position);
			istr100 >> o_position_p; 




			///in >> chr >> strand >> pos >> error>>path;
			///in >> strand >> chr >> pos;
			in >> strand >> chr >> pos >> tmp >> path; 




			getline(in, line);
			
			/**
			cout<<o_chr<<endl;
			cout<<o_position_p<<endl;
			cout<<chr<<endl;
			cout<<strand<<endl;
			cout<<pos<<endl;
			cout<<error<<endl;
			**/
			



			if(o_chr == chr && (labs(o_position_p - pos) <= 20))
			{

				total_equal_chhy++;


			}
			else
			{
				/**
				if(o_chr == chr && (labs(o_position_p - pos) <= 10))
				{
					cerr<<"id: "<<id_string<<endl;
					cerr<<"o_chr: "<<o_chr<<endl;
					cerr<<"o_position_p: "<<o_position_p<<endl;
					

					cerr<<"chr: "<<chr<<endl;
					cerr<<"pos: "<<pos<<endl;
					cerr<<"strand: "<<strand<<endl;
					cerr<<"path: "<<path<<endl;

					cerr<<"*******************************************"<<endl;
				}
				**/
				
				
			}

		}
		else
		{
			
			getline(in, line);
		}

		in >> id_string;

	
	}//while
	in.close(); in.clear();


	cout << "\n total_equal_chhy = " << total_equal_chhy << endl;
	cout << "\n total_good_chhy = " << total_good_chhy << endl;
	cout << "\n total_less_chhy = " << total_less_chhy << endl;


	return 1;






	/**
    31-mer
    bitmapper_BS相同或者更好: 0.95818380651052173359097955110101
    BSMAP更好: 0.03297927554532125694876886358956

    25-mer
    bitmapper_BS相同或者更好: 0.9607743016169390700933679473422
    BSMAP更好: 0.02905124891419699466518255070743

    20-mer
    bitmapper_BS相同或者更好: 0.96334454979977767890255748274444
    BSMAP更好: 0.0256442101946717010724449740469
    **/




    /**
    17735453 - 17709467 = 25986


    match_read = 17709467
    unmatch_read = 1252993
    empty_read = 68253


	empty_read =42267
	matched_read =17735453
	unmatched_read =1252993


    535053
    18962460
    19030713

	bitmapper_BS > bsmap: 98.6%

    **/



	vector<_uli> mask;
	make_mask(mask);

	long int total_pos_good = 0;
	in.open(inp2.c_str(), ios::in);

	if(!in){
		cout << "ERROR: please provide output file in SAM format (option p)" << endl;
		exit(0);
	}

	long int total_wrong = 0, total_good = 0, wrong_better = 0, wrong_equal = 0, wrong_worse = 0;
	_uli flag = 0, zero = 0;
	int temp;
	string idd, qual, cigar, md, dum;
        
	long int total_chr = 0, total_strand = 0, total_pos = 0, total_pos_delta = 0;
	long int counter = 0;
	in >> idd;//id starts with 0 and pos starts with 0 
    while(!in.eof()){

    	///cerr << "idd: "<< idd << endl;

		///if(idd[0] != '@'){
    	    idd = idd.substr(1, idd.length());
			istringstream istr2(idd);
			istr2 >> id; 
			counter++;
			/**
			if(counter > total + 2)
			{
				cout << "ERROR counter > total" << endl;
				break;
			}
			**/                
			///in >> flag >> chr >>  pos >> temp >> cigar >> dum >> temp >> temp >> read >> qual ;

			in >> chr >> strand >> pos;

            /**
			cerr << "id: "<< id << endl;
        	cerr << "chrom: "<< chr << endl;
       	 	cerr << "pos: "<<pos << endl;
        	cerr << "strand: "<<strand << endl;
            **/
        	
        	//if(reads[id].chr == chr && reads[id].strand == strand && reads[id].pos == pos)
        	if(reads[id].chr == chr && 
        		(labs(reads[id].pos - pos) < 5))
			{

				total_good++;

			}
			else
			{
				cerr << reads[id].pos - pos <<endl;
				/**
				cerr << "***********************************************"<< endl;
				cerr << "origin chrom: "<< reads[id].chr << endl;
				cerr << "origin strand: "<< reads[id].strand << endl;
				cerr << "origin pos: "<< reads[id].pos << endl;

				cerr << "chrom: "<< chr << endl;
				cerr << "strand: "<<strand << endl;
       	 		cerr << "pos: "<<pos << endl;
       	 		**/
        		
			}	
			

			//in >> flag >> chr >>  pos >> temp >> cigar >> dum >> temp >> temp >> read >> qual ;
            
            /**
			getline(in, dum);
			adjust_pos(pos, cigar);		

			strand = '+';
			if((mask[4] & flag) > zero)
				strand = '-'; 
			if(reads[id].chr == chr && reads[id].strand == strand && reads[id].pos == pos)
			{		
				total_good++;
				if(strand == '+')
					total_pos_good++;	
			}//if	
			else
			{	total_wrong++;
				if(reads[id].chr != chr)
					total_chr++;
				else if(reads[id].strand != strand)
					total_strand++;
				else{			
				
					long int minus_pos = pos - thresh;
					long int pos_plus = pos + thresh;
					if(reads[id].pos >= minus_pos && reads[id].pos <= pos_plus)
						total_pos_delta++;
					else 
						total_pos++;
				}
			}//else
            **/
		/**
		}//if not @
		
		else{
			string line;
			getline(in, line);
		}
		**/
		in >> idd;

    }//while
	cout << "Total of correctly mapped reads = " << total_good << ", of them positive: " 
		<< total_pos_good <<  endl;

    long int total_reads = total_good + total_wrong;
    long int total_correct2 = total_reads - (total_chr + total_strand + total_pos);
	long int total_wrong2 = total_chr + total_strand + total_pos;
    double map_acc1 = (0.0 + total_good)/(0.0 + total_good + total_wrong);
    double map_acc2 = (0.0 + total_correct2)/(0.0 + total_reads);

    cout << "Mapping accuracy = " << (100.0 * map_acc2) << endl;
 
	in.close();
}
