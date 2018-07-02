#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <time.h>

using namespace std;

///g++ -std=c++14 -O3 -DNDEBUG -msse4.2 -I ~/include -L ~/lib csa-index-hyb-8.cpp -o csa-index-hyb-8.x -lsdsl -ldivsufsort -ldivsufsort64 -msse4.2

int main(int argc, char** argv)
{

    string input_file_name  = string(argv[1]);

    ifstream input_file(argv[1]);

    if(!input_file)
    {  
        cout << "Unable to open input_file";  
        exit(1); // terminate with error  
    }





    char buffer;

    long long i = 0;



    while (!input_file.eof())  
    {  
        input_file.get(buffer); 

        if (buffer!=' ' && buffer!='\n' && buffer != '\0')
         {

            if (buffer!='A'&&
                buffer!='C'&&
                buffer!='G'&&
                buffer!='T'&&
                buffer!='a'&&
                buffer!='c'&&
                buffer!='g'&&
                buffer!='t')
            {
                cout<<"Text has non-A-C-G-T characters!"<<endl;
                cout<<"character = "<<buffer<<endl;
                cout<<"ASCII code = "<<(unsigned int)buffer<<endl;
                cout<<"i = "<<i<<endl;

                return 1;
            }

            if (buffer == 'C' 
                ||
                buffer == 'c')
            {
                buffer = 'T';
            }

             
             cerr<<buffer;
         } 

         i++;
    }    

    
}





