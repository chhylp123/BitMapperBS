SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3 -march=native

all: psascan
psascan: main.cpp
	$(CC) $(CFLAGS) -o psascan ./psascan_src/utils.cpp main.cpp -ldivsufsort -ldivsufsort64 -fopenmp
clean:
	/bin/rm -f *.o
nuclear:
	/bin/rm -f psascan *.o
