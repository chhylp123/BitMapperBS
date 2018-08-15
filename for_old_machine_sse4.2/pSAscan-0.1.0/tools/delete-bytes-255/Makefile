SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -march=native

all: delete255

delete255:
	$(CC) $(CFLAGS) -o delete255 main.cpp

clean:
	/bin/rm -f delete255 *.o
