CC=g++
CFLAGS = -w -c -msse4.2 -mpopcnt -fomit-frame-pointer -W -Wall -Winline -O3
#LDFLAGS = -lz -lm 
LDFLAGS = -lm 
SOURCES = saca-k.cpp bwt.cpp uint40.h Bitmapper_main.cpp Process_CommandLines.cpp Auxiliary.cpp Index.cpp Schema.cpp Process_sam_out.cpp Process_Reads.cpp Ref_Genome.cpp Levenshtein_Cal.cpp SAM_queue.cpp
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = bitmapperBS


.ONESHELL:
all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o; cd libdivsufsort-2.0.1/; cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="/usr/local" .; make; sudo make install; cd ../pSAscan-0.1.0/src; make; cp psascan ../../; chmod 777 psascan
	#rm -rf *.o
	#cd libdivsufsort-2.0.1/
	#cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="/usr/local" .
	#make
	#sudo make install
	#cd ../pSAscan-0.1.0/src 
	#make
	#cp psascan ../../

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) -lpthread -O3 -mpopcnt -msse4.2

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ \#* bitmapper_BS
#relase:
	#chmod 777 psascan


