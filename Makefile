CC=g++
CFLAGS = -w -c -mavx2 -mpopcnt -fomit-frame-pointer -W -Wall -Winline -O3
#LDFLAGS = -lz -lm 
LDFLAGS = -lm 
SOURCES = saca-k.cpp bwt.cpp uint40.h Bitmapper_main.cpp Process_CommandLines.cpp Auxiliary.cpp Index.cpp Schema.cpp Process_sam_out.cpp Process_Reads.cpp Ref_Genome.cpp Levenshtein_Cal.cpp SAM_queue.cpp
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = bitmapperBS


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o
	chmod 777 psascan

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) -lpthread -O3 -mpopcnt -mavx2

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ \#* bitmapper_BS
#relase:
	#chmod 777 psascan


