BIN = bin
MAXK = 31

ifdef MAC
define n


endef
$(warning $n$nNOTE On Mac OS X, you must compile using GCC. If you have not already done so, please change the line CC= in this Makefile to point to your GCC compiler. $nThis must be a full path, as Apple aliases gcc to point to it's own compiler.$n)

# Change the following to point to your GCC binary
#CC=gcc
CC=/usr/local/Cellar/gcc/4.9.2_1/bin/gcc-4.9

# On older versions of XCode, it was necessary to include the following
#MACFLAG = -fnested-functions -L/opt/local/lib/ 
endif

ifeq ($(MAXK),31)
   BITFIELDS = 1
endif

ifeq ($(MAXK),63)
   BITFIELDS = 2
endif

ifeq ($(MAXK),95)
   BITFIELDS = 3
endif

ifeq ($(MAXK),127)
   BITFIELDS = 4
endif

OPT	= $(ARCH) $(MACFLAG) -Wall -O3 -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -pthread -g

CFLAGS_KONTAMINANT = -Iinclude

KONTAMINANT_OBJ = obj/kontaminant.o obj/hash_table.o obj/hash_value.o obj/logger.o obj/binary_kmer.o obj/element.o obj/kmer_reader.o obj/cmd_line.o obj/seq.o obj/kmer_stats.o obj/kmer_build.o

all:remove_objects $(KONTAMINANT_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) -o $(BIN)/kontaminant $(KONTAMINANT_OBJ) -lm

clean:
	rm obj/*
	rm -rf $(BIN)/kontaminant

remove_objects:
	rm -rf obj

obj/%.o : src/%.c
	mkdir -p obj; $(CC) $(CFLAGS_KONTAMINANT) $(OPT) -c $< -o $@

