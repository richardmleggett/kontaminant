BIN = bin
MAXK = 31
FLAGBITS = 4
CFIELDS = 0

# On older versions of XCode, it was necessary to include the following: -fnested-functions -L/opt/local/lib/

# Calculate bitfields needed to store kmer
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

OPT	= -Wall -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -DFLAG_BITS_USED=$(FLAGBITS) -DCONTAMINANT_FIELDS=$(CFIELDS) -pthread -O3

KONTAMINANT_OBJ = obj/kontaminant.o obj/hash_table.o obj/hash_value.o obj/logger.o obj/binary_kmer.o obj/element.o obj/kmer_reader.o obj/cmd_line.o obj/seq.o obj/kmer_stats.o obj/kmer_build.o

all:remove_objects $(KONTAMINANT_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) -lm -o $(BIN)/kontaminant $(KONTAMINANT_OBJ)

clean:
	rm obj/*
	rm -rf $(BIN)/kontaminant

remove_objects:
	rm -rf obj

obj/%.o : src/%.c
	mkdir -p obj; $(CC) -Iinclude $(OPT) -c $< -o $@

