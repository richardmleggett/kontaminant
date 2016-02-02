#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <assert.h>
#include "global.h"
#include "binary_kmer.h"
#include "element.h"
#include "hash_table.h"
#include "cmd_line.h"
#include "kmer_stats.h"
#include "kmer_reader.h"

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void load_reads_into_table(CmdLine* cmd_line,  HashTable* kmer_hash)
{
    KmerFileReaderArgs fra;
    fra.bad_reads = 0;
    fra.colour = 0;
    fra.fastq_ascii_offset = cmd_line->quality_score_offset;
    
    KmerFileReaderWrapperArgs fria;
    fria.format = cmd_line->format;
    fria.kmer_size = kmer_hash->kmer_size;
    fria.max_read_length = 1000;
    fria.new_entry = true;
    fria.full_entry = false;

    fra.input_filename = cmd_line -> input_filename_one;
    fra.quality_cut_off = cmd_line->quality_score_threshold;
    fra.insert = true;
    fra.max_read_length = 1000;
    fra.maximum_ocupancy = 75;
    fra.KmerHash = kmer_hash;
    long long loaded_kmers = 0;
    
    loaded_kmers = load_seq_into_kmers_hash(&fra, &fria);
    
    printf("Loaded %'lld kmers (bad reads %'lld)", loaded_kmers, fra.bad_reads);
    
}

typedef struct {
    long long kmers_dumped;
    FILE* fout;
    int kmer_size;
} PrintNodeBinaryStruct;

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void print_node_binary(Element* node, void* data) {
    PrintNodeBinaryStruct* pnb = (PrintNodeBinaryStruct*)data;
    db_node_print_binary(pnb->fout, node, pnb->kmer_size);
    pnb->kmers_dumped++;
    //printf("k=%d\n", pnb->kmer_size);
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void dump_kmer_hash(CmdLine* cmd_line, HashTable * kmer_hash)
{
    char* output_filename = malloc(strlen(cmd_line->input_filename_one) + 16);
    KmerLibraryHeader* header = calloc(1, sizeof(KmerLibraryHeader));
    PrintNodeBinaryStruct pnb;
    
    pnb.kmers_dumped = 0;
    
    if (!header) {
        printf("Error: can't allocate room for header\n");
        exit(1);
    }
    
    strcpy(header->header_word, "KONTAMINANT");
    strcpy(header->footer_word, "KONTAMINANT");
    header->version = BINVERSION;
    header->kmer_size = kmer_hash->kmer_size;
    header->num_bitfields = NUMBER_OF_BITFIELDS_IN_BINARY_KMER;
    header->num_kmers = (uint32_t)kmer_hash->unique_kmers;

    sprintf(output_filename, "%s.%d.kmers", cmd_line->input_filename_one, cmd_line->kmer_size);
    
    printf("\nDumping hash table to file: %s\n", output_filename);

    pnb.fout = fopen(output_filename, "wb");
	if (pnb.fout == NULL) {
		fprintf(stderr, "Error: cannot open %s", output_filename);
		exit(1);
	}
    
    fwrite(header, sizeof(KmerLibraryHeader), 1, pnb.fout);
    
    hash_table_traverse_with_data(&print_node_binary, (void*)&pnb, kmer_hash);
    
	fclose(pnb.fout);
    fflush(stdout);
	printf("%'lld kmers dumped\n", pnb.kmers_dumped);
}
