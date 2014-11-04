/*----------------------------------------------------------------------*
 * File:    kmer_reader.c                                               *
 * Purpose: Read kmers from binary files and FASTQ                      *
 * Authors: Richard Leggett                                             *
 *          Ricardo Ramirez-Gonzalez                                    *
 *          The Genome Analysis Centre (TGAC), Norwich, UK              *
 *          richard.leggett@tgac.ac.uk    								*
 *----------------------------------------------------------------------*/

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
#include <limits.h>
#include "global.h"
#include "binary_kmer.h"
#include "element.h"
#include "hash_table.h"
#include "cmd_line.h"
#include "kmer_stats.h"
#include "kmer_reader.h"

struct print_binary_args {
    boolean(*condition) (Element* node);
    long long count;
    FILE * fout;
    short kmer_size;
};


/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void initialise_kmer_counts(int n, KmerCounts* counts)
{
    int i;
        
    counts->n_contaminants = n;
    counts->kmers_loaded = 0;
    counts->contaminants_detected = 0;
    for (i=0; i<MAX_CONTAMINANTS; i++) {
        counts->kmers_from_contaminant[i] = 0;
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_hash_load_sliding_windows(Element **previous_node, HashTable* kmer_hash, boolean prev_full_entry, KmerFileReaderArgs* fra, short kmer_size, KmerSlidingWindowSet *windows, int read, KmerStats* stats, KmerCounts* counts)
{
    Element *current_node = NULL;
    BinaryKmer tmp_kmer;
    int i;
    int j;
    
    // For each window...
    for (i = 0; i < windows->nwindows; i++) {
        KmerSlidingWindow *current_window = &(windows->window[i]);
        
        // For each kmer in window...
        for (j = 0; j < current_window->nkmers; j++) {
            boolean found;
            
            // Try and find kmer
            Key key = element_get_key(&(current_window->kmer[j]), kmer_size, &tmp_kmer);
            
            // Add this kmer? Not if screening/filtering, but yes if making library
            if (fra->insert) {
                current_node = hash_table_find_or_insert(key, &found, kmer_hash);
            } else {
                current_node = hash_table_find(key, kmer_hash);
            }
            
            // If we found kmer...
            if (current_node != NULL) {
                if (!(i == 0 && j == 0 && prev_full_entry == false && current_node == *previous_node)) {	// otherwise is the same old last entry
                    int c;
                    
                    if (stats != NULL) {
                        for (c=0; c<counts->n_contaminants; c++) {
                            if (element_get_contaminant_bit(current_node, c) > 0) {
                                if (counts->kmers_from_contaminant[c] == 0) {
                                    counts->contaminants_detected++;
                                }
                                
                                counts->kmers_from_contaminant[c]++;
                                
                                if (current_node->coverage[read] == 0) {
                                    if ((current_node->coverage[0] + current_node->coverage[1]) == 0) {
                                        stats->both_reads->contaminant_kmers_seen[c]++;
                                    }
                                    stats->read[read]->contaminant_kmers_seen[c]++;
                                }
                            }
                        }
                    }
                    
                    current_node->coverage[read]++;
                    counts->kmers_loaded++;
                    
                }
            }
            
            *previous_node = current_node;
            
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
int file_reader_wrapper(KmerFileReaderWrapperArgs* wargs)
{
    int length = 0;
    int offset;
    
    assert(wargs->input_fp != NULL);
    
    switch(wargs->format) {
        case FASTA:
    		offset = 0;
            // If reading long FASTA sequence, have to take into account that buffer may be too small.
            // So, if this is a continuation, we need to shift last kmer to start and offset by k.
    		if (wargs->new_entry == false) {
                shift_last_kmer_to_start_of_sequence(wargs->seq, wargs->seq->length, wargs->kmer_size);
    			offset = wargs->kmer_size;
    		}
    		length =  read_sequence_from_fasta(wargs->input_fp, wargs->seq, wargs->max_read_length, wargs->new_entry, &wargs->full_entry, offset);
            break;
        case FASTQ:
    		length = read_sequence_from_fastq(wargs->input_fp, wargs->seq, wargs->max_read_length);
            wargs->full_entry = true;
            break;
        default:
            fprintf(stderr, "Error: File format not supported yet %d\n", wargs->format);
		    exit(1);
            break;
    }
    
    // Did we get a full entry? In which case it's a new entry next time...
    if (wargs->full_entry) {
        wargs->new_entry = true;
    } else {
        wargs->new_entry = false;
    }
    
    return length;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
KmerFileReaderWrapperArgs* get_kmer_file_reader_wrapper(short kmer_size, KmerFileReaderArgs* fra)
{
    KmerFileReaderWrapperArgs* frw;
    char* filename = fra->input_filename;
    int max_read_length = fra->max_read_length;
    FileFormat format = fra->format;
    char offset = fra->fastq_ascii_offset;
    
    if (filename == NULL) {
        return NULL;
    }

    frw = calloc(1, sizeof(KmerFileReaderWrapperArgs));

    // Read from file or stdin
    if (strcmp(filename, "-") == 0 ) {
        frw->input_fp = stdin;
    } else {
        frw->input_fp = fopen(filename, "r");
    }
    
    frw->format = format;
    frw->full_entry = false;
    frw->new_entry = true;
    frw->kmer_size = kmer_size;
    frw->max_read_length = max_read_length;
    frw->seq = sequence_new(max_read_length, max_read_length, offset);
    
    if (frw->seq == NULL) {
        fprintf(stderr, "Error: Unable to allocate Sequence.\n");
        exit(-1);
    }
    
    if (frw->input_fp == NULL) {
        fprintf(stderr, "Error: Unable to open file %s\n", filename);
        exit(-1);
    }
    
    return frw;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
long long screen_kmers_from_file(KmerFileReaderArgs* fra, CmdLine* cmd_line, KmerStats* stats)
{
    HashTable* kmer_hash;
	long long seq_length = 0;
	int entry_length;
	Element *previous_node = NULL;
    boolean keep_reading = true;
	boolean prev_full_entry = true;
    KmerCounts counts;
    KmerFileReaderWrapperArgs* frw;
    KmerSlidingWindowSet* windows;
    KmerStatsReadCounts tmp_stats;
    time_t time_previous = 0;
    time_t time_now = 0;

    assert(fra != NULL);
    assert(fra->KmerHash != NULL);
    
    kmer_hash = fra->KmerHash;
    frw = get_kmer_file_reader_wrapper(kmer_hash->kmer_size, fra);
    initialise_kmer_counts(stats->n_contaminants, &counts);
    
    // Allocate new sliding window structure
    windows = binary_kmer_sliding_window_set_new_from_read_length(kmer_hash->kmer_size, fra->max_read_length);
    
    // Keep reading...
	while ((entry_length = file_reader_wrapper(frw)) && keep_reading)
	{
		int nkmers;
        
        // Update length read
		seq_length += (long long)(entry_length - (prev_full_entry == false ? kmer_hash->kmer_size : 0));
        
        // Get sliding windows
		nkmers = get_sliding_windows_from_sequence(frw->seq->seq, frw->seq->qual, entry_length, fra->quality_cut_off, kmer_hash->kmer_size, windows, windows->max_nwindows, windows->max_kmers, false, 0);
		if (nkmers == 0) {
            // Bad read
            fra->bad_reads++;
		} else {
            // Load kmers
            kmer_hash_load_sliding_windows(&previous_node, kmer_hash, prev_full_entry, fra, kmer_hash->kmer_size, windows, 0, stats, &counts);
        }
        
        if (frw->full_entry == false) {
            // If we didn't get a full entry, then shift last kmer to start of sequence...
            shift_last_kmer_to_start_of_sequence(frw->seq, entry_length, kmer_hash->kmer_size);
        } else {
            hash_table_add_number_of_reads(1, kmer_hash);
            update_stats(0, &counts, stats, &tmp_stats);
            initialise_kmer_counts(stats->n_contaminants, &counts);
        }
        
        prev_full_entry = frw->full_entry;
        
        if (cmd_line->write_progress_file) {
            time(&time_now);
            if (difftime(time_now, time_previous) > cmd_line->progress_delay) {
                kmer_stats_write_progress(stats, cmd_line);
                time_previous = time_now;
            }
        }
    }

    if (cmd_line->write_progress_file) {
        kmer_stats_write_progress(stats, cmd_line);
    }
    
    free_sequence(&frw->seq);
    frw->seq = NULL;
    binary_kmer_free_kmers_set(&windows);
    
    return seq_length;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
long long screen_or_filter_paired_end(CmdLine* cmd_line, KmerFileReaderArgs* fra_1, KmerFileReaderArgs* fra_2, KmerStats* stats)
{
    HashTable* kmer_hash;
	long long seq_length[2];
	int entry_length[2];
	Element *previous_node = NULL;
    boolean keep_reading = true;
    KmerCounts counts[2];
    KmerFileReaderArgs* fra[2];
    KmerFileReaderWrapperArgs* frw[2];
    KmerSlidingWindowSet* windows[2];
    int number_of_files = 1;
    boolean filter_read = false;
    int i;
    time_t time_previous = 0;
    time_t time_now = 0;

    // Get hash table and initialise kmer counts structure
    kmer_hash = fra_1->KmerHash;
    
    // Put file reader args and stats into array
    assert(fra_1 != 0);
    fra[0] = fra_1;
    fra[1] = fra_2;

    // And for second...
    if (fra[1] != 0) {
        number_of_files = 2;
        assert(fra_1->KmerHash == fra_2->KmerHash);
    }
    
    // Allocate...
    for (i=0; i<number_of_files; i++) {
        frw[i] = get_kmer_file_reader_wrapper(kmer_hash->kmer_size, fra[i]);
        
        if (fra[i]->output_filename) {
            frw[i]->output_fp = fopen(fra[i]->output_filename, "w");
            if (!frw[i]->output_fp) {
                printf("Error: can't open output file %s\n", fra[i]->output_filename);
                exit(3);
            } else {
                printf("Opened output %s\n", fra[i]->output_filename);
            }
        }

        if (fra[i]->removed_filename) {
            frw[i]->removed_fp = fopen(fra[i]->removed_filename, "w");
            if (!frw[i]->removed_fp) {
                printf("Error: can't open removed output file %s\n", fra[i]->removed_filename);
                exit(3);
            } else {
                printf("Opened removed %s\n", fra[i]->removed_filename);
            }
        }

        windows[i] = binary_kmer_sliding_window_set_new_from_read_length(kmer_hash->kmer_size, fra[i]->max_read_length);
    }

    for (i=0; i<2; i++) {
        seq_length[i] = 0;
        entry_length[i] = 0;
    }
    
    // Keep reading...
	while (keep_reading)
	{
        KmerStatsReadCounts* stats_both_reads;
        stats_both_reads=calloc(1, sizeof(KmerStatsReadCounts));
        if (!stats_both_reads) {
            printf("Error: Can't get memory for stats\n");
            exit(4);
        }
        
        for (i=0; i<number_of_files; i++) {
            int nkmers;

            initialise_kmer_counts(stats->n_contaminants, &(counts[i]));
            
            // Get next read
            entry_length[i] = file_reader_wrapper(frw[i]);
            if (entry_length[i] == 0) {
                keep_reading = false;
            }

            // Update length read
            seq_length[i] += (long long)entry_length[i];
        
            // Get sliding windows
            nkmers = get_sliding_windows_from_sequence(frw[i]->seq->seq, frw[i]->seq->qual, entry_length[i], fra[i]->quality_cut_off, kmer_hash->kmer_size, windows[i], windows[i]->max_nwindows, windows[i]->max_kmers, false, 0);

            if (frw[i]->full_entry == false) {
                // If we didn't get a full entry then error
                printf("Error: Line length too long.\n");
                return 0;
            }

            if (nkmers == 0) {
                // Bad read
                fra[i]->bad_reads++;
            } else {
                // Load kmers
                kmer_hash_load_sliding_windows(&previous_node, kmer_hash, true, fra[i], kmer_hash->kmer_size, windows[i], i, stats, &(counts[i]));
                hash_table_add_number_of_reads(1, kmer_hash);
                update_stats(i, &(counts[i]), stats, stats_both_reads);
            }
        }
        
        
        if (number_of_files == 2) {
            // Check for not getting both reads
            if (((entry_length[0] == 0) && (entry_length[1] > 0)) ||
                ((entry_length[1] == 0) && (entry_length[0] > 0))) {
                printf("Error: differing number of entries in files.\n");
                return 0;
            }

            // Update both read stats if we got two reads
            if ((entry_length[0] > 0) && (entry_length[1] > 0)) {
                stats->both_reads->number_of_reads++;
                update_stats_for_both(stats, stats_both_reads);
            }
        }
        
        // Keep or filter?
        filter_read=false;
        for (i=0; i<number_of_files; i++) {
            if (counts[i].kmers_loaded > 0) {
                filter_read = true;
            }
        }
        
        // Output reads
        for (i=0; i<number_of_files; i++) {
            if (entry_length[i] > 0) {
                FILE* fp_out = filter_read ? frw[i]->removed_fp:frw[i]->output_fp;
                
                if (fp_out) {
                    char temp_string[frw[i]->seq->length + 1];
                    fprintf(fp_out, "%s%s\n+\n%s\n", frw[i]->seq->id_string, frw[i]->seq->seq, sequence_get_quality_string(frw[i]->seq, temp_string));
                }
            }
        }

        if (cmd_line->write_progress_file) {
            time(&time_now);
            if (difftime(time_now, time_previous) > cmd_line->progress_delay) {
                kmer_stats_write_progress(stats, cmd_line);
                time_previous = time_now;
            }
        }
        
        free(stats_both_reads);
    }
    
    if (cmd_line->write_progress_file) {
        kmer_stats_write_progress(stats, cmd_line);
    }
    
    for (i=0; i<number_of_files; i++) {
        free_sequence(&(frw[i]->seq));
        frw[i]->seq = NULL;
        binary_kmer_free_kmers_set(&(windows[i]));
        fclose(frw[i]->input_fp);
        
        if (frw[i]->output_fp) {
            fclose(frw[i]->output_fp);
        }
        if (frw[i]->removed_fp) {
            fclose(frw[i]->removed_fp);
        }
    }
    
    return seq_length[0] + seq_length[1];
}

/*
void load_fasta_file(char* filename, HashTable* kmer_hash, CmdLine* cmd_line)
{
    char buffer[1024];
    int kmer_size = cmd_line.kmer_size;
    
    fp = fopen(filename, "r");

    
    while (!feof(fp)) {
        if (fgets(buffer, 1024, fp)) {
            if (buffer[0] == '>') {
                printf("Got ID %s\n", buffer);
            } else {
                int o;
                
                for (o=0; o<strlen(buffer) - kmer_size; o++) {
                    char kmer[kmer_size + 1];
                    strncpy(kmer, buffer + o, kmer_size);
                    kmer[kmer_size] = 0;
                    printf("Got kmer %s\n", kmer);
                }
            }
        }
        
    }
}
 */
     
/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
long long load_seq_into_kmers_hash(KmerFileReaderArgs* fra, KmerFileReaderWrapperArgs* fria)
{
    long long *bad_reads = &fra->bad_reads;
    int fastq_ascii_offset = fra->fastq_ascii_offset;
    char quality_cut_off = fra->quality_cut_off;
    int max_read_length = fra->max_read_length;
    HashTable * kmer_hash = fra->KmerHash;
    KmerCounts counts;
	long long seq_length = 0;
	long long seq_count = 0;
    long long contaminated_reads = 0;
	short kmer_size = kmer_hash->kmer_size;
	int entry_length;
    int kmers_loaded;
	boolean prev_full_entry = true;
	Element *previous_node = NULL;
    boolean keep_reading = true;
    
    // Open file
    fria->input_fp =  fopen(fra->input_filename, "r");
    if (fria->input_fp == NULL) {
		fprintf(stderr, "Error: can't open file %s\n", fra->input_filename);
		exit(1);
	}
    
	// Preallocate the memory used to read the sequences
	Sequence *seq = malloc(sizeof(Sequence));
	if (seq == NULL) {
		fputs("Out of memory trying to allocate Sequence\n", stderr);
		exit(1);
	}
    fria->seq = seq;
	alloc_sequence(seq, max_read_length, LINE_MAX, fastq_ascii_offset);
    KmerSlidingWindowSet * windows = binary_kmer_sliding_window_set_new_from_read_length(kmer_size,max_read_length);
    
    // Read file
	while ((entry_length = file_reader_wrapper(fria)) && keep_reading)
	{
        int nkmers;
		
        seq_length += (long long)(entry_length -  (prev_full_entry == false ? kmer_size : 0));
		nkmers = get_sliding_windows_from_sequence(seq->seq, seq->qual, entry_length, quality_cut_off, kmer_size, windows, windows->max_nwindows, windows->max_kmers, false, 0);

		if (nkmers == 0) {
			(*bad_reads)++;
		} else {
            kmer_hash_load_sliding_windows(&previous_node, kmer_hash, prev_full_entry, fra, kmer_size, windows, 0, 0, &counts);
        }
        
        if (fria->full_entry == false) {
            shift_last_kmer_to_start_of_sequence(seq, entry_length, kmer_size);
        } else {
            hash_table_add_number_of_reads(1, kmer_hash);

            if(kmers_loaded > 0){
                contaminated_reads++;
            }
            
            seq_count++;
            kmers_loaded=0;
        }
        
        prev_full_entry = fria->full_entry ;
        
        if (seq_count % 10000) {
            if (hash_table_percentage_occupied(kmer_hash) > fra->maximum_ocupancy) {
                keep_reading = false;
            }
        }
    }
    
    free_sequence(&seq);
    fria->seq = NULL;
    binary_kmer_free_kmers_set(&windows);
    return seq_length;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
/*
void kmer_print_binary_signature(FILE * fp, uint32_t kmer_size, uint32_t num_cols, uint32_t mean_read_len, uint64_t total_seq)
{
    char magic_number[6];
    uint32_t version = BINVERSION;
    uint32_t num_bitfields = NUMBER_OF_BITFIELDS_IN_BINARY_KMER;
    
    magic_number[0]='K';
    magic_number[1]='M';
    magic_number[2]='E';
    magic_number[3]='R';
    magic_number[4]='S';
    magic_number[5]=' ';
    
    fwrite(magic_number,sizeof(char),6,fp);
    fwrite(&version,sizeof(uint32_t),1,fp);
    fwrite(&kmer_size,sizeof(uint32_t),1,fp);
    fwrite(&num_bitfields, sizeof(uint32_t),1,fp);
    fwrite(&num_cols, sizeof(uint32_t), 1, fp);
    fwrite(&mean_read_len, sizeof(uint32_t), 1, fp);
    fwrite(&total_seq, sizeof(uint64_t), 1, fp);
    fwrite(magic_number,sizeof(char),6,fp);
}
*/

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
/*
static void print_node_binary(Element* node, void * args) {
    struct print_binary_args * arguments = (struct print_binary_args * ) args;
    
    if (arguments->condition(node)) {
        arguments->count++;
        db_node_print_binary(arguments->fout, node, arguments->kmer_size);
    }
}
*/

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
/*
void kmer_hash_dump_binary(char *filename, boolean(*condition) (Element* node), HashTable* db_graph)
{
	FILE *fout;
	int mean_read_len = 0;
    long long total_seq = 0;
    struct print_binary_args args;
    void* pass[1];
    
	fout = fopen(filename, "wb");
	if (fout == NULL) {
		fprintf(stderr, "cannot open %s", filename);
		exit(1);
	}
    
	total_seq = hash_table_get_number_of_reads(db_graph);
	kmer_print_binary_signature(fout,db_graph->kmer_size,1, mean_read_len, total_seq);
    
    args.condition = condition;
    args.count = 0;
    args.fout = fout;
    args.kmer_size = db_graph->kmer_size;
    
    pass[0] = (void*)&args;
    
    hash_table_traverse_with_args(&print_node_binary, (void **) &pass, db_graph);
    
	fclose(fout);
    
	printf("%'lld kmers dumped\n", args.count);
}
*/

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
static boolean check_binary_signature_kmers(FILE* fp, uint32_t kmer_size, uint32_t binary_version, uint32_t* number_of_colours_in_binary, uint32_t* mean_read_len, uint64_t* total_seq)
{
    //size_t read;
    //char magic_number[6];
    //boolean ret = false;
    KmerLibraryHeader header;
    
    if (fread(&header, sizeof(KmerLibraryHeader), 1, fp) <= 0) {
        printf("Error: couldn't read kmer library file header\n");
        return false;
    }
    
    if (strncmp(header.header_word, "KONTAMINANT", 11) != 0) {
        printf("Error: bad header in kmer library file\n");
        return false;
    }

    if (strncmp(header.footer_word, "KONTAMINANT", 11) != 0) {
        printf("Error: bad header in kmer library file\n");
        return false;
    }
    
    if (header.version != BINVERSION) {
        printf("Error: incompatible kmer library file version\n");
        return false;
    }
    
    if (header.kmer_size != kmer_size) {
        printf("Error: kmer library file has different kmer size\n");
        return false;
    }
    
    if (header.num_kmers < 1) {
        printf("Error: bad number of kmers in kmer library file\n");
        return false;
    }
    
    
    
    /*
    read = fread(magic_number,sizeof(char),6,fp);
    if (read>0 &&
        magic_number[0]=='K' &&
        magic_number[1]=='M' &&
        magic_number[2]=='E' &&
        magic_number[3]=='R' &&
        magic_number[4]=='S' &&
        magic_number[5]==' ' )
    {
        uint32_t version;
        read = fread(&version,sizeof(uint32_t),1,fp);
        if (read>0 && version==BINVERSION)
        {
            uint32_t file_kmer_size;
            read = fread(&file_kmer_size,sizeof(uint32_t),1,fp);
            if ((read>0) && (file_kmer_size == kmer_size) )
            {
                uint32_t num_bitfields;
                read = fread(&num_bitfields,sizeof(uint32_t),1,fp);
                if ( (read>0) && (num_bitfields==NUMBER_OF_BITFIELDS_IN_BINARY_KMER) )
                {
                    uint32_t num_cols;
                    read = fread(&num_cols,sizeof(uint32_t),1,fp);
                    
                    if ( (read>0) && (num_cols==1)  )
                    {
                        *number_of_colours_in_binary = num_cols;
                        read = fread(mean_read_len,sizeof(uint32_t),1,fp);
                        if (read>0)
                        {
                            read = fread(total_seq,sizeof(uint64_t),1,fp);
                            if (read>0)
                            {
                                magic_number[0]='\0';
                                magic_number[1]='\0';
                                magic_number[2]='\0';
                                magic_number[3]='\0';
                                magic_number[4]='\0';
                                magic_number[5]='\0';
                                read = fread(magic_number,sizeof(char),6,fp);
                                if ((read>0) &&
                                    magic_number[0]=='K' &&
                                    magic_number[1]=='M' &&
                                    magic_number[2]=='E' &&
                                    magic_number[3]=='R' &&
                                    magic_number[4]=='S' &&
                                    magic_number[5]==' ' )
                                {
                                    ret = true;
                                }
                                else
                                {
                                    printf("Binary header is missing the final magic number; we read %s, mean read len %d and total seq %'qd\n", magic_number, *mean_read_len, (long long int)*total_seq);
                                }
                                
                            }
                            else
                            {
                                printf("Binary header does not contain total seq\n");
                            }
                        }
                        else
                        {
                            printf("Binary header does not contain a mean read length\n");
                        }
                        
                    }
                    else
                    {
                        printf("You are loading  binary with %d colours into a graph with %d colours - incompatible\n", num_cols, NUMBER_OF_COLOURS);
                    }
                }
                else
                {
                    printf("Kmer of binary matches the current graph. However this binary was dumped with a different max_kmer size to that of the current graph\n");
                    printf("This binary uses %d bitfields, and current graph uses %d\n", num_bitfields, NUMBER_OF_BITFIELDS_IN_BINARY_KMER);
                }
            }
            else
            {
                printf("You are loading a binary with kmer=%d into a graph with kmer=%d - incompatible\n", file_kmer_size, kmer_size);
            }
        }
        else
        {
            printf("Binary versions do not match.\n");
        }
    }
    else
    {
        printf("Binary does not have magic number in header. Corrupt, or not a KMERS binary\n");
    }
    
    return ret; */
    
    return true;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
boolean read_kmer_from_file(FILE* fp, short kmer_size, Element* node)
{
	BinaryKmer kmer;
	//Edges edges = 0;
	//uint32_t coverage = 0;
	int read;
    
	read = fread(&kmer, sizeof(bitfield_of_64bits), NUMBER_OF_BITFIELDS_IN_BINARY_KMER, fp);
    
	if (read > 0) {
		element_initialise(node, &kmer, kmer_size);
		//read = fread(&coverage, sizeof(uint32_t), 1, fp);
		//if (read == 0) {
		//	puts("Error: couldn't read coverage\n");
		//	exit(1);
		//}
        
		//read = fread(&edges, sizeof(Edges), 1, fp);
		//if (read == 0) {
		//	puts("Error: couldn't read edges\n");
		//	exit(1);
		//}
	} else {
		return false;
	}
    
	return true;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
uint32_t load_kmer_library(char* filename, int n, int k, HashTable* contaminant_hash)
{
	FILE* fp_bin;
    uint32_t num_colours_in_binary;
    uint32_t mean_read_len;
    uint64_t total_seq;
    boolean all_entries_are_unique = false;
	Element node_from_file;
	boolean found;
	long long count = 0;
	BinaryKmer tmp_kmer;
    
    fp_bin = fopen(filename, "rb");
	if (fp_bin == NULL) {
		printf("Error: Cannot open file [%s]\n", filename);
        printf("Error string: %s\n", strerror(errno));
		exit(1);
	}
    
	if ( !check_binary_signature_kmers(fp_bin, k, BINVERSION, &num_colours_in_binary, &mean_read_len, &total_seq)) {
		exit(1);
	}
    
	//Go through all the entries in the binary file
	while (read_kmer_from_file(fp_bin, k, &node_from_file)) {
		count++;
        
		Element *current_node = NULL;
		if (!all_entries_are_unique) {
			current_node =  hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file), k, &tmp_kmer), &found, contaminant_hash);
		} else {
			current_node = hash_table_insert(element_get_key(element_get_kmer(&node_from_file), k, &tmp_kmer), contaminant_hash);
		}
        
        element_set_contaminant_bit(current_node, n);
        
        current_node->coverage[0] = 0;
        current_node->coverage[1] = 0;
	}
    
	fclose(fp_bin);
    
    return count;
}
