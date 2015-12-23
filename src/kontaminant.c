/*----------------------------------------------------------------------*
 * File:    kontaminant.c                                               *
 * Purpose: kmer based screening and filtering                          *
 * Author:  Richard Leggett                                             *
 *          Ricardo Ramirez-Gonzalez                                    *
 *          The Genome Analysis Centre (TGAC), Norwich, UK              *
 *          richard.leggett@tgac.ac.uk    								*
 *----------------------------------------------------------------------*/

/*
   TO DO:
   - binary files include n and b values in header. Scan all files before
     starting to decide on memory usage.
   - Filtering (and screening?) paired end.
   - Read/write .fasta.gzip files.
   - Write stats to plain text file.
 */

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
#include "kmer_build.h"

/*----------------------------------------------------------------------*
 * Constants
 *----------------------------------------------------------------------*/
#define VERSION "2.0.8"

/*----------------------------------------------------------------------*
 * Function:   chomp
 * Purpose:    Remove hidden characters from end of line
 * Parameters: str -> String to change
 * Returns:    None
 *----------------------------------------------------------------------*/
void chomp(char* str)
{
    int i = strlen(str) - 1;
    
    while ((i > 0) && (str[i] < ' ')) {
        str[i--] = 0;
    }    
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
int strcmp_i(char* a, char* b)
{
    int len_a = strlen(a);
    int len_b = strlen(b);
    int length = len_a > len_b ? len_b:len_a;
    int i = 0;
    int diff = 0;
    
    while ((i<length) && (diff == 0)) {
        if (tolower(a[i]) != tolower(b[i])) {
            diff = a[i] - b[i];
        }
        i++;
    }
    
    if (diff == 0) {
        if (len_a > len_b) {
            diff = 1;
        } else if (len_b > len_a) {
            diff = -1;
        }
    }
    
    return diff;
}


/*----------------------------------------------------------------------*
 * Function:   create_hash_table
 * Purpose:    Create hash table for storing PCR duplicates
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
HashTable* create_hash_table(CmdLine* cmdline, int kmer_size)
{
    HashTable* hash;
    long int entries = pow(2.0, (double)cmdline->bucket_bits)*cmdline->bucket_size;
    long int memory = entries*sizeof(Element);

    printf("Creating hash table for kmer storage...\n");
    printf("                n: %d\n", cmdline->bucket_bits);
    printf("                b: %d\n", cmdline->bucket_size);
    printf("          Entries: %ld\n", entries);
    printf("       Entry size: %ld\n", sizeof(Element));
    printf("  Memory required: %ld MB\n\n", memory/(1024*1024));
    
    hash = hash_table_new(cmdline->bucket_bits, cmdline->bucket_size, 25, 1);
    
    if (hash == NULL) {
        printf("Error: No memory for hash table\n");
        exit(101);
    }
    
    hash->kmer_size = kmer_size;
    
    hash_table_print_stats(hash);
    
    return hash;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void load_contamints_from_file(HashTable* contaminant_hash, KmerStats* stats, CmdLine* cmdline)
{
    char filename[MAX_PATH_LENGTH];
    char con[1024];
    FILE* fp = fopen(cmdline->contaminants_file, "r");
    
    if (!fp) {
        printf("Error: can't open file %s\n", cmdline->contaminants_file);
        exit(1);
    }
    
    stats->n_contaminants = 0;
    
    printf("\nContaminants file\n");
    while (!feof(fp)) {
        if (fgets(con, 1024, fp)) {
            chomp(con);
            if (strlen(con) > 1) {
                printf("Loading contaminant %s\n", con);
                sprintf(filename, "%s/%s.fasta.%d.kmers", cmdline->contaminant_dir, con, cmdline->kmer_size);
                printf("      from filename %s\n", filename);
                
                stats->contaminant_ids[stats->n_contaminants] = malloc(strlen(con) + 1);
                
                if (stats->contaminant_ids[stats->n_contaminants]) {
                    strcpy(stats->contaminant_ids[stats->n_contaminants], con);
                } else {
                    printf("Error: can't allocate memory for string!");
                    exit(1);
                }
                
                stats->contaminant_kmers[stats->n_contaminants] = load_kmer_library(filename, stats->n_contaminants, cmdline->kmer_size, contaminant_hash);
                
                stats->n_contaminants++;
            }
        }
    }
    
    fclose(fp);
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void load_contamints(HashTable* contaminant_hash, KmerStats* stats, CmdLine* cmdline)
{
    char* con = strtok(cmdline->contaminants, ",");
    char filename[MAX_PATH_LENGTH];
    
    if (cmdline->contaminants_file != 0) {
        load_contamints_from_file(contaminant_hash, stats, cmdline);
    } else {
        stats->n_contaminants = 0;

        printf("\n");
        while (con != NULL) {
            printf("Loading contaminant %s\n", con);
            sprintf(filename, "%s/%s.fasta.%d.kmers", cmdline->contaminant_dir, con, cmdline->kmer_size);
            printf("      from filename %s\n", filename);
            
            stats->contaminant_ids[stats->n_contaminants] = malloc(strlen(con) + 1);
            
            if (stats->contaminant_ids[stats->n_contaminants]) {
                strcpy(stats->contaminant_ids[stats->n_contaminants], con);
            } else {
                printf("Error: can't allocate memory for string!");
                exit(1);
            }
            
            stats->contaminant_kmers[stats->n_contaminants] = load_kmer_library(filename, stats->n_contaminants, cmdline->kmer_size, contaminant_hash);
            
            stats->n_contaminants++;
            con = strtok(NULL, ",");
        }
    }
}
            
char* get_leafname(char* pathname)
{
    char* ptr = pathname + strlen(pathname) - 1;

    while ((ptr > pathname) && (*ptr != '/')) {
        ptr--;
    }
    
    if (*ptr == '/') {
        ptr++;
    }
        
    return ptr;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void filter_or_screen(char* filename_1, char* filename_2, HashTable* contaminant_hash, KmerStats* kmer_stats, CmdLine* cmdline)
{
    char* filenames[2];
    KmerFileReaderArgs* fra[2];
    int i;
    int n_files = 0;
    
    assert(filename_1 != 0);
    
    filenames[0] = filename_1;
    filenames[1] = filename_2;
    
    for (i=0; i<2; i++) {
        if (filenames[i] != 0) {
            n_files++;
            fra[i] = calloc(1, sizeof(KmerFileReaderArgs));
            if (!fra[i]) {
                printf("Error: Can't get memory for KmerFileReaderArgs!\n");
                exit(2);
            }

            if ((cmdline->output_prefix) && (cmdline->run_type == DO_FILTER)) {
                fra[i]->output_filename = malloc(strlen(filenames[i]) + strlen(cmdline->output_prefix) + 1);
                if (!fra[i]->output_filename) {
                    printf("Error: Can't get memory for output filenames\n");
                    exit(3);
                }
            } else {
                fra[i]->output_filename = 0;
            }

            if ((cmdline->removed_prefix) && (cmdline->run_type == DO_FILTER)) {
                fra[i]->removed_filename = malloc(strlen(filenames[i]) + strlen(cmdline->removed_prefix) + 1);

                if (!fra[i]->removed_filename) {
                    printf("Error: Can't get memory for output filenames\n");
                    exit(3);
                }
            } else {
                fra[i]->removed_filename = 0;
            }
            
            fra[i]->bad_reads = 0;
            fra[i]->colour = 0;
            fra[i]->format = cmdline->format;
            fra[i]->fastq_ascii_offset = 33;
            fra[i]->input_filename = filenames[i];
            fra[i]->quality_cut_off = 0;
            fra[i]->insert = false;
            fra[i]->max_read_length = 1000;
            fra[i]->maximum_ocupancy = 75;
            fra[i]->KmerHash = contaminant_hash;
        
            if (fra[i]->output_filename) {
                sprintf(fra[i]->output_filename, "%s%s", cmdline->output_prefix, get_leafname(filenames[i]));
            }
            if (fra[i]->removed_filename) {
                sprintf(fra[i]->removed_filename, "%s%s", cmdline->removed_prefix, get_leafname(filenames[i]));
            }
        } else {
            fra[i] = 0;
        }
    }
        
    printf("\n");
    if (n_files == 1) {
        printf("  Source file %s\n", fra[0]->input_filename);
        if (cmdline->output_prefix) {
            printf("Filtered file %s\n", fra[0]->output_filename);
        }
        if (cmdline->removed_prefix) {
            printf("Removed reads %s\n", fra[0]->removed_filename);
        }
    } else {
        printf("  Source pair %s\n          and %s\n", fra[0]->input_filename, fra[1]->input_filename);
        if ((cmdline->output_prefix) && (cmdline->run_type == DO_FILTER)) {
            printf("Filtered pair %s\n          and %s\n", fra[0]->output_filename, fra[1]->output_filename);
        }
        if ((cmdline->removed_prefix) && (cmdline->run_type == DO_FILTER)) {
            printf("Removed reads %s\n          and %s\n", fra[0]->removed_filename, fra[1]->removed_filename);
        }
    }
    printf("\n");
    
    kmer_stats->number_of_files = n_files;

    if (cmdline->format == FASTA) {
        if (n_files == 1) {
            screen_kmers_from_file(fra[0], cmdline, kmer_stats);
        } else {
            printf("Error: paired FASTA files not supported.\n");
            exit(3);
        }
    } else if (cmdline->format == FASTQ) {
        screen_or_filter_paired_end(cmdline, fra[0], fra[1], kmer_stats);
    } else {
        printf("Error: Format not supported.\n");
    }
    
    hash_table_print_stats(contaminant_hash);
    
    //return loaded_kmers;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void process_files(HashTable* contaminant_hash, KmerStats* kmer_stats, CmdLine* cmdline)
{
    if (cmdline->file_of_files != 0) {
        printf("Need to implement file of file code!\n");
        exit(1);
        //FILE* fp = fopen(cmdline->input_filename_one, "r");
        //if (fp) {
        //    while (!feof(fp)) {
        //        char filename[MAX_PATH_LENGTH];
        //        if (fgets(filename, MAX_PATH_LENGTH, fp)) {
        //            if (strlen(filename) > 1) {
        //                filter_or_screen(filename, contaminant_hash, kmer_stats, cmdline);
        //            }
        //        }
        //    }
        //    fclose(fp);
        //}
    } else {
        filter_or_screen(cmdline->input_filename_one, cmdline->input_filename_two, contaminant_hash, kmer_stats, cmdline);
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void initialise_output_files(CmdLine* cmdline, KmerStats* stats)
{
    if (cmdline->read_summary_file != 0) {
        FILE* fp = fopen(cmdline->read_summary_file, "w");
        int i;
        
        if (fp) {
            printf("\nOpened %s\n", cmdline->read_summary_file);
            fprintf(fp, "ID\tnCons\tnKs");
            for (i=0; i<stats->n_contaminants; i++) {
                fprintf(fp, "\t%s", stats->contaminant_ids[i]);
            }
            fprintf(fp, "\n");
            fclose(fp);
        } else {
            printf("Error: can't open %s\n", cmdline->read_summary_file);
            cmdline->read_summary_file = 0;
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:   main
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    time_t start, end;
    double seconds;
    HashTable* contaminant_hash = NULL;
    CmdLine cmdline;
    KmerStats kmer_stats;
    
    //Element e;
    //int i;
    
    //for (i=0; i<MAX_CONTAMINANTS; i++) {
    //    element_set_contaminant_bit(&e, i);
    //}
    //exit(0);
    
    time(&start);
    
    printf("\nkONTAMINANT v%s\n\n", VERSION);
    
    printf("Max contaminants: %d\n\n", MAX_CONTAMINANTS);
    
    initialise_cmdline(&cmdline);
    parse_command_line(argc, argv, &cmdline);
    kmer_stats_initialise(&kmer_stats, &cmdline);
    contaminant_hash = create_hash_table(&cmdline, cmdline.kmer_size);

    if ((cmdline.run_type == DO_INDEX)) {
        // Index file
        load_reads_into_table(&cmdline, contaminant_hash);
        printf("\n");
        hash_table_print_stats(contaminant_hash);
        dump_kmer_hash(&cmdline, contaminant_hash);
    } else if ((cmdline.run_type == DO_SCREEN) || (cmdline.run_type == DO_FILTER)) {
        load_contamints(contaminant_hash, &kmer_stats, &cmdline);
        initialise_output_files(&cmdline, &kmer_stats);
        printf("\n");
        hash_table_print_stats(contaminant_hash);
        kmer_stats_compare_contaminant_kmers(contaminant_hash, &kmer_stats, &cmdline);

        time(&end);
        seconds = difftime(end, start);
        printf("\nProcessing read files (after %.0f seconds)...\n", seconds);
        
        process_files(contaminant_hash, &kmer_stats, &cmdline);

        time(&end);
        seconds = difftime(end, start);
        printf("\nCalculating stats (after %.0f seconds)...\n", seconds);
        
        kmer_stats_calculate(&kmer_stats);
        kmer_stats_report_to_screen(&kmer_stats, &cmdline);
    }
    
    time(&end);
    seconds = difftime(end, start);

    printf("\nDone. Completed in %.0f seconds.\n\n", seconds);
    
        

    return 0;
}
