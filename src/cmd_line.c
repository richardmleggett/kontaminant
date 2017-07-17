/*----------------------------------------------------------------------*
 * File:    cmd_line.c                                                   *
 * Purpose: Handle command line                                         *
 * Author:  Richard Leggett                                             *
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
#include <pthread.h>
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
void initialise_cmdline(CmdLine* c)
{
    c->bucket_bits = 20;
    c->bucket_size = 100;
    c->run_type = 0;
    c->kmer_size = 21;
    c->keep_contaminated_reads = false;
    c->input_filename_one = 0;
    c->input_filename_two = 0;
    c->file_of_files = 0;
    c->progress_dir = 0;
    c->contaminants = 0;
    c->contaminants_file = 0;
    c->contaminant_dir = malloc(MAX_PATH_LENGTH);
    strcpy(c->contaminant_dir, "/Users/leggettr/Documents/Projects/kscreen/contaminants");
    c->output_prefix = 0;
    c->subsample_ratio = 1;
    c->quality_score_offset = 33;
    c->quality_score_threshold = 0;
    c->format = FASTQ;
    c->max_read_length = 200000;
    c->kmer_threshold_overall = 10;
    c->kmer_threshold_read = 1;
    c->output_prefix = malloc(8);
    strcpy(c->output_prefix, "kout_");
    c->removed_prefix = 0;
    c->write_progress_file = false;
    c->progress_delay = 60;
    c->read_summary_file = 0;
    c->filter_unique = false;
    c->numthreads = 1;
    c->ratio = 1.0;
}

/*----------------------------------------------------------------------*
 * Function:   usage
 * Purpose:    Report program usage.
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void usage(void)
{
    printf("k-mer based screening and filtering of reads\n" \
           "\nSyntax: kontaminant <-s|-f|-i> [options]\n" \
           "\nWhere:\n" \
           "    [-s | --screen] invokes screening.\n" \
           "    [-f | --filter] invokes filtering.\n" \
           "    [-i | --index] indexes a reference.\n" \
           "Kmer options:\n" \
           "    [-k | --kmer_size] Kmer size (default 21).\n" \
           "    [-t | --threshold] Kmer threshold for both reads (default 10).\n" \
           "    [-l | --readthreshold] Kmer threshold for individual reads (default 1).\n" \
           "    [-y | --subsample] Ratio of reads to sample >0 <=1 (default 1).\n" \
           "    [-u | --unique] Count only unique kmers (default off).\n" \
           "Input options:\n" \
           "    [-1 | --input_one] Input R1 file (or reference FASTA for indexing).\n" \
           "    [-2 | --input_two] Input R2 file.\n" \
           "    [-g | --file_format] Input file format FASTA or FASTQ (default FASTQ).\n" \
           "    [-z | --file_of_files] Input file of files (for batch processing - instead of -1 and -2).\n" \
           "Output options:\n" \
           "    [-j | --read_summary] Read summary file.\n" \
           "    [-o | --output_prefix] Output prefix (default: 'kout_').\n" \
           "    [-p | --progress] Name of directory for streaming progress page.\n"
           "    [-r | --removed_prefix] Removed reads prefix (filtering only).\n" \
           "    [-x | --keep_contaminated_reads] Save contaminated reads into separate file.\n" \
           "Contaminant options:\n" \
           "    [-d | --contaminant_dir] Contaminant library directory.\n" \
           "    [-c | --contaminants] List of contaminants to screen/filter, OR\n" \
           "    [-e | --contaminants_file] Filename of file containing list of contaminants to screen/filer.\n" \
           "Memory options:\n" \
           "    [-b | --mem_width] Size of hash table buckets (default 100).\n" \
           "    [-n | --mem_height] Number of buckets in hash table in bits (default 20, this is a power of 2, ie 2^mem_height).\n" \
           "\nComments/suggestions to richard.leggett@tgac.ac.uk\n" \
           "\n");
}

/*----------------------------------------------------------------------*
 * Function:   parse_command_line
 * Purpose:    Parse command line options
 * Parameters: argc = number of arguments
 *             argv -> array of arguments
 * Returns:    None
 *----------------------------------------------------------------------*/
void parse_command_line(int argc, char* argv[], CmdLine* c)
{
    static struct option long_options[] = {
        {"input_one", required_argument, NULL, '1'},
        {"input_two", required_argument, NULL, '2'},
        {"mem_width", required_argument, NULL, 'b'},
        {"contaminants", required_argument, NULL, 'c'},
        {"contaminant_dir", required_argument, NULL, 'd'},
        {"contaminants_file", required_argument, NULL, 'e'},
        {"filter", no_argument, NULL, 'f'},
        {"file_format", required_argument, NULL, 'g'},
        {"help", no_argument, NULL, 'h'},
        {"index", no_argument, NULL, 'i'},
        {"read_summary", required_argument, NULL, 'j'},
        {"kmer_size", required_argument, NULL, 'k'},
        {"readthreshold", required_argument, NULL, 'l'},
        {"mem_height", required_argument, NULL, 'n'},
        {"mem_height", required_argument, NULL, 'n'},
        {"numthreads", required_argument, NULL, 'N'},
        {"progress", required_argument, NULL, 'p'},
        {"removed_prefix", required_argument, NULL, 'r'},
        {"ratio", required_argument, NULL, 'R'},
        {"screen", no_argument, NULL, 's'},
        {"threshold", required_argument, NULL, 't'},
        {"unique", no_argument, NULL, 'u'},
        {"progress_interval", required_argument, NULL, 'w'},
        {"keep_contaminated_reads", no_argument, NULL, 'x'},
        {"subsample", required_argument, NULL, 'y'},
        {"file_of_files", required_argument, NULL, 'z'},
        {0, 0, 0, 0}
    };
    int opt;
    int longopt_index;
    
    if (argc == 1) {
        usage();
        exit(0);
    }
    
    while ((opt = getopt_long(argc, argv, "1:2:b:c:d:e:fg:hij:k:l:n:N:o:p:r:R:st:uw:xy:z:", long_options, &longopt_index)) > 0)
    {
        switch(opt) {
            case '1':
                if (optarg==NULL) {
                    printf("Error: [-a | --input] option requires an argument.\n");
                    exit(1);
                }
                c->input_filename_one = malloc(strlen(optarg) + 1);
                if (c->input_filename_one) {
                    strcpy(c->input_filename_one, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case '2':
                if (optarg==NULL) {
                    printf("Error: [-b | --input] option requires an argument.\n");
                    exit(1);
                }
                c->input_filename_two = malloc(strlen(optarg) + 1);
                if (c->input_filename_two) {
                    strcpy(c->input_filename_two, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case 'b':
                if (optarg == NULL) {
                    printf("[-b | --mem_width] option requires int argument [hash table bucket size]");
                    exit(1);
                }
                c->bucket_size = atoi(optarg);
                if (c->bucket_size == 0 || c->bucket_size > 100000) {
                    printf("[-b | --mem_width] option requires 'short' argument bigger than 0");
                    exit(1);
                }
                break;
            case 'c':
                if (optarg==NULL) {
                    printf("Error: [-c | --contaminants] option requires an argument.\n");
                    exit(1);
                }
                c->contaminants = malloc(strlen(optarg) + 1);
                if (c->contaminants) {
                    strcpy(c->contaminants, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case 'd':
                if (optarg==NULL) {
                    printf("Error: [-d | --contaminant_dir] option requires an argument.\n");
                    exit(1);
                }
                c->contaminant_dir = malloc(strlen(optarg) + 1);
                if (c->contaminant_dir) {
                    strcpy(c->contaminant_dir, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case 'e':
                if (optarg==NULL) {
                    printf("Error: [-e | --contaminants_file] option requires an argument.\n");
                    exit(1);
                }
                c->contaminants_file = malloc(strlen(optarg) + 1);
                if (c->contaminants_file) {
                    strcpy(c->contaminants_file, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case 'f':
                if (c->run_type == 0) {
                    c->run_type = DO_FILTER;
                } else {
                    printf("Error: You must specify either screening, filtering or indexing.\n");
                    exit(1);
                }
                break;
            case 'g':
                if (optarg==NULL) {
                    printf("Error: [-g | --file_format] option requires an argument FASTA or FASTQ.\n");
                    exit(1);
                }
                if (strcmp(optarg, "FASTA") == 0) {
                    c->format = FASTA;
                } else if (strcmp(optarg, "FASTQ") == 0) {
                    c->format = FASTQ;
                } else {
                    printf("Error: [-g | --file_format] option requires an argument FASTA or FASTQ.\n");
                    exit(1);
                }
                break;
            case 'h':
                usage();
                exit(0);
                break;
            case 'i':
                if (c->run_type == 0) {
                    c->run_type = DO_INDEX;
                } else {
                    printf("Error: You must specify either screening, filtering or indexing.\n");
                    exit(1);
                }
                break;
            case 'j':
                if (optarg==NULL) {
                    printf("Error: [-j | --read_summary] option requires an argument.\n");
                    exit(1);
                }
                c->read_summary_file = malloc(strlen(optarg) + 1);
                if (c->read_summary_file) {
                    strcpy(c->read_summary_file, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case 'k':
                if (optarg==NULL) {
                    printf("Error: [-k | --kmer_size] option requires an argument.\n");
                    exit(1);
                }
                c->kmer_size = atoi(optarg);
                break;
            case 'l':
                if (optarg==NULL) {
                    printf("Error: [-l | --readthreshold] option requires an argument.\n");
                    exit(1);
                }
                c->kmer_threshold_read = atoi(optarg);
                break;
            case 'n':
                if (optarg == NULL) {
                    printf("[-n | --mem_height] option requires int argument [hash table number of buckets in bits]");
                    exit(1);
                }
                c->bucket_bits = atoi(optarg);
                break;
            case 'N':
                if (optarg == NULL) {
                    printf("[-N | --numthreads] option requires int argument");
                    exit(1);
                }
                c->numthreads = atoi(optarg);
                if (c->numthreads > 32) {
                    c->numthreads = 32;
                }
                break;
            case 'o':
                if (optarg==NULL) {
                    printf("Error: [-o | --output_prefix] option requires an argument.\n");
                    exit(1);
                }
                c->output_prefix = malloc(strlen(optarg) + 1);
                if (c->output_prefix) {
                    strcpy(c->output_prefix, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case 'p':
                if (optarg==NULL) {
                    printf("Error: [-p | --progress_file] option requires an argument.\n");
                    exit(1);
                }
                c->progress_dir = malloc(strlen(optarg) + 1);
                if (c->progress_dir) {
                    strcpy(c->progress_dir, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                c->write_progress_file = true;
                break;
            case 'r':
                if (optarg==NULL) {
                    printf("Error: [-r | --removed_prefix] option requires an argument.\n");
                    exit(1);
                }
                c->removed_prefix = malloc(strlen(optarg) + 1);
                if (c->removed_prefix) {
                    strcpy(c->removed_prefix, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            case 'R':
                if (optarg==NULL) {
                    printf("Error: [-R | --ratio] option requires an argument.\n");
                    exit(1);
                }
                c->ratio = atof(optarg);
                break;
            case 's':
                if (c->run_type == 0) {
                    c->run_type = DO_SCREEN;
                } else {
                    printf("Error: You must specify either screening, filtering or indexing.\n");
                    exit(1);
                }
                break;
            case 't':
                if (optarg==NULL) {
                    printf("Error: [-t | --threshold] option requires an argument.\n");
                    exit(1);
                }
                c->kmer_threshold_overall = atoi(optarg);
                break;
            case 'u':
                c->filter_unique = true;
                break;
            case 'w':
                if (optarg == NULL) {
                    printf("[-u | --progress_interval] option requires int argument [delay in seconds]");
                    exit(1);
                }
                c->progress_delay = atoi(optarg);
                break;
            case 'x':
                c->keep_contaminated_reads = true;
                break;
            case 'y':
                if (optarg==NULL) {
                    printf("Error: [-k | --kmer_size] option requires an argument.\n");
                    exit(1);
                }
                c->subsample_ratio = atof(optarg);
                break;
            case 'z':
                if (optarg==NULL) {
                    printf("Error: [-z | --file_of_files] option requires an argument.\n");
                    exit(1);
                }
                c->file_of_files = malloc(strlen(optarg) + 1);
                if (c->file_of_files) {
                    strcpy(c->file_of_files, optarg);
                } else {
                    printf("Error: can't allocate memory for string.\n");
                    exit(1);
                }
                break;
            default:
                printf("Error: Unknown option %c\n", opt);
                exit(1);
                break;
        }
    }
    
    if (c->file_of_files == 0) {
        if (c->input_filename_one == 0) {
            printf("Error: you must specify an input filename.\n");
            exit(1);
        }
    } else {
        if ((c->input_filename_one !=0) || (c->input_filename_two !=0)) {
            printf("Error: if you specify a file of files, you must not specify other input filenames.\n");
        }
    }
    
    if (c->contaminant_dir == 0) {
        printf("Error: you must specify a contaminant dir\n");
        exit(1);
    }
    
    if ((c->run_type == DO_SCREEN) || (c->run_type == DO_FILTER)) {
        if ((c->contaminants == 0) && (c->contaminants_file == 0)) {
            printf("Error: you must specify a contaminant list\n");
            exit(1);
        }
    }
    
    if (c->output_prefix == 0) {
        printf("Error: you must specify an output prefix\n");
        exit(1);
    }
    
    if ((c->kmer_threshold_read * 2) > c->kmer_threshold_overall) {
        printf("NOTE: by specifying a read threshold of %d, you require at least %d kmers over both reads, which has the effect of increasing your overall threshold past the specified value (%d). This isn't a problem, as long as you are aware.\n\n", c->kmer_threshold_read, c->kmer_threshold_read*2, c->kmer_threshold_overall);
    }
}
