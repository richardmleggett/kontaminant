/*----------------------------------------------------------------------*
 * File:    kmer_stats.c                                                *
 * Purpose: Handle calculation and display of stats                     *
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
void kmer_stats_read_counts_initialise(KmerStatsReadCounts *r)
{
    int i;
    
    r->number_of_reads = 0;
    r->k1_contaminated_reads = 0;
    r->kn_contaminated_reads = 0;
    
    for (i=0; i<MAX_CONTAMINANTS; i++) {
        r->contaminant_kmers_seen[i] = 0;
        r->k1_contaminated_reads_by_contaminant[i] = 0;
        r->k1_unique_contaminated_reads_by_contaminant[i] = 0;
    }
    
    for (i=0; i<MAX_READ_LENGTH; i++) {
        r->contaminated_kmers_per_read[i] = 0;
    }
    
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_both_reads_initialise(KmerStatsBothReads *r)
{
    int i;
    
    r->number_of_reads = 0;
    r->both_k1_contaminated_reads = 0;
    r->either_k1_contaminated_reads = 0;
    r->both_kn_contaminated_reads = 0;
    r->either_kn_contaminated_reads = 0;
    
    for (i=0; i<MAX_CONTAMINANTS; i++) {
        r->k1_contaminated_reads_by_contaminant[i] = 0;
        r->k1_unique_contaminated_reads_by_contaminant[i] = 0;
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_initialise(KmerStats* stats, CmdLine* cmd_line)
{
    int i;
    int j;
    int r;

    stats->n_contaminants = 0;
    stats->number_of_files = 0;
    stats->kmer_threshold = cmd_line->kmer_threshold;
    
    for (i=0; i<MAX_CONTAMINANTS; i++) {
        stats->contaminant_kmers[i] = 0;
    }

    for (i=0; i<MAX_CONTAMINANTS; i++) {
        stats->unique_kmers[i] = 0;
        for (j=0; j<MAX_CONTAMINANTS; j++) {
            stats->kmers_in_common[i][j] = 0;
        }
    }
    
    for (r=0; r<2; r++) {
        stats->read[r] = calloc(1, sizeof(KmerStatsReadCounts));
        if (stats->read[r] == 0) {
            printf("Error: Can't allocate space for KmerStatsReadCounts!\n");
            exit(2);
        }
        kmer_stats_read_counts_initialise(stats->read[r]);
    }
    
    stats->both_reads = calloc(1, sizeof(KmerStatsBothReads));
    if (stats->both_reads == 0) {
        printf("Error: Can't allocate space for KmerStatsBothReads!\n");
        exit(2);
    }
    kmer_stats_both_reads_initialise(stats->both_reads);
}

/*----------------------------------------------------------------------*
 * Function:   update_stats
 * Purpose:    Update overall read stats from a KmerCounts read structure
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void update_stats(int r, KmerCounts* counts, KmerStats* stats, KmerStatsReadCounts* both_stats)
{
    int i;
    int largest_contaminant = 0;
    int largest_kmers = 0;
    int unique_largest_contaminant = 0;
    int unique_largest_kmers = 0;
    
    stats->read[r]->number_of_reads++;
    
    // We allow up to the maximum read length. The last element is the cumulative of the reads/contigs that have more than the space we have allocated
    stats->read[r]->contaminated_kmers_per_read[counts->kmers_loaded < MAX_READ_LENGTH? counts->kmers_loaded:MAX_READ_LENGTH]++;
    
    if (counts->kmers_loaded > 0) {
        for (i=0; i<stats->n_contaminants; i++) {
            if (counts->kmers_from_contaminant[i] > 0) {
                if (counts->kmers_from_contaminant[i] > largest_kmers) {
                    largest_kmers = counts->kmers_from_contaminant[i];
                    largest_contaminant = i;
                }

                if (counts->unique_kmers_from_contaminant[i] > unique_largest_kmers) {
                    unique_largest_kmers = counts->unique_kmers_from_contaminant[i];
                    unique_largest_contaminant = i;
                }
                
                stats->read[r]->k1_contaminated_reads_by_contaminant[i]++;
                both_stats->k1_contaminated_reads_by_contaminant[i]++;
                if (counts->contaminants_detected == 1) {
                    stats->read[r]->k1_unique_contaminated_reads_by_contaminant[i]++;
                    both_stats->k1_unique_contaminated_reads_by_contaminant[i]++;
                }
            }
        }
        
        stats->read[r]->k1_contaminated_reads++;
        both_stats->k1_contaminated_reads++;
    }
    
    if (largest_kmers == 0) {
        stats->read[r]->reads_unclassified++;
        counts->assigned_contaminant = -1;
    } else {
        stats->read[r]->reads_with_highest_contaminant[largest_contaminant]++;
        counts->assigned_contaminant = largest_contaminant;
    }
    
    if (unique_largest_kmers == 0) {
        counts->unique_assigned_contaminant = -1;
    } else {
        counts->unique_assigned_contaminant = unique_largest_contaminant;
    }
    
    if (counts->kmers_loaded > stats->kmer_threshold) {
        for (i=0; i<stats->n_contaminants; i++) {
            if (counts->kmers_from_contaminant[i] > 0) {
                stats->read[r]->kn_contaminated_reads_by_contaminant[i]++;
                both_stats->kn_contaminated_reads_by_contaminant[i]++;
                if (counts->contaminants_detected == 1) {
                    stats->read[r]->kn_unique_contaminated_reads_by_contaminant[i]++;
                    both_stats->kn_unique_contaminated_reads_by_contaminant[i]++;
                }
            }
        }
        
        stats->read[r]->kn_contaminated_reads++;
        both_stats->kn_contaminated_reads++;
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void update_stats_for_both(KmerStats* stats, KmerStatsReadCounts* both_stats)
{
    int i;
        
    for (i=0; i<stats->n_contaminants; i++) {
        if (both_stats->k1_contaminated_reads_by_contaminant[i] == 2) {
            stats->both_reads->k1_contaminated_reads_by_contaminant[i]++;
        }

        if (both_stats->k1_unique_contaminated_reads_by_contaminant[i] == 2) {
            stats->both_reads->k1_unique_contaminated_reads_by_contaminant[i]++;
        }

        if (both_stats->kn_contaminated_reads_by_contaminant[i] == 2) {
            stats->both_reads->kn_contaminated_reads_by_contaminant[i]++;
        }
        
        if (both_stats->kn_unique_contaminated_reads_by_contaminant[i] == 2) {
            stats->both_reads->kn_unique_contaminated_reads_by_contaminant[i]++;
        }
    }

    if (both_stats->k1_contaminated_reads == 2) {
        stats->both_reads->both_k1_contaminated_reads++;
    } else if (both_stats->k1_contaminated_reads == 1) {
        stats->both_reads->either_k1_contaminated_reads++;
    }

    if (both_stats->kn_contaminated_reads == 2) {
        stats->both_reads->both_kn_contaminated_reads++;
    } else if (both_stats->kn_contaminated_reads == 1) {
        stats->both_reads->either_kn_contaminated_reads++;
    }        

}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_calculate_read(KmerStats* stats, KmerStatsReadCounts* read)
{
    int i;
    
    read->k1_contaminaned_reads_pc = (100.0 * (double)read->k1_contaminated_reads) / (double)read->number_of_reads;
    read->kn_contaminaned_reads_pc = (100.0 * (double)read->kn_contaminated_reads) / (double)read->number_of_reads;
    
    for (i=0; i<stats->n_contaminants; i++) {
        read->k1_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)read->k1_contaminated_reads_by_contaminant[i]) / (double)read->number_of_reads;
        read->k1_unique_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)read->k1_unique_contaminated_reads_by_contaminant[i]) / (double)read->number_of_reads;
        read->kn_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)read->kn_contaminated_reads_by_contaminant[i]) / (double)read->number_of_reads;
        read->kn_unique_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)read->kn_unique_contaminated_reads_by_contaminant[i]) / (double)read->number_of_reads;
        read->contaminant_kmers_seen_pc[i] = (100.0 * (double)read->contaminant_kmers_seen[i]) / (double)stats->contaminant_kmers[i];
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_calculate_both(KmerStats* stats)
{
    int i;
    
    stats->both_reads->both_k1_contaminaned_reads_pc = (100.0 * (double)stats->both_reads->both_k1_contaminated_reads) / (double)stats->both_reads->number_of_reads;
    stats->both_reads->both_kn_contaminaned_reads_pc = (100.0 * (double)stats->both_reads->both_kn_contaminated_reads) / (double)stats->both_reads->number_of_reads;
    stats->both_reads->either_k1_contaminaned_reads_pc = (100.0 * (double)stats->both_reads->either_k1_contaminated_reads) / (double)stats->both_reads->number_of_reads;
    stats->both_reads->either_kn_contaminaned_reads_pc = (100.0 * (double)stats->both_reads->either_kn_contaminated_reads) / (double)stats->both_reads->number_of_reads;
    
    for (i=0; i<stats->n_contaminants; i++) {
        stats->both_reads->k1_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)stats->both_reads->k1_contaminated_reads_by_contaminant[i]) / (double)stats->both_reads->number_of_reads;
        stats->both_reads->k1_unique_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)stats->both_reads->k1_unique_contaminated_reads_by_contaminant[i]) / (double)stats->both_reads->number_of_reads;
        stats->both_reads->kn_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)stats->both_reads->kn_contaminated_reads_by_contaminant[i]) / (double)stats->both_reads->number_of_reads;
        stats->both_reads->kn_unique_contaminated_reads_by_contaminant_pc[i] = (100.0 * (double)stats->both_reads->kn_unique_contaminated_reads_by_contaminant[i]) / (double)stats->both_reads->number_of_reads;
        stats->both_reads->contaminant_kmers_seen_pc[i] = (100.0 * (double)stats->both_reads->contaminant_kmers_seen[i]) / (double)stats->contaminant_kmers[i];
    }
}


/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_calculate(KmerStats* stats)
{
    kmer_stats_calculate_read(stats, stats->read[0]);
    kmer_stats_calculate_read(stats, stats->read[1]);
    kmer_stats_calculate_both(stats);
 }

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_report_read_stats(KmerStats* stats, KmerStatsReadCounts* read)
{
    int i;
    char str[1024];
    
    printf("Overall statistics\n\n");
    printf("%64s: %d\n", "Number of reads", read->number_of_reads);
    printf("%64s: %d\n", "Number of reads with 1+ kmer contamination", read->k1_contaminated_reads);
    printf("%64s: %.2f\n", "%% of reads with 1+ kmer contamination", read->k1_contaminaned_reads_pc);
    sprintf(str, "Number of reads with %d+ kmer contamination", stats->kmer_threshold);
    printf("%64s: %d\n", str, read->kn_contaminated_reads);
    sprintf(str, "%% of reads with %d+ kmer contamination", stats->kmer_threshold);
    printf("%64s: %.2f\n", str, read->kn_contaminaned_reads_pc);
    
    printf("\nPer-contaminant statistics\n\n");
    
    printf("%-30s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n", "Contaminant", "nKmers", "kFound", "%%kFound", "ReadsW1k", "%%ReadsW1k", "UniqW1k", "%%UniqW1k", "ReadsWnk", "%%ReadsWnk", "UniqWnk", "%%UniqWnk");
           
    for (i=0; i<stats->n_contaminants; i++) {
        printf("%-30s %-10d %-10d %-10.2f %-10d %-10.2f %-10d %-10.2f %-10d %-10.2f %-10d %-10.2f\n",
               stats->contaminant_ids[i],
               stats->contaminant_kmers[i],
               read->contaminant_kmers_seen[i],
               read->contaminant_kmers_seen_pc[i],
               read->k1_contaminated_reads_by_contaminant[i],
               read->k1_contaminated_reads_by_contaminant_pc[i],
               read->k1_unique_contaminated_reads_by_contaminant[i],
               read->k1_unique_contaminated_reads_by_contaminant_pc[i],
               read->kn_contaminated_reads_by_contaminant[i],
               read->kn_contaminated_reads_by_contaminant_pc[i],
               read->kn_unique_contaminated_reads_by_contaminant[i],
               read->kn_unique_contaminated_reads_by_contaminant_pc[i]);
    }
    
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_report_both_stats(KmerStats* stats)
{
    int i;
    char str[1024];
    
    printf("Overall statistics\n\n");
    printf("%64s: %d\n", "Number of pairs", stats->both_reads->number_of_reads);
    printf("%64s: %d\n", "Number of pairs with 1+ kmer contamination in only R1 or R2", stats->both_reads->either_k1_contaminated_reads);
    printf("%64s: %.2f\n", "%% of pairs with 1+ kmer contamination in only R1 or R2", stats->both_reads->either_k1_contaminaned_reads_pc);
    sprintf(str, "Number of pairs with %d+ kmer contamination in only R1 or R2", stats->kmer_threshold);
    printf("%64s: %d\n", str, stats->both_reads->either_kn_contaminated_reads);
    sprintf(str, "%% of pairs with %d+ kmer contamination in only R1 or R2", stats->kmer_threshold);
    printf("%64s: %.2f\n", str, stats->both_reads->either_kn_contaminaned_reads_pc);
    printf("%64s: %d\n", "Number of pairs with 1+ kmer contamination in both R1 and R2", stats->both_reads->both_k1_contaminated_reads);
    printf("%64s: %.2f\n", "%% of pairs with 1+ kmer contamination in both R1 and R2", stats->both_reads->both_k1_contaminaned_reads_pc);
    sprintf(str, "Number of pairs with %d+ kmer contamination in both R1 and R2", stats->kmer_threshold);
    printf("%64s: %d\n", str, stats->both_reads->both_kn_contaminated_reads);
    sprintf(str, "%% of pairs with %d+ kmer contamination in both R1 and R2", stats->kmer_threshold);
    printf("%64s: %.2f\n", str, stats->both_reads->both_kn_contaminaned_reads_pc);
    
    printf("\nPer-contaminant statistics\n\n");

    printf("%-30s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n", "Contaminant", "nKmers", "kFound", "%%kFound", "ReadsW1k", "%%ReadsW1k", "UniqW1k", "%%UniqW1k", "ReadsWnk", "%%ReadsWnk", "UniqWnk", "%%UniqWnk");
    
    for (i=0; i<stats->n_contaminants; i++) {
        printf("%-30s %-10d %-10d %-10.2f %-10d %-10.2f %-10d %-10.2f %-10d %-10.2f %-10d %-10.2f\n",
               stats->contaminant_ids[i],
               stats->contaminant_kmers[i],
               stats->both_reads->contaminant_kmers_seen[i],
               stats->both_reads->contaminant_kmers_seen_pc[i],
               stats->both_reads->k1_contaminated_reads_by_contaminant[i],
               stats->both_reads->k1_contaminated_reads_by_contaminant_pc[i],
               stats->both_reads->k1_unique_contaminated_reads_by_contaminant[i],
               stats->both_reads->k1_unique_contaminated_reads_by_contaminant_pc[i],
               stats->both_reads->kn_contaminated_reads_by_contaminant[i],
               stats->both_reads->kn_contaminated_reads_by_contaminant_pc[i],
               stats->both_reads->kn_unique_contaminated_reads_by_contaminant[i],
               stats->both_reads->kn_unique_contaminated_reads_by_contaminant_pc[i]);
    }
    
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_report_to_screen(KmerStats* stats)
{
    int r;
    
    for (r=0; r<stats->number_of_files; r++) {
        printf("\n========== Statistics for Read %d ===========\n\n", r+1);
        kmer_stats_report_read_stats(stats, stats->read[r]);
    }

    if (stats->number_of_files == 2) {
        printf("\n========== Statistics for both reads ===========\n\n");
        kmer_stats_report_both_stats(stats);
    }
    
    printf("\n========== Key ==========\n\n");
    printf("nKmers    - Number of kmers in contaminant reference\n");
    printf("kFound    - Number of unique contaminant kmers found in reads\n");
    printf("%%kFound   - Percentage of contaminant kmers found in reads\n");
    printf("ReadsW1k  - Reads containing 1 or more kmer from the contaminant\n");
    printf("%%ReadsW1k - Percentage of reads containing 1 or more kmer from the contaminant\n");
    printf("UniqW1k   - Reads containing 1 or more kmer from the contaminant and not any other\n");
    printf("%%UniqW1k  - Percentage of reads containing 1 or more kmer from the contaminant and not any other\n");
    printf("ReadsWnk  - Reads containing n or more kmer from the contaminant (n=%d)\n", stats->kmer_threshold);
    printf("%%ReadsWnk - Percentage of reads containing n or more kmer from the contaminant (n=%d)\n", stats->kmer_threshold);
    printf("UniqWnk   - Reads containing n or more kmer from the contaminant and not any other (n=%d)\n", stats->kmer_threshold);
    printf("%%UniqWnk  - Percentage of reads containing n or more kmer from the contaminant and not any other (n=%d)\n", stats->kmer_threshold);
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void kmer_stats_compare_contaminant_kmers(HashTable* hash, KmerStats* stats, CmdLine* cmd_line)
{
    int i, j;
    FILE* fp_abs;
    FILE* fp_pc;
    FILE* fp_abs_unique;
    FILE* fp_pc_unique;
    char* filename_abs;
    char* filename_pc;
    char* filename_pc_unique;
    char* filename_abs_unique;
    
    if (stats->n_contaminants < 2) {
        return;
    }
    
    printf("\nComparing contaminant kmers...\n");
    
    filename_abs = malloc(strlen(cmd_line->output_prefix)+32);
    if (!filename_abs) {
        printf("Error: No room to store filename!\n");
        exit(1);
    }
    sprintf(filename_abs, "%skmer_similarity_absolute.txt", cmd_line->output_prefix);

    filename_pc = malloc(strlen(cmd_line->output_prefix)+32);
    if (!filename_pc) {
        printf("Error: No room to store filename!\n");
        exit(1);
    }
    sprintf(filename_pc, "%skmer_similarity_pc.txt", cmd_line->output_prefix);

    filename_abs_unique = malloc(strlen(cmd_line->output_prefix)+32);
    if (!filename_abs_unique) {
        printf("Error: No room to store filename!\n");
        exit(1);
    }
    sprintf(filename_abs_unique, "%skmer_unique_absolute.txt", cmd_line->output_prefix);
    
    filename_pc_unique = malloc(strlen(cmd_line->output_prefix)+32);
    if (!filename_pc_unique) {
        printf("Error: No room to store filename!\n");
        exit(1);
    }
    sprintf(filename_pc_unique, "%skmer_unique_pc.txt", cmd_line->output_prefix);
    
    fp_abs = fopen(filename_abs, "w");
    if (!fp_abs) {
        printf("Error: can't open %s\n", filename_abs);
        exit(1);
    } else {
        printf("Opened %s\n", filename_abs);
    }

    fp_pc = fopen(filename_pc, "w");
    if (!fp_pc) {
        printf("Error: can't open %s\n", filename_pc);
        exit(1);
    } else {
        printf("Opened %s\n", filename_pc);
    }

    fp_abs_unique = fopen(filename_abs_unique, "w");
    if (!fp_abs_unique) {
        printf("Error: can't open %s\n", filename_abs_unique);
        exit(1);
    } else {
        printf("Opened %s\n", filename_abs_unique);
    }
    
    fp_pc_unique = fopen(filename_pc_unique, "w");
    if (!fp_pc_unique) {
        printf("Error: can't open %s\n", filename_pc_unique);
        exit(1);
    } else {
        printf("Opened %s\n", filename_pc_unique);
    }
    
    void check_kmers_in_common(Element* node) {
        int i, j;
        for (i=0; i<(stats->n_contaminants); i++) {
            if (element_get_contaminant_bit(node, i) > 0) {
                for (j=i; j<stats->n_contaminants; j++) {
                    if (element_get_contaminant_bit(node, j) > 0) {
                        stats->kmers_in_common[i][j]++;
                        if (i != j) {
                            stats->kmers_in_common[j][i]++;
                        }
                    }
                }
            }
        }
    }
    
    hash_table_traverse(&check_kmers_in_common, hash);
    
    void check_unique_kmers(Element* node) {
        int i;
        int count = 0;
        int index = 0;
        
        for (i=0; i<(stats->n_contaminants); i++) {
            if (element_get_contaminant_bit(node, i) > 0) {
                index = i;
                count++;
                if (count > 1) {
                    break;
                }
            }
        }
        
        if (count == 1) {
            stats->unique_kmers[index]++;
        }
    }

    hash_table_traverse(&check_unique_kmers, hash);
    
    printf("\n%15s ", "");
    fprintf(fp_abs, "Contaminant");
    fprintf(fp_pc, "Contaminant");
    for (i=0; i<stats->n_contaminants; i++) {
        printf(" %15s", stats->contaminant_ids[i]);
        fprintf(fp_abs, "\t%s", stats->contaminant_ids[i]);
        fprintf(fp_pc, "\t%s", stats->contaminant_ids[i]);
        fprintf(fp_abs_unique, "\t%s", stats->contaminant_ids[i]);
        fprintf(fp_pc_unique, "\t%s", stats->contaminant_ids[i]);
    }
    printf("\n");
    fprintf(fp_abs, "\n");
    fprintf(fp_pc, "\n");
    fprintf(fp_abs_unique, "\n");
    fprintf(fp_pc_unique, "\n");
    
    for (i=0; i<stats->n_contaminants; i++) {
        printf("%15s", stats->contaminant_ids[i]);
        fprintf(fp_abs, "%s", stats->contaminant_ids[i]);
        fprintf(fp_pc, "%s", stats->contaminant_ids[i]);
        for (j=0; j<stats->n_contaminants; j++) {
            double pc = 0;
            printf(" %15d", stats->kmers_in_common[i][j]);
            fprintf(fp_abs, "\t%d", stats->kmers_in_common[i][j]);

            if (stats->kmers_in_common[i][j] > 0) {
                pc = (100.0 * (double)stats->kmers_in_common[i][j]) / (double)stats->contaminant_kmers[i];
            }

            fprintf(fp_pc, "\t%.2f", pc);
        }
        printf("\n");
        fprintf(fp_abs, "\n");
        fprintf(fp_pc, "\n");
    }

    for (i=0; i<stats->n_contaminants; i++) {
        double pc = 0;
        
        if (stats->unique_kmers[i] > 0) {
            pc = (100.0 * (double)stats->unique_kmers[i]) / (double)stats->contaminant_kmers[i];
        }
        
        if (i > 0) {
            fprintf(fp_abs_unique, "\t");
            fprintf(fp_pc_unique, "\t");
        }
        fprintf(fp_abs_unique, "%d", stats->unique_kmers[i]);
        fprintf(fp_pc_unique, "%.2f", pc);
    }
    
    fclose(fp_abs_unique);
    fclose(fp_pc_unique);
    fclose(fp_abs);
    fclose(fp_pc);
}

void kmer_stats_write_progress(KmerStats* stats, CmdLine* cmd_line)
{
    char* filename;
    FILE* fp;
    int r;
    
    printf("Updating...\n");
    
    filename = malloc(strlen(cmd_line->progress_dir) + 64);
    if (!filename) {
        printf("Error: no room for filename\n");
        exit(1);
    }
    
    for (r=0; r<stats->number_of_files; r++) {
        sprintf(filename, "%s/data_overall_r%d.txt", cmd_line->progress_dir, r+1);
        fp = fopen(filename, "w");
        if (fp) {
            printf("Opening %s\n", filename);
            fprintf(fp, "name\tvalue\n");
            fprintf(fp, "Number of reads\t%d\n", stats->read[r]->number_of_reads);
            fprintf(fp, "Number with k1 contaminants\t%d\n", stats->read[r]->k1_contaminated_reads);
            fprintf(fp, "Number with k%d contaminants\t%d\n", stats->kmer_threshold, stats->read[r]->kn_contaminated_reads);
            fclose(fp);
        } else {
            printf("Error: can't open %s\n", filename);
        }

        sprintf(filename, "%s/data_per_contaminant_r%d.txt", cmd_line->progress_dir, r+1);
        fp = fopen(filename, "w");
        if (fp) {
            int i;
            printf("Opening %s\n", filename);
            fprintf(fp, "name\tvalue\n");
            for (i=0; i<stats->n_contaminants; i++) {
                fprintf(fp, "%s\t%d\n", stats->contaminant_ids[i], stats->read[r]->kn_contaminated_reads_by_contaminant[i]);
            }
            fclose(fp);
        } else {
            printf("Error: can't open %s\n", filename);
        }

        sprintf(filename, "%s/largest_contaminant_r%d.txt", cmd_line->progress_dir, r+1);
        fp = fopen(filename, "w");
        if (fp) {
            int i;
            printf("Opening %s\n", filename);
            fprintf(fp, "name\tvalue\n");
            for (i=0; i<stats->n_contaminants; i++) {
                fprintf(fp, "%s\t%d\n", stats->contaminant_ids[i], stats->read[r]->reads_with_highest_contaminant[i]);
            }
            fprintf(fp, "Unclassified\t%d\n", stats->read[r]->reads_unclassified);
            fclose(fp);
        } else {
            printf("Error: can't open %s\n", filename);
        }
    
    }
}
