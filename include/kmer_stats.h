typedef struct {
    uint32_t number_of_reads;
    
    // Counts for k1 - reads with 1 or more kmers
    uint32_t k1_contaminated_reads;
    double   k1_contaminaned_reads_pc;
    uint32_t k1_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   k1_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t k1_unique_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   k1_unique_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];

    // Counts for kn - reads with n or more kmers (n == kmer_threshold)
    uint32_t kn_contaminated_reads;
    double   kn_contaminaned_reads_pc;
    uint32_t kn_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   kn_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t kn_unique_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   kn_unique_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];

    // Contaminant kmers seen (ie. number of reference kmers found in reads)
    uint32_t contaminant_kmers_seen[MAX_CONTAMINANTS];
    double   contaminant_kmers_seen_pc[MAX_CONTAMINANTS];
    
    // Number of reads where most kmers are the contaminant
    uint32_t reads_with_highest_contaminant[MAX_CONTAMINANTS];
    uint32_t reads_unclassified;
    
    // Coverage histogram
    uint32_t contaminated_kmers_per_read[MAX_READ_LENGTH];
} KmerStatsReadCounts;

typedef struct {
    uint32_t number_of_reads;
    
    // Counts for k1 - reads with 1 or more kmers
    uint32_t both_k1_contaminated_reads;
    double   both_k1_contaminaned_reads_pc;
    uint32_t either_k1_contaminated_reads;
    double   either_k1_contaminaned_reads_pc;
    uint32_t k1_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   k1_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t k1_unique_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   k1_unique_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t k1_both_reads_contaminated;
    double   k1_both_reads_contaminated_pc;
    
    // Counts for kn - reads with n or more kmers (n == kmer_threshold)
    uint32_t both_kn_contaminated_reads;
    double   both_kn_contaminaned_reads_pc;
    uint32_t either_kn_contaminated_reads;
    double   either_kn_contaminaned_reads_pc;
    uint32_t kn_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   kn_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t kn_unique_contaminated_reads_by_contaminant[MAX_CONTAMINANTS];
    double   kn_unique_contaminated_reads_by_contaminant_pc[MAX_CONTAMINANTS];

    // Contaminant kmers seen (ie. number of reference kmers found in reads)
    uint32_t contaminant_kmers_seen[MAX_CONTAMINANTS];
    double   contaminant_kmers_seen_pc[MAX_CONTAMINANTS];
} KmerStatsBothReads;

typedef struct {
    uint32_t kmer_threshold;
    uint32_t n_contaminants;
    char* contaminant_ids[MAX_CONTAMINANTS];
    uint32_t contaminant_kmers[MAX_CONTAMINANTS];
    uint32_t unique_kmers[MAX_CONTAMINANTS];
    uint32_t kmers_in_common[MAX_CONTAMINANTS][MAX_CONTAMINANTS];
    uint32_t number_of_files;
    KmerStatsReadCounts* read[2];
    KmerStatsBothReads* both_reads;
} KmerStats;

typedef struct {
    uint32_t n_contaminants;
    uint32_t kmers_loaded;
    uint32_t contaminants_detected;
    uint32_t kmers_from_contaminant[MAX_CONTAMINANTS];
    uint32_t unique_kmers_from_contaminant[MAX_CONTAMINANTS];
    uint32_t assigned_contaminant;
    uint32_t unique_assigned_contaminant;
} KmerCounts;

void kmer_stats_initialise(KmerStats* stats, CmdLine* cmd_line);
void kmer_stats_calculate(KmerStats* stats);
void update_stats(int r, KmerCounts* counts, KmerStats* stats, KmerStatsReadCounts* both_stats);
void update_stats_for_both(KmerStats* stats, KmerStatsReadCounts* both_stats);
void kmer_stats_report_to_screen(KmerStats* stats);
void kmer_stats_compare_contaminant_kmers(HashTable* hash, KmerStats* stats, CmdLine* cmd_line);
void kmer_stats_write_progress(KmerStats* stats, CmdLine* cmd_line);