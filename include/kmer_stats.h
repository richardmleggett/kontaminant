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

    // For parallel access
    pthread_mutex_t lock;
} KmerStatsReadCounts;

typedef struct {
    uint32_t number_of_reads;
    
    // Contaminant kmers seen (ie. number of reference kmers found in reads)
    uint32_t contaminant_kmers_seen[MAX_CONTAMINANTS];
    double   contaminant_kmers_seen_pc[MAX_CONTAMINANTS];

    uint32_t threshold_passed_reads;
    double   threshold_passed_reads_pc;
    uint32_t k1_both_reads_not_threshold;
    double   k1_both_reads_not_threshold_pc;
    uint32_t k1_either_read_not_threshold;
    double   k1_either_read_not_threshold_pc;
    
    uint32_t threshold_passed_reads_unique;
    double   threshold_passed_reads_pc_unique;
    uint32_t k1_both_reads_not_threshold_unique;
    double   k1_both_reads_not_threshold_pc_unique;
    uint32_t k1_either_read_not_threshold_unique;
    double   k1_either_read_not_threshold_pc_unique;
    
    uint32_t threshold_passed_reads_by_contaminant[MAX_CONTAMINANTS];
    double   threshold_passed_reads_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t k1_both_reads_not_threshold_by_contaminant[MAX_CONTAMINANTS];
    double   k1_both_reads_not_threshold_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t k1_either_read_not_threshold_by_contaminant[MAX_CONTAMINANTS];
    double   k1_either_read_not_threshold_by_contaminant_pc[MAX_CONTAMINANTS];
    
    uint32_t threshold_passed_reads_unique_by_contaminant[MAX_CONTAMINANTS];
    double   threshold_passed_reads_unique_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t k1_both_reads_not_threshold_unique_by_contaminant[MAX_CONTAMINANTS];
    double   k1_both_reads_not_threshold_unique_by_contaminant_pc[MAX_CONTAMINANTS];
    uint32_t k1_either_read_not_threshold_unique_by_contaminant[MAX_CONTAMINANTS];
    double   k1_either_read_not_threshold_unique_by_contaminant_pc[MAX_CONTAMINANTS];
    
    boolean filter_read;

    // For parallel access
    pthread_mutex_t lock;
} KmerStatsBothReads;

typedef struct {
    uint32_t n_contaminants;
    char* contaminant_ids[MAX_CONTAMINANTS];
    uint32_t contaminant_kmers[MAX_CONTAMINANTS];
    uint32_t unique_kmers[MAX_CONTAMINANTS];
    uint32_t kmers_in_common[MAX_CONTAMINANTS][MAX_CONTAMINANTS];
    uint32_t number_of_files;
    KmerStatsReadCounts* read[2];
    KmerStatsBothReads* both_reads;

    // For parallel access
    pthread_mutex_t lock;
} KmerStats;

typedef struct {
    uint32_t n_contaminants;
    uint32_t kmers_loaded;
    uint32_t contaminants_detected;
    uint32_t kmers_from_contaminant[MAX_CONTAMINANTS];
    uint32_t unique_kmers_from_contaminant[MAX_CONTAMINANTS];
    uint32_t assigned_contaminant;
    uint32_t unique_assigned_contaminant;

    // For parallel access
    pthread_mutex_t lock;
} KmerCounts;

void kmer_stats_initialise(KmerStats* stats, CmdLine* cmd_line);
void kmer_stats_calculate(KmerStats* stats);
void update_stats(int r, KmerCounts* counts, KmerStats* stats, CmdLine* cmd_line);
void update_stats_parallel(int r, KmerCounts* counts, KmerStats* stats, CmdLine* cmd_line);
boolean update_stats_for_both(KmerStats* stats, CmdLine* cmd_line, KmerCounts* counts_a, KmerCounts* counts_b);
boolean update_stats_for_both_parallel(KmerStats* stats, CmdLine* cmd_line, KmerCounts* counts_a, KmerCounts* counts_b);
void kmer_stats_report_to_screen(KmerStats* stats, CmdLine* cmd_line);
void kmer_stats_compare_contaminant_kmers(HashTable* hash, KmerStats* stats, CmdLine* cmd_line);
void kmer_stats_write_progress(KmerStats* stats, CmdLine* cmd_line);