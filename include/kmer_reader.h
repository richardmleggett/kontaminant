typedef enum
{
    UNSPECIFIED_FORMAT   = 0,
    FASTA                = 1,
    FASTQ                = 2,
    CTX                  = 3,
    ROCHE                = 4,
    HASH                 = 5,
    CSFASTA              = 6,
    KMERS                = 7,
    FILE_FORMAT_LAST     = 8,
} FileFormat ;

typedef struct{
    FILE* input_fp;
    FILE* output_fp;
    FILE* removed_fp;
    Sequence * seq;
    int max_read_length;
    boolean new_entry;
    boolean full_entry;
    FileFormat format;
    short kmer_size;
} KmerFileReaderWrapperArgs;

typedef struct {
    char* input_filename;
    char* output_filename;
    char* removed_filename;
    FileFormat format;
    short colour;
    long long bad_reads;
    char quality_cut_off;
    int max_read_length;
    int fastq_ascii_offset;
    float maximum_ocupancy;
    boolean stop_on_full;
    boolean insert;
    CmdLine *cmd_line;
    HashTable * KmerHash;
} KmerFileReaderArgs;

uint32_t load_kmer_library(char* filename, int n, int k, HashTable* contaminant_hash);
long long screen_kmers_from_file(KmerFileReaderArgs* fra, CmdLine* cmd_line, KmerStats* stats);
long long screen_or_filter_paired_end(KmerFileReaderArgs* fra_1, KmerFileReaderArgs* fra_2, KmerStats* stats);