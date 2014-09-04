#define MAX_PATH_LENGTH 1024
#define DO_SCREEN 1
#define DO_FILTER 2
#define DO_INDEX 3

typedef struct {
    int bucket_bits;
    int bucket_size;
    int run_type;
    int kmer_size;
    int format;
    int quality_score_offset;
    int max_read_length;
    int kmer_threshold;
    boolean keep_contaminated_reads;
    boolean write_progress_file;
    int progress_delay;
    char* input_filename_one;
    char* input_filename_two;
    char* file_of_files;
    char* progress_dir;
    char* contaminants;
    char* contaminant_dir;
    char* output_prefix;
    char* removed_prefix;
    double subsample_ratio;
} CmdLine;

void initialise_cmdline(CmdLine* c);
void parse_command_line(int argc, char* argv[], CmdLine* c);