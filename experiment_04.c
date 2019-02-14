/*
 * Fourth experiment: a variator!
 * For every variation the variator randomly decides
 * do implement it in the reference.
 *
 * Achievement:
 * - First variant generator
 */
#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>
#include "fileManager.h"

struct filemanager *fm;
struct sequence_t *seq;

int main(int argc, char **argv){
    int fa_offset = 0;
    FILE *log;
    char *str;
    // TOFIX: 1024 is totally arbitrary
    float p = 0;
    // Initialize "random" (because of the constant seed) generator
    srandom(42);

    // Filenames from argument
    if (argc != 4){
        printf("usage: e04 vcf_filename fasta_filename probability\n");
        exit(EXIT_FAILURE);
    }
    
    // LOGFILE
    log = fopen("e04.log", "w+");


    // Open VCF file
    htsFile * inf = bcf_open (argv[1], "r" );
    // TODO: Manage errors?
    
    // Read of VCF header
    bcf_hdr_t * hdr = bcf_hdr_read(inf);
    if ( hdr == NULL ){
        perror("Error reading VCF header");
        exit(EXIT_FAILURE);
    } 
    // Struct for storing each line
    bcf1_t * line = bcf_init();
    if (line == NULL){
        perror("Error creating line object");
        exit(EXIT_FAILURE);
    }

    //  Open FASTA file
    fm = filemanager_init (argv[2]);
    if (fm == NULL){
        perror("Error reading FASTA file");
        exit(EXIT_FAILURE);
    }

    // User-provided probability
    p = atof(argv[3]);

    // First read from the FASTA file
    // It's also the unique read, this is because
    // the variator doesn't support multisequences at this time
    // All the test have been executed
    // on a single chromosome
    seq = NULL;
    seq = filemanager_next_seq (fm, seq);
    // Iteration on the line set
    while (bcf_read(inf, hdr, line) == 0){
        // Unpack 'till ALT 
        if (bcf_unpack( line, BCF_UN_STR) != 0){
            printf("Unpack error");
            return -1;
        }
        
        int width = line->pos - fa_offset;
        if (width < 0){
            continue;
        }
        printf("%.*s",
                line->pos - fa_offset,
                seq->sequence + fa_offset); 
        // Probability test
        if ( ((float)random())/RAND_MAX < p ){
            // Print of a variation
            // TODO: Currently just the first one
            str = line->d.allele[1];
        }
        else{
            // Print of the reference
            str = line->d.allele[0];
        }

        printf("%s", str);
        fa_offset = line->pos + strlen(str);
    }
    // Cleanup
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    filemanager_destroy (fm);
    fclose(log);
    exit(EXIT_SUCCESS);
}
