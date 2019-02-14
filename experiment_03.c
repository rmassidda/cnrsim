/*
 * Third experiment with HTSLIB.
 * Try to read parse both VCF and FA file.
 *
 * Achievement:
 * - First "useful" toy 
 */
#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>
#include "fileManager.h"

struct filemanager *fm;
struct sequence_t *seq;

int main(int argc, char **argv){
    int fa_start = 0;
    int fa_limit = 0;

    // Filenames from argument
    if (argc != 3){
        printf("usage: e03 vcf_filename fasta_filename\n");
        exit(EXIT_FAILURE);
    }

    // Open VCF file
    htsFile * inf = bcf_open (argv[1], "r" );
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

    // Struct containing the FASTA sequence
    seq = NULL;
    // Iteration on the line set
    while (bcf_read(inf, hdr, line) == 0){
        // Unpack 'till ALT 
        if (bcf_unpack( line, BCF_UN_STR) != 0){
            printf("Unpack error");
            return -1;
        }
        // line->pos            Position of the variation
        // line->d.allele[0]    Reference
        printf("%d\t", line->pos);
        printf("%s\t", line->d.allele[0]);
        while (line->pos > fa_limit) {
          seq = filemanager_next_seq (fm, seq);
          fa_start = fa_limit;
          fa_limit += seq->sequence_size;
        }
        for(int i=0; i < strlen(line->d.allele[0]); i++){
            printf("%c", seq->sequence[line->pos + fa_start + i]);
        }
        printf("\n");
    }

    // Cleanup
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    filemanager_destroy (fm);

    exit(EXIT_SUCCESS);
}
