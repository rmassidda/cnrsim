    /*
     * First experiment with HTSLIB.
     * Try to read variations from chrY.vcf and filter the results.
     *
     * Achievement:
     * - try to understand vcf.h
     */

#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
    // Filename from argument
    if (argc != 2){
        printf("Filename required.\n");
        return -1;
    }
    // Open File
    htsFile * inf = bcf_open (argv[1], "r" );
    // Read of VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    if ( hdr == NULL ){
        return -1;
    } 
    // Struct for storing each line
    bcf1_t *line = bcf_init();
    if (line == NULL){
        return -1;
    }
    // Iteration on the line set
    while (bcf_read(inf, hdr, line) == 0){
        // Unpack 'till ALT 
        if (bcf_unpack( line, BCF_UN_STR) != 0){
            printf("Unpack error");
            return -1;
        }
        printf("%d %d:",line->pos, line->n_allele);
        for (int i = 0; i < line->n_allele; i++){
            printf("\t%s", line->d.allele[i]);
        }
        printf("\n");
    }


    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
}
