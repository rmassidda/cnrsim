/*
 * Second experiment with HTSLIB.
 * Try to read secondary informations from a VCF file.
 * This could be useful to parse genotypes, allel frequency and so on. 
 *
 * Achievement:
 * - Truly understand VCF file format
 * - Wider use of vcf.h and vcfutils.h 
 */

#include <string.h>
#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
    // float *af = NULL;
    // int af_size = 0;
    // int af_ret;
    char *freq = NULL;
    int freq_size = 0;
    int freq_ret;
    char *substr;
    // Filename from argument
    if (argc != 2){
        printf("Filename required.\n");
        return -1;
    }
    // Open File
    htsFile * inf = bcf_open (argv[1], "r" );
    // Read of VCF header
    bcf_hdr_t * hdr = bcf_hdr_read(inf);
    if ( hdr == NULL ){
        return -1;
    } 
    // Struct for storing each line
    bcf1_t * line = bcf_init();
    if (line == NULL){
        return -1;
    }
    // Iteration on the line set
    while (bcf_read(inf, hdr, line) == 0){
        // Unpack 'till ALT 
        if (bcf_unpack( line, BCF_UN_INFO) != 0){
            printf("Unpack error");
            return -1;
        }
        printf("%d", line->pos);
        // Definito dallo standard
        // af_ret = bcf_get_info_float( hdr, line, "AF", af, &af_size );
        // Usato da dbSNP
        freq_ret = bcf_get_info_string( hdr, line, "FREQ", &freq, &freq_size );
        if (freq_ret >= 0){
            substr = strtok(freq, ":");
            printf("\t%s",substr);
        }
        for (int i = 0; i < line->n_allele; i++){
            printf("\t%s", line->d.allele[i]);
            if (substr != NULL){
                substr = strtok(NULL, ",");
                printf(":%s",substr);
            }
        }
        printf("\n");
    }

    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    return 0;
}
