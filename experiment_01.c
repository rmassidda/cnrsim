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
    // Struct for storing each record
    bcf1_t *rec = bcf_init();
    if (rec == NULL){
        return -1;
    }
    int record = 1;
    // Iteration on the record set
    while (bcf_read(inf, hdr, rec) == 0){
        if (bcf_unpack(rec,record) != 0){
            printf("Unpack error");
            return -1;
        }
        printf("%d",rec->pos);
        int i = 0;
        while ( rec->d.allele[i] != NULL ){
            printf("\t%s", rec->d.allele[i]);
            i++;
        }
        printf("\n");
        record++;
    }


    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
}
