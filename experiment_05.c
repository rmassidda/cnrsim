/*
 * Fifth experiment. 
 * The variator is ready!
 * But.. it only works for one sequence at the time!
 * What if the VCF or the FASTA files containes more subsequences?
 * This is the case of a whole genome where different chromosomes
 * are separated inside a multisequence FASTA file.
 *
 * Achievement:
 * - Parse multisequence FASTA files
 * - Parse multisequence VCF files
 */

#include <htslib/synced_bcf_reader.h>
#include <stdio.h>
#include "fileManager.h"


int main( int argc, char ** argv ){
    // FASTA structures
    // struct filemanager *fm;
    // struct sequence_t *seq;
    // FASTA filename
    // char * filename = argv[1]; 
    // VCF 
    bcf_srs_t *sr =  bcf_sr_init() ;
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    bcf_sr_add_reader (sr, argv[1] );
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    bcf1_t * line; 
    // Open File
    // fm = filemanager_init (filename);
    // // FASTA parse
    // seq = filemanager_next_seq (fm, NULL);
    // while (seq != NULL) {
    //   // Print sequence name
    //   printf("%s\n",seq->label);
    //   seq = filemanager_next_seq (fm, seq);
    // };
    if ( sr->regions != NULL ){
        for ( int i = 0; i < sr->regions->nseqs; i++){
            printf("%s\t", sr->regions->seq_names[i]);
            int res = bcf_sr_seek ( sr, sr->regions->seq_names[i], 0 );    
            printf("%d\t", res);
            if(bcf_sr_next_line (sr)) { //loop through file
                line =  bcf_sr_get_line(sr, 0);  //read a line
                printf ( "%s\t%d\n", bcf_hdr_id2name( hdr, line->rid  ), line->pos );
            }
        }
    }
    else{
        int res = bcf_index_build2(argv[1], NULL, 0);
        printf( "Build index %d\n", res );
    }

    // filemanager_destroy (fm);  /*  Releasefile managemant resources  */
    return 0;
}
