/*
 * CNRSIM
 * error_profiler.c
 * Given a BAM containing real reads
 * analyzes the distribution and the
 * entity of the errors.
 *
 * @author Riccardo Massidda
 */
#include <stdlib.h>
#include <stdio.h>
#include <htslib/sam.h>
#include "align.h"
#include "fileManager.h"

int main ( int argc, char ** argv ) {
    // FASTA
    struct filemanager * fm;
    struct sequence_t * seq;
    // BAM
    htsFile *fp;
    bam_hdr_t *hdr;
    bam1_t *line;
    hts_idx_t *index;
    hts_itr_t *itr;

    if ( argc != 4 ) {
        printf ( "usage: %s fasta bam out\n", argv[0]);
        exit ( EXIT_FAILURE );
    }

    // FASTA file
    fm = filemanager_init ( argv[1] );
    if ( fm == NULL ) {
        exit ( EXIT_FAILURE );
    }

    // BAM file
    fp = hts_open ( argv[2], "r" );
    hdr = sam_hdr_read ( fp );
    line = bam_init1 ();

    // BAM index
    index = bam_index_load ( argv[2] );
    if ( index == NULL ){
        // File not indexed
        perror ( "Can't load BAM index" );
        exit ( EXIT_FAILURE );
    }
    

    // Load of the first sequence
    seq = filemanager_next_seq ( fm, NULL );
    
    // While there are sequences to read in the FASTA file
    while ( seq != NULL ) {
        // Seek on the BAM
        itr = bam_itr_querys( index, hdr, seq->label);
        if ( itr != NULL ){
            printf ( " %s found\n", seq->label );
            while( bam_itr_next( fp, itr, line ) > 0){
                int32_t pos = line->core.pos;
                char *chr = hdr->target_name[line->core.tid];
                uint32_t len = line->core.l_qseq; //length of the read.
                uint8_t *q = bam_get_seq(line); //quality string
                        
                char *qseq = (char *)malloc(len);

                for(int i=0; i< len ; i++){
                    qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //converts into IUPAC id.
                }
                printf ( "%s\t%d\n", chr, pos );
                printf("%s\t%d\t%d\t%s\n",chr,pos,len,qseq);
            }
        }		
        else{
            printf ( "%s not found\n ", seq -> label);
        }
        // Next sequence
        seq = filemanager_next_seq ( fm, seq );
    }
    
    // Cleanup
    filemanager_destroy ( fm );
    bam_destroy1( line );
    bam_hdr_destroy ( hdr );
    bam_itr_destroy ( itr );
    sam_close( fp );
    return 0;
}
