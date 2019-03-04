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
#include <math.h>
#include "align.h"
#include "fileManager.h"
#include "translate_notation.h"

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
    // Alias dictionary
    region_index_t * alias_index;
    char * alias = NULL;
    // Aligner
    aligner_t * aligner = NULL;
    char * read = NULL;
    int pos = 0;
    int len = 0;
    int start = 0;
    int end = 0;

    if ( argc != 4 ) {
        printf ( "usage: %s fasta bam dictionary\n", argv[0]);
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

    // Alias dictionary
    alias_index = tr_init ( argv[3] );
    if ( alias_index == NULL ){
        perror ( "Can't load alias dictionary" );
        exit ( EXIT_FAILURE );
    }

    // Load of the first sequence
    seq = filemanager_next_seq ( fm, NULL );
    
    // While there are sequences to read in the FASTA file
    while ( seq != NULL ) {
        // Seek on the BAM
        alias = tr_translate ( alias_index, seq->label );
        itr = bam_itr_querys( index, hdr, alias);
        if ( itr != NULL ){
            printf ( " %s found\n", alias );
            while( bam_itr_next( fp, itr, line ) > 0){
                // Read information
                pos = line->core.pos;
                len = line->core.l_qseq;
                uint8_t *q = bam_get_seq( line ); //quality string        
                
                // Interval of the reference
                // start = floor ( pos - len/4 );
                // end = floor ( pos + len/4 );
                // if ( start <= 0 ){
                //     start = 0;
                // }
                // if ( end >= seq->sequence_size){
                //     end = seq->sequence_size - 1;
                // }
                start = pos;
                end = pos + len;

                // Read string
                read = realloc ( read, sizeof ( char ) * ( len + 1 ) );
                int i;
                for ( i = 0; i < len; i++ ){
                    read[i] = seq_nt16_str[ bam_seqi ( q, i ) ];
                }
                read[i] = 0;

                // Align
                aligner = al_init ( aligner, seq->sequence, end, read );

                // Useless print
                printf ( "%s\t%d\t%d\t%d\t%d\n", alias, pos, len, start, end );
                printf ( "%.*s\n", ( end - start ), &seq->sequence[start]);
                printf ( "%s\n", build_alignment( aligner ) );
                printf("%s\n", read );
            }
        }		
        else{
            printf ( "%s not found\n ", alias );
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
    free ( read );
    tr_destroy ( alias_index );
    return 0;
}
