/*
 * CNRSIM
 * stats.c
 * Defines the structure containing
 * the statistics of a sequencer.
 *
 * @author Riccardo Massidda
 */

#include <stdlib.h>
#include <string.h>
#include "stats.h"

stats_t * stats_init ( ){
    stats_t * stats = malloc ( sizeof ( stats_t ) );
    if ( stats == NULL ) return stats;
    
    // Alphabets
    alphabet_t * alignment = NULL;
    alphabet_t * nucleotides = alphabet_init ( ( unsigned char * )"ACGTN", 5 );
    alphabet_t * quality = alphabet_init (
            ( unsigned char * ) "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~",
            94 );

    // Data sources
    stats->alignment = source_init ( alignment, alignment, 1, 1 );
    stats->mismatch = source_init ( nucleotides, nucleotides, 1, 1 );
    stats->quality = source_init ( alignment, quality, 1, 1 );

    return stats;
}

bool stats_update ( stats_t * stats, unsigned char * align, int alg_len, char * read, int read_len, char * reference ){
    int read_pos = 0;
    int ref_pos = 0;

    for ( int i = 0; i < alg_len; i ++ ){
        switch ( align[i] ){
            case MAT:
                read_pos ++;
                ref_pos ++;
                break;
            case MIS:
                source_update (
                        ( unsigned char * ) &reference[ref_pos],
                        read_pos,
                        read[read_pos],
                        stats->mismatch
                        );
                read_pos ++;
                ref_pos ++;
                break;
            case INS:
                read_pos ++;
                break;
            case DEL:
                ref_pos ++;
                break;
        }
        source_update (
                &align[i-1],
                i,
                align[i],
                stats->alignment
        );
    }
    return true;
}

void stats_dump ( FILE * file, stats_t * stats ) {}

void stats_destroy ( stats_t * stats ) {
    // Free alphabets
    free ( stats->alignment->omega );
    free ( stats->alignment->sigma );
    free ( stats->mismatch->omega );
    free ( stats->mismatch->sigma );
    // Free sources
    source_destroy ( stats->alignment ); 
    source_destroy ( stats->mismatch ); 
    // Free structure
    free ( stats );
}

