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

    // Data sources
    // Alignment: cigar -> cigar
    stats->alignment = source_init ( 4, 4, 2 );
    // Mismatch: nucleotides -> nucleotides
    stats->mismatch = source_init ( 5, 5, 1 );

    return stats;
}

static unsigned char __nucleotide ( char nucleotide ){
    switch ( nucleotide ){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }
    return 4;
}

void stats_update ( unsigned char * align, int alg_len, char * read, char * ref, stats_t * stats ){
    // Pointers
    char * ptr_read = read;
    char * ptr_ref = ref;
    int i = 0;
    unsigned char in, out;

    // Alignment string
    source_learn_word ( align, alg_len, stats->alignment );

    // Read
    for ( int z = 0; z < alg_len; z ++ ){
        switch ( align[z] ){
            case 0: i++; ptr_read ++; ptr_ref ++; break;
            case 1: i++; ptr_read ++; break;
            case 2: ptr_ref ++; break;
            case 3: {
                in = __nucleotide ( *ptr_ref );
                out = __nucleotide ( *ptr_read );
                source_update ( &in, 1, i, out, stats->mismatch );
                i++;
                ptr_ref++; 
                ptr_read ++; 
                break;
            }
        }
    }
}

void stats_dump ( FILE * file, stats_t * stats ) {}

void stats_destroy ( stats_t * stats ) {
    // Free sources
    source_destroy ( stats->alignment ); 
    source_destroy ( stats->mismatch ); 
    // Free structure
    free ( stats );
}

