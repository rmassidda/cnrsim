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

stats_t * stats_init ( ) {
    stats_t * stats = malloc ( sizeof ( stats_t ) );
    if ( stats == NULL ) return stats;

    // Data sources
    // Alignment: cigar -> cigar
    stats->alignment = source_init ( 4, 4, 2 );
    // Mismatch: nucleotides -> nucleotides
    stats->mismatch = source_init ( 5, 5, 1 );
    // Quality: cigar -> ASCII
    stats->quality = source_init ( 4, 128, 1 );
    // Distribution of errors in read
    stats->distribution = source_init ( 0, 4, 0 );

    return stats;
}

static unsigned char __nucleotide ( char nucleotide ) {
    switch ( nucleotide ) {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    }
    return 4;
}

void stats_update ( unsigned char * align, int alg_len, char * read, char * ref, unsigned char * quality, stats_t * stats ) {
    // Pointers
    char * ptr_read = read;
    char * ptr_ref = ref;
    int i = 0;
    unsigned char in, out;

    // Alignment string
    source_learn_word ( align, alg_len, stats->alignment );

    // Read
    for ( int z = 0; z < alg_len; z ++ ) {
        switch ( align[z] ) {
        case 0:
            i++;
            ptr_read ++;
            ptr_ref ++;
            break;
        case 1:
            i++;
            ptr_read ++;
            break;
        case 2:
            ptr_ref ++;
            break;
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
        if ( align[z] != 2 ) {
            source_update ( &align[z], 1, i - 1, quality[i - 1], stats->quality );
            source_update ( NULL, 0, i - 1, align[z], stats->distribution );
        }
    }
}

void stats_dump ( FILE * file, stats_t * stats ) {
    fprintf ( file, "//Alignment prefix %d\n", stats->alignment->m );
    source_dump ( file, stats->alignment );
    fprintf ( file, "//Mismatch\n" );
    source_dump ( file, stats->mismatch );
    fprintf ( file, "//Quality score\n" );
    source_dump ( file, stats->quality );
    fprintf ( file, "//Length\n" );
    source_dump ( file, stats->distribution );
}

void stats_destroy ( stats_t * stats ) {
    // Free sources
    source_destroy ( stats->alignment );
    source_destroy ( stats->mismatch );
    source_destroy ( stats->quality );
    source_destroy ( stats->distribution );
    // Free structure
    free ( stats );
}

