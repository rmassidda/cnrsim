/*
 * CNRSIM
 * stats.c
 * Defines the structure containing
 * the statistics of a sequencer.
 *
 * @author Riccardo Massidda
 */

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "stats.h"

stats_t * stats_init ( ) {
    stats_t * stats = malloc ( sizeof ( stats_t ) );
    if ( stats == NULL ) return stats;

    // Data sources
    // Alignment: cigar -> cigar
    stats->alignment = source_init ( 4, 4, 2, 1 );
    // Mismatch: nucleotides -> nucleotides
    stats->mismatch = source_init ( 5, 5, 1, 0);
    // Quality: cigar -> ASCII
    stats->quality = source_init ( 4, 128, 1, 0 );
    // Distribution of errors in read
    stats->distribution = source_init ( 1, 4, 0, 0 );

    return stats;
}

static unsigned char __nucleotide ( char nucleotide ) {
    nucleotide = toupper ( nucleotide );
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

static char __nucleotide_rev ( unsigned char nucleotide ) {
    switch ( nucleotide ) {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    }
    return 'N';
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
        }
        source_update ( NULL, 0, i - 1, align[z], stats->distribution );
    }
}

read_t * stats_generate_read ( char * ref, read_t * read, stats_t * stats ){
    int i = 0;
    int pos = 0;
    unsigned char in, out;

    if ( read == NULL ){
        read = malloc ( sizeof ( read_t ) );
        read->align = NULL;
        read->alg_len = 0;
        read->read = malloc ( sizeof ( char ) * stats->quality->n );
        read->quality = malloc ( sizeof ( char ) * stats->quality->n );
    }

    // Alignment generation
    read->align = source_generate_word ( read->align, &read->alg_len, stats->alignment );
    // Read status
    read->cut = false;
   
    // Read
    for ( int z = 0; z < read->alg_len; z ++ ) {
        // Minimum position
        pos = ( i < stats->quality->n ) ? i : stats->quality->n - 1;
        pos = ( pos < stats->mismatch->n ) ? pos : stats->mismatch->n - 1;
        if ( read->align[z] != 2 ) {
            read->quality[pos] = source_generate ( &read->align[z], 1, pos, stats->quality );
        }
        switch ( read->align[z] ) {
        case 0:
            read->read[pos] = *ref;
            i++;
            ref ++;
            break;
        case 1:
            read->read[pos] = __nucleotide_rev ( rand () % 4 );
            i++;
            break;
        case 2:
            ref ++;
            break;
        case 3:
            in = __nucleotide ( *ref );
            out = source_generate ( &in, 1, pos, stats->mismatch );
            read->read[pos] = __nucleotide_rev ( out );
            i++;
            ref++;
            break;
        }
        // Reference ended
        if ( *ref == '\0' ){
            read->cut = true;
            break;
        }
    }

    // Terminal
    pos = ( i < stats->quality->n ) ? i : stats->quality->n - 1;
    pos = ( pos < stats->mismatch->n ) ? pos : stats->mismatch->n - 1;
    read->read[pos] = '\0';
    read->quality[pos] = '\0';

    return read;
}


void stats_dump ( FILE * file, stats_t * stats ) {
    source_dump ( file, "alignment", stats->alignment );
    source_dump ( file, "mismatch", stats->mismatch );
    source_dump ( file, "quality", stats->quality );
    source_dump ( file, "distribution", stats->distribution );
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

