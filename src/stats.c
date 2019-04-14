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
    // Alignment: alignment -> aligment
    stats->alignment = source_init ( 4, 5, 2 );

    return stats;
}

void stats_update ( unsigned char * align, int alg_len, stats_t * stats ){
    // Alignment string
    source_learn_word ( align, alg_len, stats->alignment );
}

void stats_dump ( FILE * file, stats_t * stats ) {}

void stats_destroy ( stats_t * stats ) {
    // Free sources
    source_destroy ( stats->alignment ); 
    // Free structure
    free ( stats );
}

