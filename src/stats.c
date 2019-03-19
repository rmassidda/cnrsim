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
    if ( stats != NULL ){
        for ( int i = 0; i < N_STAT; i ++ ){
            stats->count[i] = NULL;
        }
        stats->size = 0;
        stats->tot = 0;
    }
    return stats;
}

bool stats_update ( stats_t * stats, unsigned char * align, int alg_len, char * read, int read_len, char * reference ){
    // Memory must be reallocated
    // alg_len >= read_len
    if ( alg_len > stats->size ){
        for ( int i = 0; i < N_STAT; i ++ ){
            unsigned long * tmp;
            tmp = realloc ( stats->count[i], sizeof ( unsigned long ) * ( alg_len ) );
            // Check realloc
            if ( tmp == NULL ) {
                return false;
            }
            stats->count[i] = tmp;
            // Clean memory
            // for ( int j = stats->size; j < alg_len; j ++ ){
            //     stats->count[i][j] = 0;
            // }
            memset ( &stats->count[i][stats->size], 0, sizeof ( unsigned long )  * ( alg_len - stats->size ) );
        }
    }

    stats->size = alg_len;
    stats->tot ++;
    stats->count[L][read_len - 1] ++;

    // Variation stats
    for ( int i = 0; i < alg_len; i ++ ){
        stats->count[align[i]][i] ++;
    }

    return true;
}

void stats_dump ( FILE * file, stats_t * stats ) {
    fprintf ( file, "%ld\n", stats->tot );
    printf ( "P\tM\tD\tI\tX\tL\n" );
    for ( int i = 0; i < stats->size; i ++ ){
        printf ( "%d\t", i );
        for ( int j = 0; j < N_STAT; j ++ ){
            printf ( "%ld\t", stats->count[j][i] );
        }
        printf ( "\n" );
    }
}

void stats_destroy ( stats_t * stats ){
    for ( int i = 0; i < N_STAT; i ++ ){
        free ( stats->count[i] );
    }
    free ( stats );
}

