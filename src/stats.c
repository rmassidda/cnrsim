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
        stats->read = NULL;
        stats->max_read = 0;
        stats->tot = 0;
    }
    return stats;
}

int _b2i ( char a ){
    switch ( a ){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
        default: return 0;
    }
}

bool stats_update ( stats_t * stats, unsigned char * align, int alg_len, char * read, int read_len, char * reference ){
    // Catch indels
    int start_in = 0;
    int start_del = 0;
    int ref_pos = 0;
    int read_pos = 0;

    // Memory must be reallocated
    if ( read_len > stats->max_read ){
        unsigned long * read_tmp;
        
        read_tmp = realloc ( stats->read, sizeof ( unsigned long ) * ( read_len ) );

        // Check realloc
        if ( read_tmp == NULL ) {
            return false;
        }
        else{
            stats->max_read = read_len;
            stats->read = read_tmp;
        }
    }

    stats->tot ++;
    stats->read[read_len - 1] ++;

    // Variation stats
    for ( int i = 0; i < alg_len; i ++ ){
        // Catch end of INDEL
        if ( start_in && align[i] != I ){
            start_in = 0;
            // Inserzione dei caratteri read[ read_pos - start_in, read_pos - 1 ] 
        }
        else if ( start_del && align[i] != D ){
            start_del = 0;
            // Delezione dei caratteri ref[ ref_pos - start_del, ref_pos - 1 ]
        }
        
        if ( align[i] == M ){
            read_pos ++;
            ref_pos ++;
        }
        else if ( align[i] == X ){
            read_pos ++;
            ref_pos ++;
        }
        else if ( align[i] == I ){
            read_pos ++;
            start_in ++;
        }
        else{
            ref_pos ++;
            start_del ++;
        }
    }

    return true;
}

void stats_dump ( FILE * file, stats_t * stats ) {
    fprintf ( file, "%ld\n", stats->tot );
}

void stats_destroy ( stats_t * stats ){
    free ( stats->read );
    free ( stats );
}

