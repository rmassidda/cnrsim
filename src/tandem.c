/*
 * CNRSIM
 * tandem.c
 * Defines the structure containing
 * a tandem repeat inside the reference
 *
 * @author Riccardo Massidda
 */

#include <stdlib.h>
#include <string.h>
#include "tandem.h"

tandem_set_t * tandem_set_init ( int length, int max_motif, int max_repetition, tandem_set_t * set ){
    int tmp_size = length / 2;
    // First initialization
    if ( set == NULL ) {
        set = malloc ( sizeof ( tandem_set_t ) );
        set->set = malloc ( sizeof ( tandem_t * ) * tmp_size );
    }
    // Sequence bigger
    else if ( set->size < tmp_size ){
        set->set = realloc ( set->set, sizeof ( tandem_t * ) * tmp_size );
    }
    // Motif size
    set->max_motif = max_motif;
    // Repetitions
    set->max_repetition = max_repetition;
    // Actual size
    set->size = tmp_size;
    // Number of valid tandem repeats in the set
    set->n = 0;
    // Current tandem
    set->i = 0;

    return set;
}

tandem_set_t * tandem_set_analyze ( char * reference, int length, tandem_set_t * set ){
    int pos = 0;
    int max = 0;
    int w[set->max_motif + 1];

    while ( pos < length ){
        // Ignore 'N' regions in a chromosome
        while ( reference[pos] == 'N' ){ pos ++; }
        // Find best repetition
        max = 1;
        // Motif size € [1, max]
        for ( int i = 1; i <= set->max_motif; i ++ ){
            // Motif appears at least one time
            w[i] = 1;
            // If the sequences are comparable
            while ( ( pos + (w[i] + 1) * i ) <= length ){
                if ( strncasecmp ( &reference[pos], &reference[pos + w[i] * i], i ) == 0 ){
                    w[i] ++;
                }
                else{
                    break;
                }
            }
            // The best tandem is the wider one with at least
            // two repetitions with the same width wins 
            // the smaller motif
            if ( w[i] > 1 && (i * w[i]) > (max * w[max]) ){
                max = i;
            }
        }
        if ( w[max] > 1 ){
            if ( set->n >= set->size ){
                return set;
            }
            // Save tandem repeat
            set->set[set->n].pos = pos;
            set->set[set->n].pat = max;
            set->set[set->n].rep = w[max];
            set->n ++;
        }
        // Update position
        pos += w[max] * max;
    }

    return set;
}

void tandem_set_destroy ( tandem_set_t * set ){
    if ( set == NULL ){
        return;
    }
    free ( set->set );
    free ( set );
}
