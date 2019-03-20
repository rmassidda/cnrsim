/*
 * CNRSIM
 * source.c
 * Finite memory sequence generator
 *
 * @author Riccardo Massidda
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "source.h"

source_t * source_init ( char * sigma, char * omega, int m, int k ){
    source_t * source = malloc ( sizeof ( source_t ) );
    if ( source == NULL ){
        return NULL;
    }

    // Memory
    source->m = m;
    // Precision
    source->k = k;
    // Matrixes
    source->n = 1;

    // Alphabets
    source->sigma = malloc ( strlen ( sigma ) );
    strcpy ( source->sigma, sigma );
    source->sigma_size = pow ( strlen ( sigma ), m );
    source->omega = malloc ( strlen ( omega ) );
    strcpy ( source->omega, omega );
    source->omega_size = strlen ( omega );

    // Data
    source->raw = realloc ( source->raw, sizeof ( unsigned long * ) );
    source->raw[0] = calloc ( 
            source->omega_size * source->sigma_size,
            sizeof ( unsigned long )
            );

    return source;
}

int source_update ( char * prefix, int pos, char out, source_t * source ){
    int index = pos % source->k;

    if ( index >= source->n ){
        unsigned long ** tmp_raw;

        // Safe realloc
        tmp_raw = realloc ( source->raw, sizeof ( unsigned long * ) * ( index + 1 ) );
        if ( tmp_raw == NULL ){
            return -1;
        }
        source->raw = tmp_raw;

        // New matrices
        for ( int i = source->n; i <= index; i ++ ){
            source->raw[i] = calloc ( 
                    source->omega_size * source->sigma_size,
                    sizeof ( unsigned long )
                    );
            if ( source->raw[i] == NULL ){
                return -1;
            }
        }
        source->n = index + 1;
    }

    // This is where the magic happens!
    int hash_prefix = 0;
    int hash_out = 0;

    //  Update stats
    source->raw[index][ hash_prefix * source->sigma_size + hash_out ] ++;
        
    return 0;
}

char source_generate ( char * prefix, int pos, source_t * source );

void source_destroy ( source_t * source ){
    for ( int i = 0; i < source->n; i ++ ){
        free ( source->raw[i] );
    }
    free ( source->raw );
    free ( source );
}
