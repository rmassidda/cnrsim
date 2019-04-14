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
#include <time.h>
#include <math.h>
#include "source.h"

source_t * source_init ( int sigma, int omega, int m ){
    source_t * source = malloc ( sizeof ( source_t ) );
    if ( source == NULL ){
        return NULL;
    }

    // Matrixes
    source->n = 0;

    // Memory
    source->m = m;

    // Alphabets
    source->sigma = sigma + 1; // epsilon for uncomplete prefix
    source->omega = omega;

    // Data
    source->raw = NULL;
    source->normalized = NULL;

    return source;
}

int source_update ( unsigned char * in, int len, int pos, unsigned char out, source_t * source ){
    int index = 0;
    int char_value = 0;

    if ( pos >= source->n ){
        source->raw = realloc ( source->raw, sizeof ( unsigned long * ) * ( pos + 1 ) );

        // New matrices
        for ( int i = source->n; i <= pos; i ++ ){
            source->raw[i] = calloc ( 
                    source->omega * source->sigma,
                    sizeof ( unsigned long )
                    );
        }
        source->n = pos + 1;
    }

    for ( int i = 0; i < source->m; i ++ ){
        if ( i < len )
            char_value = in[len - i - 1];
        else
            char_value = source->sigma - 1;

        index += char_value * pow ( source->sigma, i );
    }

    //  Update stats
    source->raw[pos][ index * source->omega + out ] ++;

    return 0;
}

void __normalize ( source_t * source ) {
    double sum;

    // Alloc matrix
    source->normalized = malloc ( sizeof ( double * ) * source->n );
    for ( int i = 0; i < source->n; i ++ ){
        source->normalized[i] = calloc ( 
                source->omega * source->sigma,
                sizeof ( double )
                );
        for ( int j = 0; j < source->sigma; j ++ ){
            sum = 0;
            for ( int k = 0; k < source->omega; k ++ ) {
                sum += source->raw[i][ j * source->omega + k ];
            }
            for ( int k = 0; k < source->omega; k ++ ) {
                if ( sum == 0 ){
                    source->normalized[i][ j * source->omega + k ] = 0;
                }
                else {
                    source->normalized[i][ j * source->omega + k ] = source->raw[i][ j * source->omega + k ] / sum;
                }
            }
        }
    }
}

unsigned char source_generate ( unsigned char in, int pos, source_t * source ){
    double * p;
    double outcome;
    double threshold;

    if ( source->normalized == NULL ){
        __normalize( source );
        srand ( time ( NULL ) );
    }

    // Random decision about the alternatives
    outcome = ( double ) rand() / RAND_MAX;
    threshold = 0;
    p = & ( source->normalized[pos][ in * source->omega ] );
    for ( int i = 0; i < source->omega; i++ ) {
        //printf ( "%d -> %d %f %f %f\n", in, i, threshold, outcome, threshold + *p );
        if ( threshold <= outcome && outcome < threshold + *p ) {
            return ( unsigned char ) i;
        } else {
            threshold += *p;
        }
        p ++;
    }
    return 0;
}

void source_destroy ( source_t * source ){
    for ( int i = 0; i < source->n; i ++ ){
        free ( source->raw[i] );
    }
    free ( source->raw );
    free ( source );
}
