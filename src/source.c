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
#include <assert.h>
#include <time.h>
#include <math.h>
#include "source.h"

source_t * source_init ( int sigma, int omega, int m, int graph ) {
    source_t * source = malloc ( sizeof ( source_t ) );
    if ( source == NULL ) {
        return NULL;
    }

    // Matrixes
    source->n = 0;

    // Memory
    source->m = m;

    // Alphabets
    source->sigma = sigma + graph; // exceptional start character
    source->omega = omega + graph; // exceptional end character
    source->prefix = ( int ) pow ( source->sigma, source->m );

    // Data
    source->raw = NULL;
    source->normalized = NULL;

    return source;
}

int __index ( unsigned char * in, int len, source_t * source ) {
    int i;
    int empty = source->m - len;
    int index = 0;
    int pow = 1;

    // Actual prefix
    i = 0;
    while ( i < len ) {
        index += ( in[len-1-i] * pow );
        pow *= source->sigma;
        i ++;
    }

    // Empty characters
    i = 0;
    while ( i < empty ) {
        // The last character in sigma is used to represent an uncomplete prefix
        index += ( source->sigma - 1 ) * pow;
        pow *= source->sigma;
        i ++;
    }

    return index;
}


int source_update ( unsigned char * in, int len, int pos, unsigned char out, source_t * source ) {
    if ( pos >= source->n ) {
        source->raw = realloc ( source->raw, sizeof ( unsigned long * ) * ( pos + 1 ) );

        // New matrices
        for ( int i = source->n; i <= pos; i ++ ) {
            source->raw[i] = calloc (
                                 source->omega * source->prefix,
                                 sizeof ( unsigned long )
                             );
        }
        source->n = pos + 1;
    }

    int index = __index ( in, len, source );

    //  Update stats
    source->raw[pos][ index * source->omega + out ] ++;

    return 0;
}

void source_learn_word ( unsigned char * w, int size, source_t * source ) {
    int m = source->m;
    int len;
    int i;
    for ( i = 0; i < size; i ++ ) {
        // Length of the sample
        len = ( i < m ) ? i : m;
        source_update ( &w[i - len], len, i, w[i], source );
    }

    // Add terminal char
    len = ( i < m ) ? i : m;
    source_update ( &w[i - len], len, i, source->omega - 1, source );
}

void __normalize ( source_t * source ) {
    unsigned long sum;

    // Alloc matrix
    source->normalized = malloc ( sizeof ( double * ) * source->n );
    for ( int i = 0; i < source->n; i ++ ) {
        source->normalized[i] = calloc (
                                    source->omega * source->prefix,
                                    sizeof ( double )
                                );
        for ( int j = 0; j < source->prefix; j ++ ) {
            sum = 0;
            for ( int k = 0; k < source->omega; k ++ ) {
                sum += source->raw[i][ j * source->omega + k ];
            }
            for ( int k = 0; k < source->omega; k ++ ) {
                if ( sum == 0 ) {
                    source->normalized[i][ j * source->omega + k ] = 0;
                } else {
                    source->normalized[i][ j * source->omega + k ] = ( double ) source->raw[i][ j * source->omega + k ] / sum;
                }
            }
        }
    }
}

unsigned char source_generate ( unsigned char * in, int len, int pos, source_t * source ) {
    double * p;
    double outcome;
    double threshold;

    if ( source->normalized == NULL ) {
        __normalize ( source );
    }

    // Random decision about the alternatives
    outcome = ( double ) rand() / RAND_MAX;
    threshold = 0;
    int index = __index ( in, len, source );
    p = & ( source->normalized[pos][ index * source->omega ] );
    for ( int i = 0; i < source->omega; i++ ) {
        if ( threshold <= outcome && outcome < threshold + *p ) {
            return ( unsigned char ) i;
        } else {
            threshold += *p;
        }
        p ++;
    }
    return 0;
}

unsigned char * source_generate_word ( unsigned char * w, int * size, source_t * source ) {
    int m = source->m;
    int len;
    int i;
    
    if ( w == NULL ){
        w = malloc ( sizeof ( unsigned char ) * source->n );
    }

    *size = source->n;

    for ( i = 0; i < *size; i ++ ){
        // Length of the sample
        len = ( i < m ) ? i : m;
        w[i] = source_generate ( &w[i - len], len, i, source );
        if ( w[i] == source->omega - 1 ){
            *size = i+1;
            return w;
        }
    }

    return w;
}

void source_dump ( FILE * file, char * source_name, source_t * source ) {
    fprintf ( file, "@%s %d %d %d\n", source_name, source->n, source->prefix, source->omega );

    for ( int i = 0; i < source->n; i ++ ) {
        for ( int j = 0; j < source->prefix; j ++ ) {
            for ( int z = 0; z < ( int ) source->omega; z ++ ) {
                fprintf ( file, "%lu ", source->raw[i][j * source->omega + z] );
            }
        }
        fprintf ( file, "\n" );
    }
}

void source_destroy ( source_t * source ) {
    if ( source == NULL ){
        return;
    }
    for ( int i = 0; i < source->n; i ++ ) {
        free ( source->raw[i] );
    }
    free ( source->raw );
    if ( source->normalized != NULL ) {
        for ( int i = 0; i < source->n; i ++ ) {
            free ( source->normalized[i] );
        }
        free ( source->normalized );
    }
    free ( source );
}
