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

alphabet_t * alphabet_init ( unsigned char * symbols, int length ){
    alphabet_t * alphabet = malloc ( sizeof ( alphabet_t ) );
    if ( alphabet == NULL ){
        return NULL;
    }
    alphabet->symbols = malloc ( sizeof ( unsigned char ) * length );
    memcpy ( alphabet->symbols, symbols, sizeof ( unsigned char ) * length );
    alphabet->length = length;
    return alphabet;
}

int alphabet_hash ( unsigned char * word, int length, alphabet_t * alphabet ){
    int c = 0;
    for ( int i = 0; i < length; i ++ ){
        int h = 0;
        for ( int j = 0; j < alphabet->length; j ++ ){
            if ( alphabet->symbols[j] == word[i] ){
                h = i;
                break;
            }
        }
        c += ( h * pow ( alphabet->length, length - i - 1 ) );
    }
    return c;
}

source_t * source_init ( alphabet_t * sigma, alphabet_t * omega, int m, int k ){
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
    source->sigma = sigma;
    source->omega = omega;
    source->column_size = pow ( source->sigma->length, m );

    // Data
    source->raw = realloc ( source->raw, sizeof ( unsigned long * ) );
    source->raw[0] = calloc ( 
            source->omega->length * source->sigma->length,
            sizeof ( unsigned long )
            );

    return source;
}

int source_update ( unsigned char * prefix, int pos, unsigned char out, source_t * source ){
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
                    source->omega->length * source->sigma->length,
                    sizeof ( unsigned long )
                    );
            if ( source->raw[i] == NULL ){
                return -1;
            }
        }
        source->n = index + 1;
    }

    int min = ( pos < source->m ) ? pos : source->m; 
    int hash_prefix = alphabet_hash ( prefix, min, source->sigma );
    int hash_out = alphabet_hash ( &out, 1, source->omega ) ;

    if ( pos < source->m ){
        int c = hash_prefix;
        for ( int i = 0; i < pow ( source->sigma->length, source->m - pos ); i ++ ){
            source->raw[index][ c * source->column_size + hash_out ] ++;
            c += pow ( source->sigma->length, pos );
        }
    }


    //  Update stats
    source->raw[index][ hash_prefix * source->column_size + hash_out ] ++;
        
    return 0;
}

unsigned char source_generate ( unsigned char * prefix, int pos, source_t * source );

void source_destroy ( source_t * source ){
    for ( int i = 0; i < source->n; i ++ ){
        free ( source->raw[i] );
    }
    free ( source->raw );
    free ( source );
}
