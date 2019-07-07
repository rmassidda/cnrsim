/*
 * CNRSIM
 * allele.c
 * Defines the structure containing
 * the sequence in an allele
 *
 * @author Riccardo Massidda
 */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "allele.h"

allele_t * allele_init ( long int size, allele_t * allele ) {
    // First initialization
    if ( allele == NULL ) {
        allele = malloc ( sizeof ( allele_t ) );
        allele->sequence = NULL;
        allele->alignment = NULL;
    }
    // Update of internal values
    allele->buffer_size = floor ( size * 1.5 );
    allele->sequence = realloc ( allele->sequence, ( sizeof ( char ) ) * allele->buffer_size );
    allele->alignment = realloc ( allele->alignment, ( sizeof ( char ) ) * allele->buffer_size );
    // Initial conditions
    allele->pos = 0;
    allele->ref = 0;
    allele->alg = 0;
    return allele;
}

allele_t * allele_point ( long int size, char * sequence, char * alignment, allele_t * allele ) {
    if ( allele == NULL ) {
        allele = malloc ( sizeof ( allele_t ) );
    }
    allele->buffer_size = size;
    allele->sequence = sequence;
    allele->alignment = alignment;
    allele->pos = 0;
    allele->ref = 0;
    allele->alg = 0;
    return allele;
}

int allele_seek ( int position, bool from_reference, allele_t * allele ) {
    long int * p;
    char target;
    if ( allele->alignment == NULL ) {
        allele->pos = position;
        allele->ref = position;
        return allele->pos;
    }

    // Set which position is to be monitored
    p = ( from_reference ) ? &(allele->ref) : &(allele->pos);
    // Set which character has to be ignored
    target = ( from_reference ) ? 'I' : 'D';

    // Up to the end of the alignment
    while ( allele->alg < allele->buffer_size ) {
      // If the desired position is found and
      // there are no more corresponding ones
      if ( *p == position && allele->alignment[allele->alg] != target ) {
        return ( from_reference ) ? allele->pos : allele->ref;
      }
      // Increase the pointers
      switch ( allele->alignment[allele->alg] ) {
        case 'I': allele->pos ++;
        case 'D': allele->ref ++;
        default: allele->pos ++; allele->ref ++;
      }
      allele->alg ++;
    }
    // Position not found
    return ( from_reference ) ? allele->pos : allele->ref;
}

int allele_variation ( char * ref, char * alt, allele_t * allele ) {
    int ref_len = strlen ( ref );
    int alt_len = strlen ( alt );
    int offset = ref_len - alt_len;
    int min = ( offset < 0 ) ? ref_len : alt_len;
    int max = ( offset < 0 ) ? alt_len : ref_len;
    char c = ( offset < 0 ) ? 'I' : 'D';

    // Application
    memcpy (
        &allele->sequence[allele->pos],
        alt,
        alt_len
    );

    // Alignment
    for ( int i = 0; i < min; i ++ ) {
        if ( alt[i] == ref[i] ) {
            allele->alignment[allele->alg + i] = '=';
        } else {
            allele->alignment[allele->alg + i] = 'X';
        }
    }
    memset (
        &allele->alignment[allele->alg + min],
        c,
        abs ( offset )
    );

    //  Update index
    allele->ref += ref_len;
    allele->pos += alt_len;
    allele->alg += max;

    return 0;
}

void allele_destroy ( allele_t * allele ) {
    free ( allele->sequence );
    free ( allele->alignment );
    free ( allele );
}
