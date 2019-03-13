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
        memset ( allele->sequence, 0, sizeof ( char ) *allele->buffer_size );
        memset ( allele->alignment, '=', sizeof ( char ) *allele->buffer_size );
        allele->pos = 0;
        allele->off = 0;
        allele->ref = 0;
        allele->alg = 0;
        return allele;
}

allele_t * allele_point ( long int size, char * sequence, char * alignment, allele_t * allele ){
    if ( allele == NULL ){
        allele = malloc ( sizeof ( allele_t ) );
    }
    allele->buffer_size = size;
    allele->sequence = sequence;
    allele->alignment = alignment;
    allele->pos = 0;
    allele->off = 0;
    allele->ref = 0;
    allele->alg = 0;
    return allele;
}

int allele_seek ( int position, allele_t * allele ){
    if ( allele->alignment == NULL ){
        allele->pos = position;
        allele->ref = position;
        return allele->pos;
    }
    while ( position > allele->ref ){
        // The structure contains the offset
        // computed up to the allele->ref position
        // in the reference
        do{
            // Next position in the alignment
            allele->alg ++;
            // Offset update
            if ( allele->alignment[allele->alg] == 'i' ){
                allele->off --;
            }
            else if ( allele->alignment[allele->alg] == 'd' ){
                allele->off ++;
            }
        // If there is a 'd' or a '=' it corresponds
        // to a new position on the reference
        }while ( allele->alignment[allele->alg] == 'i' );
        allele->ref ++;
    }
    // Update of the position
    allele->pos = allele->ref - allele->off;
    return allele->pos;
}


void allele_destroy ( allele_t * allele ){
    free ( allele->sequence );
    free ( allele->alignment );
    free ( allele );
}
