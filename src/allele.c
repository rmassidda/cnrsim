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
        memset ( allele->sequence, '\0', sizeof ( char ) *allele->buffer_size );
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
    int direction = ( position < allele->ref ) ? -1 : 1;
    while ( position != allele->ref ){
        // The structure contains the offset
        // computed up to the allele->ref position
        // in the reference
        do{
            // Next position in the alignment
            allele->alg += direction;
            // Offset update
            if ( allele->alignment[allele->alg] == 'i' ){
                allele->off -= direction;
            }
            else if ( allele->alignment[allele->alg] == 'd' ){
                allele->off += direction;
            }
        // If there is a 'd' or a '=' it corresponds
        // to a new position on the reference
        }while ( allele->alignment[allele->alg] == 'i' );
        allele->ref += direction;
    }
    // Update of the position
    allele->pos = allele->ref - allele->off;
    return allele->pos;
}

int allele_variation ( char * ref, char * alt, allele_t * allele ){
    int ref_len = strlen ( ref );
    int alt_len = strlen ( alt );
    int offset = ref_len - alt_len;
    int min = ( offset < 0 ) ? ref_len : alt_len;
    char c = ( offset < 0 ) ? 'i' : 'd';
    memcpy (
            &allele->sequence[allele->pos],
            alt,
            alt_len
           );
    memset ( 
            &allele->alignment[allele->alg],
            '=',
            min
           );
    memset ( 
            &allele->alignment[allele->alg + min],
            c,
            abs ( offset )
           );
    return allele_seek ( allele->ref + ref_len, allele );
}

void allele_dump ( allele_t * allele ){
    int i = 0;
    for ( i = 0; i < allele->ref; i ++ ) printf ( " " );
    printf ( "R\n" );
    for ( i = 0; i < allele->pos; i ++ ) printf ( " " );
    printf ( "P\n" );
    for ( i = 0; i < allele->alg; i ++ ) printf ( " " );
    printf ( "A\n" );
    printf ( "%s\n", allele->sequence );
    printf ( "%.*s\n", (int)allele->alg + 1, allele->alignment );
}

void allele_destroy ( allele_t * allele ){
    free ( allele->sequence );
    free ( allele->alignment );
    free ( allele );
}
