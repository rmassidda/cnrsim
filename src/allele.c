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
        allele->max_ref = 0;
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
            if ( allele->alignment[allele->alg] == 'I' ){
                allele->off -= direction;
            }
            else if ( allele->alignment[allele->alg] == 'D' ){
                allele->off += direction;
            }
        // If there is a 'D' or a '=' it corresponds
        // to a new position on the reference
        }while ( allele->alignment[allele->alg] == 'I' );
        allele->ref += direction;
    }
    // Update of the position
    allele->pos = allele->ref - allele->off;
    return allele->pos;
}

int allele_variation ( char * ref, char * alt, allele_t * allele ){
    int ref_len = strlen ( ref );
    int alt_len = strlen ( alt );
    char * all = &allele->sequence[allele->pos];
    int all_len = strlen ( all );
    int offset = ref_len - alt_len;
    int min = ( offset < 0 ) ? ref_len : alt_len;
    int max = ( offset < 0 ) ? alt_len : ref_len;
    char c = ( offset < 0 ) ? 'I' : 'D';

    // Variation coerency check
    if ( all_len <= ref_len ){
        if ( strncmp ( all, ref, all_len ) != 0 ){
            return -1;
        }
    }
    else{
        if ( strncmp ( all, ref, ref_len ) != 0 ){
            return -1;
        }
        // Avoid INDEL inside the allele
        else if ( alt_len != ref_len ){
            return -2;
        }
    }

    // Check if the position stil exists
    for ( int i = 0; i < max; i ++ ){
        if ( allele->alignment[allele->alg + i] == 'D' ){
            return -3;
        }
    }

    // Update reference limit
    if ( allele->ref + ref_len > allele->max_ref ){
        allele->max_ref = allele->ref + ref_len;
    }

    // Application
    memcpy (
            &allele->sequence[allele->pos],
            alt,
            alt_len
           );
    if ( alt_len < ref_len ){
        memset (
                &allele->sequence[allele->pos+alt_len],
                '\0',
                ref_len - alt_len
               );
    }
    for ( int i = 0; i < min; i ++ ){
        if ( alt[i] == ref[i] ){
            allele->alignment[allele->alg + i] = '=';
        }
        else { 
            allele->alignment[allele->alg + i] = '!';
        }
    }
    memset ( 
            &allele->alignment[allele->alg + min],
            c,
            abs ( offset )
           );
    return 0;
}

void allele_destroy ( allele_t * allele ){
    free ( allele->sequence );
    free ( allele->alignment );
    free ( allele );
}
