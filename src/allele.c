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
    // It's necessary to clean the memory
    memset ( allele->sequence, 0, sizeof ( char ) *allele->buffer_size );
    memset ( allele->alignment, 0, sizeof ( char ) *allele->buffer_size );
    allele->pos = 0;
    allele->off = 0;
    return allele;
}

void allele_destroy ( allele_t * allele ){
    free ( allele->sequence );
    free ( allele->alignment );
    free ( allele );
}
