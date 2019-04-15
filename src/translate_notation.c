/*
 * CNRSIM
 * translate_notation.c
 * Library that parses a tab-separated file
 * provided by the user containing the
 * notation dictionary.
 *
 * @author Riccardo Massidda
 */
#include <stdio.h>
#include "translate_notation.h"

region_index_t * tr_init ( char * filename ) {
    // File related variables
    FILE * file;
    size_t len = 0;
    ssize_t read = 0;
    char * line = NULL;
    // Index
    region_index_t * index = NULL;
    char * reg = NULL;
    char * translation = NULL;

    // File open
    file = fopen ( filename, "r" );
    if ( file == NULL ) {
        return NULL;
    }

    // Read line
    while ( ( read = getline ( &line, &len, file ) ) != -1 ) {
        // Remove new line
        if ( line[read - 1] == '\n' ) {
            line[read - 1] = '\0';
            read--;
        }

        reg = strtok ( line, "\t" );
        translation = strtok ( NULL, "\t" );
        if ( reg != NULL && translation != NULL ) {
            // Update dictionary
            region_index_t * entry = malloc ( sizeof ( region_index_t ) );
            entry->reg = malloc ( sizeof ( char ) * ( strlen ( reg ) + 1 ) );
            entry->alt = malloc ( sizeof ( char ) * ( strlen ( translation ) + 1 ) );
            strcpy ( entry->reg, reg );
            strcpy ( entry->alt, translation );
            HASH_ADD_KEYPTR (
                hh,
                index,
                entry->reg,
                strlen ( entry->reg ),
                entry
            );
        }
    }
    // Cleanup
    free ( line );
    fclose ( file );
    return index;
}


char * tr_translate ( region_index_t * index, char * label ) {
    region_index_t * result;
    // Search in the Hash Table
    HASH_FIND_STR ( index, label, result );
    if ( result ) {
        return result->alt;
    }
    return NULL;
}

void tr_destroy ( region_index_t * index ) {
    region_index_t * entry;
    region_index_t * tmp;
    // Free of the entries in the hash table
    HASH_ITER ( hh, index, entry, tmp ) {
        HASH_DEL ( index, entry );
        free ( entry->reg );
        free ( entry->alt );
        free ( entry );
    }
    // Free of the set
    free ( index );
}
