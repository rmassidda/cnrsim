#include <stdio.h>
#include "user_variation.h"

variation_set_t * udv_init ( char * filename ){
    // File related variables
    FILE * udv_file;
    size_t len = 0;
    ssize_t read = 0;
    char * line = NULL;
    int i = 0;
    // Variation set
    variation_set_t * udv_set;
    // Variation
    variation_t * udv_var;
    // Entries
    struct entry * udv_entry;
    char * tmp_region = NULL;

    // File open
    udv_file = fopen ( filename, "r" );
    if ( udv_file == NULL ){
        return NULL;
    }

    // Allocate set
    udv_set = malloc ( sizeof(variation_set_t) );
    // Initialize values
    udv_set->n = 0;
    udv_set->current_region = NULL;
    udv_set->next_variation = 0;
    udv_set->elements = NULL;
    udv_set->region_index = NULL;


    // Read line
    while ( ( read = getline ( &line, &len, udv_file ) ) != -1 ){
        // Remove new line
        if ( line[read - 1] == '\n'){
            line[read - 1] = '\0';
            read--;
        }
        udv_var = udv_parse( line );
        // Variation correctly parsed    
        if ( udv_var != NULL ){
            // Add to the set
            i = udv_add ( udv_set, udv_var );
            if ( i < 0 ){
                fclose ( udv_file );
                free ( line );
                return NULL;
            }            
            // Is it a new region?
            if ( tmp_region == NULL || strcmp ( tmp_region, udv_set->elements[i]->region ) != 0 ){
                tmp_region = udv_set->elements[i]->region;
                udv_entry = malloc ( sizeof ( struct entry ) );
                udv_entry->region = udv_set->elements[i]->region;        
                udv_entry->position = i;
                HASH_ADD_KEYPTR (
                       hh, // uthash refernce
                       udv_set->region_index, // index pointer
                       udv_entry->region, // id
                       strlen (udv_entry->region), // id length 
                       udv_entry // structure
                       );
            }
        }
    }
    // Cleanup 
    free ( line );
    fclose ( udv_file );
    return udv_set;
}

variation_t * udv_parse( char * line ){
    // Allocate variation
    variation_t * var = malloc ( sizeof ( variation_t ) );
    // Parse region
    char * region = strtok ( line, "\t" );
    if ( region == NULL ){
        return NULL;
    }
    var->region = malloc ( sizeof( char ) * ( strlen ( region ) + 1 ) );
    strcpy ( var->region, region );
    // Parse position
    char * pos = strtok ( NULL, "\t" );
    if ( pos == NULL ){
        return NULL;
    }
    var->pos = atoi ( pos );
    // Parse reference
    char * ref = strtok ( NULL, "\t" );
    if ( ref == NULL ){
        return NULL;
    }
    var->ref = malloc ( sizeof( char ) * ( strlen ( ref ) + 1 ) );
    strcpy ( var->ref, ref );
    // Parse first allele
    char * all1 = strtok ( NULL, "\t" );
    if ( all1 == NULL ){
        return NULL;
    }
    var->all1 = malloc ( sizeof( char ) * ( strlen ( all1 ) + 1 ) );
    strcpy ( var->all1, all1 );
    // Parse second allele
    char * all2 = strtok ( NULL, "\t" );
    if ( all2 == NULL ){
        return NULL;
    }
    var->all2 = malloc ( sizeof( char ) * ( strlen ( all2 ) + 1 ) );
    strcpy ( var->all2, all2 );
    return var;
}

int udv_add ( variation_set_t * set, variation_t * var ){
    int i;
    // Update number of elements
    set->n ++;
    set->elements = realloc ( set->elements, sizeof( variation_t * ) * set->n );
    // Check for errors
    if ( set->elements == NULL ){
        return -1;
    }
    i = set->n - 1;
    set->elements[i] = var;
    return i;
}

bool udv_seek ( variation_set_t * set, char * label ){
    struct entry * result;
    // Search in the Hash Table
    HASH_FIND_STR ( set->region_index, label, result );
    if (result){
        set->current_region = result->region;
        set->next_variation = result->position;       
        return true;
    }
    return false;
}

bool udv_next_line(variation_set_t * set){
    // Reached end of variation set
    if ( set->next_variation < set->n ){
        // The variation refers to a new region
        char * next_region = set->elements[set->next_variation]->region;
        if ( strcmp ( set->current_region, next_region ) == 0){
            return true;
        }
        else{
            return false;
        }
    }
    return false;
}

variation_t * udv_get_line(variation_set_t * set){
    // Returns the variation and increment the counter
    return set->elements[set->next_variation++];
}

void udv_destroy (variation_set_t * set){
    struct entry *udv_entry, *tmp;
    // Free of the entries in the hash table
    HASH_ITER( hh, set->region_index, udv_entry, tmp ){
        HASH_DEL( set->region_index, udv_entry );
        free ( udv_entry );
    }
    // Free of the parsed lines
    for ( int i = 0; i < set->n; i ++ ){
        free ( set->elements[i]->region );
        free ( set->elements[i]->ref );
        free ( set->elements[i]->all1 );
        free ( set->elements[i]->all2 );
        free ( set->elements[i] );
    }
    // Free of the array
    free ( set->elements );
    // Free of the set
    free ( set );
}

void udv_print ( variation_t * var ){
    printf (
            "%s\t%d\t%s\t%s\t%s\n", 
            var->region,
            var->pos,
            var->ref,
            var->all1,
            var->all2 );
}

// Main function to test the library
int main ( int argc, char ** argv ){
    // Variation set
    variation_set_t * set;

    // Initialize set
    set = udv_init ( argv[1] );

    // Test: read all the lines
    while ( set->next_variation < set->n ){
        udv_print ( udv_get_line ( set ) );
    }

    // Read region chr7
    if ( udv_seek ( set, "chr7" ) ){
        while ( udv_next_line ( set ) ){
           udv_print ( udv_get_line ( set ) ); 
        }
    }

    // Read region chr1
    if ( udv_seek ( set, "chr1" ) ){
        while ( udv_next_line ( set ) ){
           udv_print ( udv_get_line ( set ) ); 
        }
    }
    
    // Read region chrY
    if ( udv_seek ( set, "chrY" ) ){
        while ( udv_next_line ( set ) ){
           udv_print ( udv_get_line ( set ) ); 
        }
    }

    // Read region chrX
    if ( udv_seek ( set, "chrX" ) ){
        while ( udv_next_line ( set ) ){
           udv_print ( udv_get_line ( set ) ); 
        }
    }

    udv_destroy ( set );
    return 0;
}
