/*
 * CNRSIM
 * translate_notation.h
 * Library that parses a tab-separated file
 * provided by the user containing the
 * notation dictionary.
 *
 * @author Riccardo Massidda
 */
#ifndef TRANSLATE_NOTATION
#define TRANSLATE_NOTATION

#include <uthash.h>
#include <stdbool.h>

typedef struct region_index_t region_index_t;

struct region_index_t {
    char * reg; // region label
    char * alt; // alternate notation
    UT_hash_handle hh;
};

/*
 * Initialize a structure containing
 * the dictionary.
 *
 * @param filename path to the file to be parsed
 * @returns pointer to the initialized dictionary, NULL otherwise
 */
region_index_t * tr_init ( char * filename );

/*
 * Translates the notation
 *
 * @param index pointer to the index
 * @param label string that contains the notation to translate
 * @returns string containing alternate notation, NULL if there isn't any
 */
char * tr_translate ( region_index_t * index, char * label );

/*
 * Free allocated memory
 *
* @param index pointer to the index
 */
void tr_destroy ( region_index_t * index );

#endif
