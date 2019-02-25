/*
 * CNRSIM
 * user_variation.h
 * Library that parses a tab-separated file
 * provided by the user containing genomic
 * variants.
 *
 * @author Riccardo Massidda
 */
#ifndef USER_VARIATION
#define USER_VARIATION

#include <uthash.h>
#include <stdbool.h>

#ifndef ALL_N
#define ALL_N 2
#endif

typedef struct variation_t variation_t;
typedef struct variation_set_t variation_set_t;

struct variation_t{
    char * region;
    int pos;
    char * ref;
    char * all[ALL_N];
};

struct entry{
    char * region; // region label
    int position; // first position in the variation array
    UT_hash_handle hh;
};

struct variation_set_t{
    int n; // number of elements
    char * current_region; // last region visited
    int next_variation; // index of the next line
    variation_t ** elements; // variation parsed
    struct entry * region_index;
};

/*
 * Initialize a structure containing
 * the set of user defined variations.
 *
 * @param filename path to the file to be parsed
 * @returns pointer to the initialized set, NULL otherwise
 */
variation_set_t * udv_init ( char * filename );

/*
 * Parse a line of the user defined variations
 * elements must be tab separated.
 * 
 * @param line line containing the variation
 * @returns pointer to the structure, NULL in case of errors
 */
variation_t * udv_parse( char * line );

/*
 * Add a variation to the set
 *
 * @param set pointer to the variation set
 * @param var object containing the variation
 * @returns the position in the array if correctly added, -1 otherwise
 */
int udv_add ( variation_set_t * set, variation_t * var );

/*
 * Sets the position of the reader
 * to the start point of the region.
 *
 * @param set pointer to the variation set
 * @param region label of the desired region
 * @returns false if there isn't such region, true otherwise
 */
bool udv_seek ( variation_set_t * set, char * label );

/*
 * Checks if there are any lines available
 * in the current region.
 *
 * @param set pointer to the variation set
 * @returns true or false
 */
bool udv_next_line(variation_set_t * set);

/*
 * Getter for the next available line.
 *
 * @param set pointer to the variation set
 * @returns pointer to the structure, NULL if there isn't a next line
 */
variation_t * udv_get_line(variation_set_t * set);

/*
 * Free allocated memory
 *
 * @param set pointer to the variation set
 */
void udv_destroy (variation_set_t * set);

#endif
