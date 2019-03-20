/*
 * CNRSIM
 * source.h
 * Finite memory sequence generator
 *
 * @author Riccardo Massidda
 */
#ifndef SOURCE_H
#define SOURCE_H

typedef struct source_t source_t;

struct source_t {
    int m; // memory of the source
    int k; // sampling finite precision
    int n; // number of matrixes
    unsigned long ** raw; // data
    char * sigma; // input alphabet
    int sigma_size; // input alphabet size
    char * omega; // output alphabet
    int omega_size; // output alphabet size
};

/*
 * Initialize the source
 *
 * @param       sigma   input alphabet
 * @param       omega   output alphabet
 * @param       m       finite memory
 * @param       k       sampling finite precision
 * @returns     the initialized structure
 */
source_t * source_init ( char * sigma, char * omega, int m, int k );

/*
 * Updates the data with a new example
 *
 * @param       prefix  string containing the prefix
 * @param       out     output character
 * @param       pos     position of the example
 * @param       source  source to be updated
 * @returns     0 on success, -1 otherwise 
 */
int source_update ( char * prefix, int pos, char out, source_t * source ); 

/*
 * Generates a character given a prefix
 * and a position
 *
 * @param       prefix  string containing the prefix
 * @param       pos     position of the example
 * @param       source  source to be used
 * @returns     output character given the learned probabilities
 */
char source_generate ( char * prefix, int pos, source_t * source );

/*
 * Deallocates a source
 *
 * @param       source  source to be deallocated
 */
void source_destroy ( source_t * source );

#endif
