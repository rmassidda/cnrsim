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
    int n; // number of matrixes
    unsigned long ** raw; // data
    double ** normalized; // normalized data
    int m; // memory
    int sigma; // input alphabet
    int omega; // output alphabet
};

/*
 * Initialize the source
 *
 * @param       sigma   input alphabet
 * @param       omega   output alphabet
 * @param       m       memory
 * @returns     the initialized structure
 */
source_t * source_init ( int sigma, int omega, int m );

/*
 * Updates the data with a new example
 *
 * @param       in      input prefix
 * @param       out     output character
 * @param       len     lenghth of the prefix
 * @param       pos     position of the example
 * @param       source  source to be updated
 * @returns     0 on success, -1 otherwise 
 */
int source_update ( unsigned char * in, int len, int pos, unsigned char out, source_t * source ); 

/*
 * Generates a character given a prefix
 * and a position
 *
 * @param       prefix  string containing the prefix
 * @param       len     lenghth of the prefix
 * @param       pos     position of the example
 * @param       source  source to be used
 * @returns     output character given the learned probabilities
 */
unsigned char source_generate ( unsigned char * in, int len, int pos, source_t * source );

/*
 * Deallocates a source
 *
 * @param       source  source to be deallocated
 */
void source_destroy ( source_t * source );

#endif
