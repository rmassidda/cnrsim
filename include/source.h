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
typedef struct alphabet_t alphabet_t;

struct alphabet_t {
    unsigned char * symbols;
    int length;
};

struct source_t {
    int m; // memory of the source
    int k; // sampling finite precision
    int n; // number of matrixes
    unsigned long ** raw; // data
    alphabet_t * sigma; // input alphabet
    alphabet_t * omega; // output alphabet
    int column_size; // column
};

/*
 * Initialize an alphabet
 *
 * @param       symbols
 * @param       length
 * @returns     the initialized alphabet
 */
alphabet_t * alphabet_init ( unsigned char * symbols, int length );

/*
 * Contiguos hash
 *
 * @param       word
 * @param       length
 * @param       alphabet
 * @returns     hash value
 */
int alphabet_hash ( unsigned char * word, int length, alphabet_t * alphabet );

/*
 * Initialize the source
 *
 * @param       sigma   input alphabet
 * @param       omega   output alphabet
 * @param       m       finite memory
 * @param       k       sampling finite precision
 * @returns     the initialized structure
 */
source_t * source_init ( alphabet_t * sigma, alphabet_t * omega, int m, int k );

/*
 * Updates the data with a new example
 *
 * @param       prefix  string containing the prefix
 * @param       out     output character
 * @param       pos     position of the example
 * @param       source  source to be updated
 * @returns     0 on success, -1 otherwise 
 */
int source_update ( unsigned char * prefix, int pos, unsigned char out, source_t * source ); 

/*
 * Generates a character given a prefix
 * and a position
 *
 * @param       prefix  string containing the prefix
 * @param       pos     position of the example
 * @param       source  source to be used
 * @returns     output character given the learned probabilities
 */
unsigned char source_generate ( unsigned char * prefix, int pos, source_t * source );

/*
 * Deallocates a source
 *
 * @param       source  source to be deallocated
 */
void source_destroy ( source_t * source );

#endif
