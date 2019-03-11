/*
 * CNRSIM
 * allele.h
 * Defines the structure containing
 * the sequence in an allele
 *
 * @author Riccardo Massidda
 */
#ifndef VARIATOR_H
#define VARIATOR_H

typedef struct allele_t allele_t;

struct allele_t {
    char * sequence; // pointer to the mutated sequence
    char * alignment; // pointer to the alignment of the sequence
    long int buffer_size; // size of the sequence
    long int pos; // current position
    long int off; // offset from reference position
};

/*
 * Initialize an allele
 *
 * @param       size  size of the reference, slightly bigger to allow variations
 * @param       allele        allele to be initialized
 * @returns     the initialized structure
 */
allele_t * allele_init ( long int size, allele_t * allele );

/*
 * Deallocates an allele
 *
 * @param       allele  allele to be deallocated
 */
void allele_destroy ( allele_t * allele );

#endif
