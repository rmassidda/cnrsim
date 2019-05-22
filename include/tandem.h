/*
 * CNRSIM
 * tandem.h
 * Defines the structure containing
 * a tandem repeat inside the reference
 *
 * @author Riccardo Massidda
 */
#ifndef TANDEM_H
#define TANDEM_H

typedef struct tandem_t tandem_t;
typedef struct tandem_set_t tandem_set_t;

struct tandem_t {
    unsigned int pos; // position inside of the reference
    unsigned short pat; // size of the pattern
    unsigned short rep; // number of repetitions
};

struct tandem_set_t {
    tandem_t * set; // set of tandem repeats in a sequence
    int i; // index of the current tandem
    int n; // number of tr
    int size; // size of allocated memory
    int max_motif; // size of the maximum motif analyzed
    int max_repetition; // number of maximum allowed repetitions
};


/*
 * Initialize a set
 *
 * @param       size            size of the reference
 * @param       max_motif       size of the maximum motif
 * @param       set             set to be initialized
 * @returns     the initialized structure
 */
tandem_set_t * tandem_set_init ( int length, int max_motif, int max_repetition, tandem_set_t * set );

/*
 * Populates a set analyzing a reference sequence
 *
 * @param       reference       reference sequence
 * @param       length          size of the reference
 * @param       set             empty set
 * @returns     set containing the repetitions
 */
tandem_set_t * tandem_set_analyze ( char * reference, int length, tandem_set_t * set );


/*
 * Destroy a set
 *
 * @param       set  set to be destroyed
 */
void tandem_set_destroy ( tandem_set_t * set );

#endif
