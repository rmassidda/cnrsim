/*
 * CNRSIM
 * stats.h
 * Defines the structure containing
 * the statistics of a sequencer.
 *
 * @author Riccardo Massidda
 */

#ifndef STATS
#define STATS
#include <stdbool.h>
#include <stdio.h>
#include "source.h"

enum cigar{
    MAT = 0,
    DEL = 1,
    INS = 2,
    MIS = 3,
};

typedef struct stats_t stats_t;

struct stats_t {
    source_t * alignment;
};

/*
 * Initalize the structure
 *
 * @ret initialized structure, NULL if error
 */
stats_t * stats_init ( );

/*
 * Updates the statistics
 *
 * @param align         alignment string
 * @param alg_len       length of the alignment
 * @param read          readen nucleotides
 * @param read_len      length of the read
 * @param reference     reference sequence
 * @return false if the read is been ignored, true otherwise
 */
bool stats_update ( stats_t * stats, unsigned char * align, int alg_len, char * read, int read_len, char * reference );

/*
 * Frees the memory
 *
 * @param al pointer to the aligner
 */
void stats_destroy ( stats_t * stats );

/*
 * Dumps the known statistics in a file
 *
 * @param file  pointer to the file
 * @param stats pointer to the statistics
 */
void stats_dump ( FILE * file, stats_t * stats );

#endif
