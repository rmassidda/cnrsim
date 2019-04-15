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
#include "allele.h"

typedef struct stats_t stats_t;

struct stats_t {
    source_t * alignment;
    source_t * mismatch;
    source_t * quality;
    source_t * distribution;
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
 * @param ref           reference sequence
 * @param quality       quality score
 */
void stats_update ( unsigned char * align, int alg_len, char * read, char * ref, unsigned char * quality, stats_t * stats );

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
