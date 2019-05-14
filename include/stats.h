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
typedef struct read_t read_t;

struct stats_t {
    source_t * alignment;
    source_t * mismatch;
    source_t * quality;
    source_t * distribution;
};

struct read_t {
    unsigned char * align;
    int alg_len;
    char * read;
    unsigned char * quality;
    bool cut;
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
 * Generates a read using the internal statistics
 *
 * @param ref   reference sequence
 * @param read  pointer to the read to be filled
 * @param stats pointer to the statistics
 * @returns     pointer to the generated structure
 */
read_t * stats_generate_read ( char * ref, read_t * read, stats_t * stats );

/*
 * Frees the memory
 *
 * @param stats pointer to the statistics
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
