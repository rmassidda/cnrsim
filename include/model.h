/*
 * CNRSIM
 * model.h
 * Defines the structure containing
 * the statistics of a sequencer.
 *
 * @author Riccardo Massidda
 */

#ifndef MODEL
#define MODEL
#include <stdbool.h>
#include <stdio.h>
#include "source.h"
#include "stats.h"

typedef struct model_t model_t;

struct model_t {
    stats_t * single;
    stats_t * pair;
    source_t * amplification;
    int max_insert_size;
    source_t * insert_size;
    source_t * orientation;
};

/*
 * Initalize the model
 *
 * @ret initialized structure, NULL if error
 */
model_t * model_init ( int max_repetition, int max_insert_size, int size_granularity );

/*
 * Initalize the model from a file
 *
 * @ret initialized structure, NULL if error
 */
model_t * model_parse ( FILE * file );

/*
 * Frees the memory
 *
 * @param model pointer to the statistics
 */
void model_destroy ( model_t * model );

/*
 * Dumps the known statistics in a file
 *
 * @param file  pointer to the file
 * @param model pointer to the statistics
 */
void model_dump ( FILE * file, model_t * model );

#endif
