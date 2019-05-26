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
    int insert_size;
    int _read_number;
};

/*
 * Initalize the model
 *
 * @ret initialized structure, NULL if error
 */
model_t * model_init ( int max_repetition );

/*
 * Initalize the model from a file
 *
 * @ret initialized structure, NULL if error
 */
model_t * model_parse ( FILE * file, model_t * model );

void update_insert_size ( int value, model_t * model );

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
