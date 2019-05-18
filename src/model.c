/*
 * CNRSIM
 * model.c
 * Defines the structure containing
 * the statistics of a sequencer.
 *
 * @author Riccardo Massidda
 */

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "stats.h"
#include "model.h"

model_t * model_init (){
    model_t * model = malloc ( sizeof ( model_t ) );
    if ( model == NULL ) return model;

    model->single = stats_init ();
    model->pair = stats_init ();
    model->insert_size = 0;
    model->_read_number = 0;
    model->amplification = NULL;

    return model;
}

model_t * model_parse ( FILE * file, model_t * model ){
    return model;
}

void model_destroy ( model_t * model ){
    stats_destroy ( model->single );
    stats_destroy ( model->pair );
    source_destroy ( model->amplification );
}

void update_insert_size ( int value, model_t * model ){
    model->insert_size += value;
    model->_read_number ++;
}

void model_dump ( FILE * file, model_t * model ){
    if ( model->_read_number != 0 ){
        model->insert_size = model->insert_size / model->_read_number;
    }
    fprintf ( file, "#single 0\n" );
    stats_dump ( file, model->single );
    fprintf ( file, "#pair %d\n", model->insert_size );
    stats_dump ( file, model->pair );
}
