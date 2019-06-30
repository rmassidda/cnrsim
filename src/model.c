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

model_t * model_init ( int max_motif, int max_repetition, int max_insert_size, int size_granularity ){
    model_t * model = malloc ( sizeof ( model_t ) );
    if ( model == NULL ) return model;

    model->max_motif = max_motif;
    model->max_repetition = max_repetition;
    model->size_granularity = size_granularity;
    model->max_insert_size = max_insert_size;
    model->single = stats_init ();
    model->pair = stats_init ();
    model->amplification = source_init ( max_repetition, max_repetition, 1, 0 );
    model->insert_size = source_init ( 1, size_granularity, 0, 0 );
    model->orientation = source_init ( 1, 4, 0, 0 );

    return model;
}

model_t * model_parse ( FILE * file ){
    // Read line
    char * line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    int n_line = 0;
    char * token;
    // Pointer to the working data
    model_t * model = NULL;
    stats_t * curr_end = NULL;
    source_t * curr_source = NULL;
    int max_repetition = 0;
    int max_insert_size = 0;
    int max_motif = 0;
    int size_granularity = 0;

    // Model Parsing
    while ( ( read = getline ( &line, &len, file ) ) != -1 ) {
        // Remove new line
        line[read - 1] = '\0';
        // Parse first token
        token = strtok ( line, " " );    
        if ( token == NULL ) {
          continue;
        }

        if ( token[0] == '$' ){
            strtok ( NULL, " " );
            if ( strcmp ( token, "$max_repetition" ) == 0 ){
                max_repetition = atoi ( token );
            }
            else if ( strcmp ( token, "$max_insert_size" ) == 0 ){
                max_insert_size = atoi ( token );
            }
            else if ( strcmp ( token, "$size_granularity" ) == 0 ){
                size_granularity = atoi ( token );
            }
            else if ( strcmp ( token, "$max_motif" ) == 0 ){
                max_motif = atoi ( token );
            }
            else{
                fprintf ( stderr, "%s not parsable.\n", token );
                exit ( EXIT_FAILURE );
            }
        }
        else if ( token[0] == '#' ){
            if ( model == NULL ) {
                model = model_init ( max_motif, max_repetition, max_insert_size, size_granularity );
            }
            if ( strcmp ( token, "#single" ) == 0 ){
                curr_end = model->single;
            }
            else if ( strcmp ( token, "#pair" ) == 0 ){
                curr_end = model->pair;
            }
            else if ( strcmp ( token, "#amplification" ) == 0 ){
                curr_end = NULL;
            }
            else{
                fprintf ( stderr, "%s not parsable.\n", token );
                exit ( EXIT_FAILURE );
            }
        }
        else if ( token[0] == '@' ){
            // Update parse status
            if ( strcmp ( token, "@alignment" ) == 0 ) {
                curr_source = curr_end->alignment;
            }
            else if ( strcmp ( token, "@mismatch" ) == 0 ){
                curr_source = curr_end->mismatch;
            }
            else if ( strcmp ( token, "@quality" ) == 0 ){
                curr_source = curr_end->quality;
            }
            else if ( strcmp ( token, "@distribution" ) == 0 ){
                curr_source = curr_end->distribution;
            }
            else if ( strcmp ( token, "@tandem" ) == 0 ){
                curr_source = model->amplification;
            }
            else if ( strcmp ( token, "@insert_size" ) == 0 ){
                curr_source = model->insert_size;
            }
            else if ( strcmp ( token, "@orientation" ) == 0 ){
                curr_source = model->orientation;
            }
            else{
                fprintf ( stderr, "%s not parsable.\n", token );
                exit ( EXIT_FAILURE );
            }
            // Get length
            token = strtok ( NULL, " " );    
            curr_source->n = atoi ( token );
            // Pre-allocate matrixes
            curr_source->normalized = malloc ( sizeof ( double * ) * curr_source->n );
            curr_source->raw = malloc ( sizeof ( unsigned long * ) * curr_source->n );
            for ( int i = 0; i < curr_source->n; i ++ ) {
                curr_source->normalized[i] = NULL;
                curr_source->raw[i] = malloc ( 
                        curr_source->omega * curr_source->prefix * sizeof ( unsigned long )
                        );
            }
            n_line = 0;
        }
        // Load data
        else{
            for ( int i = 0; i < curr_source->prefix; i ++ ) {
                for ( int j = 0; j < curr_source->omega; j ++ ) {
                    curr_source->raw[n_line][ i * curr_source->omega + j] = strtoul ( token, NULL, 10 );
                    token = strtok ( NULL, " " );    
                }
            }
            n_line ++;
        }
    }

    free ( line );

    return model;
}

void model_destroy ( model_t * model ){
    if ( model == NULL ){
        return;
    }
    stats_destroy ( model->single );
    stats_destroy ( model->pair );
    source_destroy ( model->amplification );
    source_destroy ( model->insert_size );
    source_destroy ( model->orientation );
    free ( model );
}

void model_dump ( FILE * file, model_t * model ){
    fprintf ( file, "$max_insert_size %d\n", model->max_insert_size );
    fprintf ( file, "$max_repetition %d\n", model->max_repetition );
    fprintf ( file, "$max_motif %d\n", model->max_motif );
    fprintf ( file, "$size_granularity %d\n", model->size_granularity );
    fprintf ( file, "#single\n" );
    stats_dump ( file, model->single );
    fprintf ( file, "#pair\n" );
    stats_dump ( file, model->pair );
    source_dump ( file, "insert_size", model->insert_size );
    source_dump ( file, "orientation", model->orientation );
    fprintf ( file, "#amplification\n" );
    source_dump ( file, "tandem", model->amplification );
}
