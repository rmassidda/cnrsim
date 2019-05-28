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

model_t * model_init ( int max_repetition ){
    model_t * model = malloc ( sizeof ( model_t ) );
    if ( model == NULL ) return model;

    model->single = stats_init ();
    model->pair = stats_init ();
    model->insert_size = 0;
    model->_read_number = 0;
    model->amplification = source_init ( max_repetition, max_repetition, 1, 0 );

    return model;
}

model_t * model_parse ( FILE * file, model_t * model ){
    // Parsing
    size_t len = 0;
    ssize_t read = 0;
    char * line = NULL;
    int n_line = 0;
    double w;
    stats_t * curr_end = NULL;
    source_t * curr_source = NULL;
    char * token;

    // Model Parsing
    curr_end = model->single;
    while ( ( read = getline ( &line, &len, file ) ) != -1 ) {
        // Remove new line
        line[read - 1] = '\0';
        // Parse first token
        token = strtok ( line, " " );    

        // Start of new stats
        if ( token[0] == '#' ){
            // Read end
            if ( strcmp ( token, "#single" ) == 0 ){
                curr_end = model->single;
            }
            else if ( strcmp ( token, "#pair" ) == 0 ){
                curr_end = model->pair;
                // Insert size
                token = strtok ( NULL, " " );    
                model->insert_size = atoi ( token );
            }
            else if ( strcmp ( token, "#amplification" ) == 0 ){
                curr_end = NULL;
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
            else{
                printf ( "%s, not parsable.\n", token );
                exit ( EXIT_FAILURE );
            }
            // Get length
            token = strtok ( NULL, " " );    
            curr_source->n = atoi ( token );
            // Allocate matrix
            curr_source->normalized = malloc ( sizeof ( unsigned long * ) * curr_source->n );
            curr_source->raw = malloc ( sizeof ( double * ) * curr_source->n );
            for ( int i = 0; i < curr_source->n; i ++ ) {
                curr_source->raw[i] = NULL;
                curr_source->normalized[i] = malloc ( 
                        curr_source->omega * curr_source->prefix * sizeof ( double )
                        );
            }
            n_line = 0;
        }
        // Load data
        else{
            w = atof ( token );
            for ( int i = 0; i < curr_source->prefix; i ++ ) {
                for ( int j = 0; j < curr_source->omega; j ++ ) {
                    curr_source->normalized[n_line][ i * curr_source->omega + j] = w;
                    token = strtok ( NULL, " " );    
                    if ( token != NULL ) {
                        w = atof ( token );
                    }
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
    free ( model );
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
    fprintf ( file, "#amplification\n" );
    source_dump ( file, "tandem", model->amplification );
}
