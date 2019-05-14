/*
 * CNRSIM
 * variator.c
 * Generatos two alleles from a FASTA reference
 * and a VCF file containing known variations.
 *
 * @author Riccardo Massidda
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include "stats.h"
#include "source.h"

// Init kseq structure
KSEQ_INIT ( gzFile, gzread );

void usage ( char * name ) {
    fprintf ( stderr, "Usage: %s error_model fastq [fastq ...]\n", name );
}

int main ( int argc, char ** argv ) {
    // Parsing
    size_t len = 0;
    ssize_t read = 0;
    char * line = NULL;
    int n_line = 0;
    int status = -1;
    double w;
    stats_t * curr_end = NULL;
    source_t * curr_source = NULL;
    char * token;
    // Input
    int ploidy;
    char * model_name;
    char * fastq;
    // Error
    FILE * model;
    stats_t * single;
    stats_t * pair;
    // FASTA
    gzFile * fp;
    kseq_t ** seq;

    // Non optional arguments
    if ( argc - optind < 2 ) {
        usage ( argv[0] );
        exit ( EXIT_FAILURE );
    }

    // Error Model
    model_name = argv[optind++];
    // File open
    model = fopen ( model_name, "r" );
    if ( model == NULL ) {
        fprintf ( stderr, "File %s not found.\n", model_name );
        exit ( EXIT_FAILURE );
    }
    // Init stats
    single = stats_init();
    pair = stats_init();

    // Model Parsing
    curr_end = single;
    while ( ( read = getline ( &line, &len, model ) ) != -1 ) {
        // Remove new line
        if ( line[read - 1] == '\n' ) {
            line[read - 1] = '\0';
            read--;
        }
        // Start of new source
        if ( line[0] == '@' ){
            // Update parse status
            status ++;
            switch ( status ) {
                case 0: curr_source = curr_end->alignment; break;
                case 1: curr_source = curr_end->mismatch; break;
                case 2: curr_source = curr_end->quality; break;
                case 3: curr_source = curr_end->distribution; break;
                case 4: curr_end = pair; curr_source = curr_end->alignment; status = 0; break;
                default: fprintf ( stderr, "Unexpected source\n" ); exit ( EXIT_FAILURE );
            }
            // Get length
            token = strtok ( line, " " );    
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
            token = strtok ( line, " " );    
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


    // Input sequences
    ploidy = argc - optind;
    fp = malloc ( sizeof ( gzFile ) * ploidy );
    seq = malloc ( sizeof ( kseq_t * ) * ploidy );

    // Open input files
    for ( int i = 0; i < ploidy; i ++ ){
        fastq = argv[optind++];
        fp[i] = gzopen ( fastq, "r" );
        if ( fp[i] == NULL ) {
            fprintf ( stderr, "File %s not found.\n", fastq );
            exit ( EXIT_FAILURE );
        }
        seq[i] = kseq_init ( fp[i] );
    }

    for ( int i = 0; i < ploidy; i ++ ){
        while ( kseq_read ( seq[i] ) >= 0 ) {
            for ( int j = 0; j < seq[i]->seq.l; j ++ ){
            }
        }
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i++ ) {
        gzclose ( fp[i] );
        kseq_destroy ( seq[i] );
    }
    fclose ( model );
    free ( fp );
    free ( seq );
    free ( line );
    stats_destroy ( single );
    stats_destroy ( pair );
    exit ( EXIT_SUCCESS );
}
