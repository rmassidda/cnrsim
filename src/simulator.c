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
    int insert_size = 0;
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
    read_t * generated = NULL;
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
        line[read - 1] = '\0';
        // Parse first token
        token = strtok ( line, " " );    

        // Start of new stats
        if ( token[0] == '#' ){
            // Read end
            if ( strcmp ( token, "#single" ) == 0 ){
                curr_end = single;
            }
            else if ( strcmp ( token, "#pair" ) == 0 ){
                curr_end = pair;
            }
            // Insert size
            token = strtok ( NULL, " " );    
            insert_size = atoi ( token );
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
                generated = stats_generate_read ( &seq[i]->seq.s[j], generated, single );
                if ( ! generated->cut ){
                    printf ( ">%s %d\n", seq[i]->name.s, j );
                    printf ( "%s\n", generated->read );
                    printf ( "+\n" );
                    printf ( "%s\n\n", generated->quality );
                }
                // Generate pair
                if ( j + insert_size < seq[i]->seq.l ){
                    generated = stats_generate_read ( &seq[i]->seq.s[j+insert_size], generated, pair );
                    if ( ! generated->cut ){
                        printf ( ">%s %d\n", seq[i]->name.s, j + insert_size );
                        printf ( "%s\n", generated->read );
                        printf ( "+\n" );
                        printf ( "%s\n\n", generated->quality );
                    }
                }
            }
        }
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i++ ) {
        gzclose ( fp[i] );
        kseq_destroy ( seq[i] );
    }
    fclose ( model );
    free ( generated->read );
    free ( generated->align );
    free ( generated->quality );
    free ( generated );
    free ( fp );
    free ( seq );
    free ( line );
    stats_destroy ( single );
    stats_destroy ( pair );
    exit ( EXIT_SUCCESS );
}
