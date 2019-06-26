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
#include "model.h"
#include "stats.h"
#include "source.h"

// Init kseq structure
KSEQ_INIT ( gzFile, gzread );

void usage ( char * name ) {
    fprintf ( stderr, "Usage: %s error_model fastq [fastq ...]\n", name );
}

int main ( int argc, char ** argv ) {
    // Input
    int ploidy;
    char * model_name;
    char * fastq;
    // Error
    FILE * model_fp;
    model_t * model;
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
    model_fp = fopen ( model_name, "r" );
    if ( model_fp == NULL ) {
        fprintf ( stderr, "File %s not found.\n", model_name );
        exit ( EXIT_FAILURE );
    }

    // Init model
    model = model_parse ( model_fp );
    fclose ( model_fp );

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

    // TODO: simulation
    for ( int i = 0; i < ploidy && false; i ++ ){
        while ( kseq_read ( seq[i] ) >= 0 ) {
            int j = 0;
            while ( j < seq[i]->seq.l ){
                generated = stats_generate_read ( &seq[i]->seq.s[j], generated, model->single );
                if ( ! generated->cut ){
                    printf ( ">%s %d\n", seq[i]->name.s, j );
                    printf ( "%s\n", generated->read );
                    printf ( "+\n" );
                    printf ( "%s\n\n", generated->quality );
                }
                // TODO: arbitrary value, this is only for testing
                j += 20;
            }
        }
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i++ ) {
        gzclose ( fp[i] );
        kseq_destroy ( seq[i] );
    }
    if ( generated != NULL ){
        free ( generated->read );
        free ( generated->align );
        free ( generated->quality );
    }
    free ( generated );
    free ( fp );
    free ( seq );
    model_destroy ( model );
    exit ( EXIT_SUCCESS );
}
