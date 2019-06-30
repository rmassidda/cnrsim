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
#include "tandem.h"

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
    char * amplified_seq = NULL;
    tandem_set_t * tandem = NULL;

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
    for ( int i = 0; i < ploidy; i ++ ){
        while ( kseq_read ( seq[i] ) >= 0 ) {
            // Analysis of the repetitions in the original sequence
            tandem = tandem_set_init ( seq[i]->seq.l, model->max_motif, model->max_repetition, tandem );
            tandem = tandem_set_analyze ( seq[i]->seq.s, seq[i]->seq.l, tandem );
            amplified_seq = realloc ( amplified_seq, sizeof ( char ) * ( seq[i]->seq.l * 2 ) );

            // Index for the FASTA sequence
            int seq_p = 0;
            // Index for the amplified sequence
            int aseq_p = 0;
            // Statistics
            int amp = 0;
            int deamp = 0;
            for ( int t = 0; t < tandem->n; t ++ ) {
              int gap = tandem->set[t].pos - seq_p;
              // Copy of the nucleotides between different tandems
              if ( gap > 0 ) {
                memcpy (
                    &amplified_seq[aseq_p],
                    &seq[i]->seq.s[seq_p],
                    sizeof ( char ) * gap );
                aseq_p += gap; 
                seq_p += gap;
              }
              else if ( gap < 0 ) {
                fprintf ( stderr, "The tandem set isn't ordered.\n" );
                exit ( EXIT_FAILURE );
              }
              unsigned char in = tandem->set[t].rep;
              unsigned char out = source_generate (
                  &in,
                  1,
                  tandem->set[t].pat,
                  model->amplification );
              ( in < out ) ? ++amp : ++deamp;
              int rep = out;
              memcpy (
                  &amplified_seq[aseq_p],
                  &seq[i]->seq.s[seq_p],
                  sizeof ( char ) * rep * tandem->set[t].pat );
              seq_p += ( tandem->set[t].rep * tandem->set[t].pat );
              aseq_p += ( rep * tandem->set[t].pat );
            }
            amplified_seq[aseq_p] = '\0';
            fprintf ( stderr, "original sequence:\t%d\n", seq_p );
            fprintf ( stderr, "amplified sequence:\t%d\n", aseq_p );
            fprintf ( stderr, "amplification:\t\t%d\t%.3f\n", amp, (amp*100.0/(amp+deamp)) );
            fprintf ( stderr, "deamplification:\t%d\t%.3f\n", deamp, (deamp*100.0/(amp+deamp)) );

            int j = 0;
            while ( j < seq[i]->seq.l && false ){
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
    free ( amplified_seq );
    free ( seq );
    tandem_set_destroy ( tandem );
    model_destroy ( model );
    exit ( EXIT_SUCCESS );
}
