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
#include <time.h>
#include "model.h"
#include "stats.h"
#include "source.h"
#include "tandem.h"

// Init kseq structure
KSEQ_INIT ( gzFile, gzread );

void usage ( char * name ) {
    fprintf ( stderr, "Usage: %s coverage error_model fastq [fastq ...]\n", name );
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
    // Amplification
    char * amplified_seq = NULL;
    tandem_set_t * tandem = NULL;
    // Coverage
    int coverage;
    int sequenced;
    // Generation
    bool single_only = true;
    int insert_size = 0;
    int orientation = 0;
    int length;
    int reverse;
    int pos;
    stats_t * curr_end;

    // Init pseudorandom generator
    srand ( time ( NULL ) );

    // Non optional arguments
    if ( argc - optind < 3 ) {
        usage ( argv[0] );
        exit ( EXIT_FAILURE );
    }

    // Coverage
    coverage = atoi ( argv[optind++] );

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

    // Check if there are pair reads
    single_only = ( model->pair->alignment->n == 0 );

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

    // Simulated read generation
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
              unsigned char in = tandem->set[t].rep;
              if ( in >= model->max_repetition ) {
                continue;
              }
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
            // Copy of the remaining sequence
            memcpy (
                &amplified_seq[aseq_p],
                &seq[i]->seq.s[seq_p],
                sizeof ( char ) * ( seq[i]->seq.l - seq_p ) );
            aseq_p += ( seq[i]->seq.l - seq_p );
            seq_p = seq[i]->seq.l;
            amplified_seq[aseq_p] = '\0';
            fprintf ( stderr, "%s\n", seq[i]->name.s );
            fprintf ( stderr, "\t(amplified):\t%d\t%d\t%.3f\n", aseq_p, seq_p, (aseq_p*100.0/seq_p));

            // Initial conditions
            sequenced = 0;
            pos = 0;
            curr_end = model->single;

            // Reach the coverage
            while ( coverage > ( sequenced / aseq_p ) ) {
              if ( curr_end == model->single ) {
                // Two bits: ++,+-,-+,--
                orientation = source_generate ( NULL, 0, 0, model->orientation );
                // Not sequenced nucleotides between pairs
                insert_size = source_generate ( NULL, 0, 0, model->insert_size );
                int lo_bound = insert_size * ( model->max_insert_size / model->size_granularity );
                int up_bound = ( insert_size + 1 ) * ( model->max_insert_size / model->size_granularity );
                insert_size = rand () % ( up_bound - lo_bound + 1 );
                insert_size += lo_bound;
              }
              else {
                // There is no insert size after mate pair
                insert_size = 0;
              }

              // Generate new read
              generated = stats_generate_read ( &amplified_seq[pos], generated, curr_end );
              length = strlen ( generated->read );
              reverse = ( curr_end == model->single ) ? ( orientation & 2 ) : ( orientation & 1 );

              // The read reached the limit of the sequence
              if ( !generated->cut ) {
                // Reverse the order of the nucleotides
                if ( reverse ) {
                  char swap;
                  for ( int j = 0; j < (length/2); j ++ ) {
                    swap = generated->read[j];
                    generated->read[j] = generated->read[length-j-1];
                    generated->read[length-j-1] = swap;
                    swap = generated->quality[j];
                    generated->quality[j] = generated->quality[length-j-1];
                    generated->quality[length-j-1] = swap;
                  }
                }

                // Adjust quality score for visualization
                for ( int j = 0; j < length; j ++ ) {
                  generated->quality[j] += 33;
                }

                // Print result to file
                printf ( "@%s %d %c\n", seq[i]->name.s, pos, ( reverse ) ? '-' : '+' );
                printf ( "%s\n", generated->read );
                printf ( "+\n" );
                printf ( "%s\n\n", generated->quality );

                // Update sequenced bases
                sequenced += length;

                // Change read-end
                curr_end = ( curr_end == model->single ) ? model->pair : model->single;
                curr_end = single_only ? model->single : curr_end;
              }

              // Update start position
              pos += ( insert_size + length );
              if ( pos >= aseq_p ) {
                // Start from the beginning
                pos = 0;
                curr_end = model->single;
              }
            }
            fprintf ( stderr, "\t(sequenced):\t%d\t%d\t%.3f\n", sequenced, aseq_p, (sequenced/aseq_p)*100.0);
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
