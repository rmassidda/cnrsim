/*
 * CNRSIM
 * variator.c
 * Generatos two alleles from a FASTA reference
 * and a VCF file containing known variations.
 *
 * @author Riccardo Massidda
 */

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <libgen.h>
#include <zlib.h>
#include <htslib/kseq.h>
#include <time.h>
#include "allele.h"
#include "parse_frequency.h"
#include "wrapper.h"

// Init kseq structure
KSEQ_INIT ( gzFile, gzread );

void usage ( char * name ) {
    fprintf ( stderr, "Usage: %s [-n number of alleles] [-u udv_file] [-o output_name] fasta_file vcf_file\n", name );
}

int main ( int argc, char ** argv ) {
    // Parsing
    int opt;
    int ploidy = 2;
    // Filenames
    char * fasta_fn = NULL;
    char * udv_fn = NULL;
    char * vcf_fn = NULL;
    char * out_fn = NULL;
    // FASTA
    gzFile fp;
    kseq_t * seq;
    // Wrapper
    wrapper_t * w;
    // Alleles
    allele_t ** allele;
    int gap;
    // Output
    FILE ** output;
    char * str;
    // Statistics
    bool stats = false;
    unsigned long int done = 0;
    unsigned long int igno = 0;
    unsigned long int udv_collision = 0;
    unsigned long int vcf_collision = 0;

    // Init pseudorandom generator
    srand ( time ( NULL ) );

    while ( ( opt = getopt ( argc, argv, "sn:u:o:" ) ) != -1 ) {
        switch ( opt ) {
        case 's':
            stats = true;
            break;
        case 'n':
            ploidy = atoi ( optarg );
            break;
        case 'u':
            udv_fn = optarg;
            break;
        case 'o':
            out_fn = optarg;
            break;
        case '?':
            if ( optopt == 'u' || optopt == 'o' )
                fprintf ( stderr, "Option -%c requires an argument.\n", optopt );
            else if ( isprint ( optopt ) )
                fprintf ( stderr, "Unknown option `-%c'.\n", optopt );
            else
                fprintf ( stderr, "Unknown option character `\\x%x'.\n", optopt );
            exit ( EXIT_FAILURE );
        default:
            usage ( argv[0] );
            exit ( EXIT_FAILURE );
        }
    }
    // Non optional arguments
    if ( argc - optind < 2 ) {
        usage ( argv[0] );
        exit ( EXIT_FAILURE );
    }
    fasta_fn = argv[optind++];
    vcf_fn = argv[optind];

    // Allocate
    allele = malloc ( sizeof ( allele_t * ) * ploidy );
    output = malloc ( sizeof ( FILE * ) * 2 * ploidy );


    // Initialize wrapper
    w = wr_init ( vcf_fn, udv_fn, ploidy );

    // FASTA file
    fp = gzopen ( fasta_fn, "r" );
    if ( fp == NULL ) {
        fprintf ( stderr, "File %s not found.\n", fasta_fn );
        exit ( EXIT_FAILURE );
    }
    seq = kseq_init ( fp );

    /*
     * Output files, one per allele
     * [filename_0.fa, filename_N.fa)
     */
    if ( out_fn == NULL ) {
        out_fn = basename ( fasta_fn );
        out_fn = strtok ( out_fn, "." );
    }
    str = malloc ( sizeof ( char ) * ( strlen ( out_fn ) + 20 ) );
    for ( int i = 0; i < ploidy; i++ ) {
        sprintf ( str, "%s_%d.fq", out_fn, i );
        output[i] = fopen ( str, "w+" );
    }
    for ( int i = 0; i < ploidy; i++ ) {
        sprintf ( str, "%s_%d.fa", out_fn, i );
        output[i+ploidy] = fopen ( str, "w+" );
    }

    // Initialize alleles
    for ( int i = 0; i < ploidy; i++ ) {
        allele[i] = allele_init ( 0, NULL );
    }

    // While there are sequences to read in the FASTA file
    while ( kseq_read ( seq ) >= 0 ) {
        // Resize allele
        for ( int i = 0; i < ploidy; i++ ) {
            allele[i] = allele_init ( seq->seq.l, allele[i] );
        }
        // Label separated by white space
        if ( stats )
            printf ( "%s\n", seq->name.s );
        // Seek to the desired region
        if ( wr_seek ( w, seq->name.s ) ) {
            // Up to the end of the region
            while ( wr_region ( w ) ) {
                if ( wr_update_wrapper ( w ) ) {
                    // Per allele
                    for ( int i = 0; i < ploidy; i++ ) {
                        // Don't "apply" reference
                        if ( w->alt_index[i] < 0 ) {
                            igno ++;
                            continue;
                        }

                        // Gap between variations
                        gap = w->pos - allele[i]->ref;

                        // Avoid collision
                        if ( gap < 0 ){
                            //  Keep track of the collision number
                            if ( stats ){
                                vcf_collision ++;
                            }
                            continue;
                        }

                        // Distance between the reference and the variation pointers
                        if ( gap > 0 ) {
                            /*
                             * The variation starts far from the current
                             * reference position, what is in between can
                             * be copied without any mutation.
                             */
                            memcpy (
                                &allele[i]->sequence[allele[i]->pos],
                                &seq->seq.s[allele[i]->ref],
                                sizeof ( char ) * gap
                            );
                            memset ( 
                                &allele[i]->alignment[allele[i]->alg],
                                '=',
                                sizeof ( char ) * gap
                            );
                            // Update position
                            allele[i]->pos += gap;
                            allele[i]->alg += gap;
                            allele[i]->ref += gap;
                        }

                        // Alternative
                        allele_variation (
                                        w->ref,
                                        w->alt[w->alt_index[i]],
                                        allele[i] );

                        done ++;
                    }
                } else if ( stats ) {
                    udv_collision++;
                }
            }
        }
        // Copy of the remaining part of the sequence
        for ( int i = 0; i < ploidy; i++ ) {
            allele_variation (
                &seq->seq.s[allele[i]->ref],
                &seq->seq.s[allele[i]->ref],
                allele[i]
            );
            // End of the sequence
            allele[i]->sequence[allele[i]->pos] = '\0';
            allele[i]->alignment[allele[i]->alg] = '\0';
        }
        // Write of the sequence on file
        for ( int i = 0; i < ploidy; i++ ) {
            fprintf ( output[i], ">%s\n", seq->name.s );
            fprintf ( output[i], "%s\n", allele[i]->sequence );
            fprintf ( output[i], "+\n" );
            fprintf ( output[i], "%s\n", allele[i]->alignment );
        }
        for ( int i = 0; i < ploidy; i++ ) {
            fprintf ( output[i+ploidy], ">%s\n", seq->name.s );
            fprintf ( output[i+ploidy], "%s\n", allele[i]->sequence );
        }
    }
    if ( stats ) {
        unsigned long int sum = done + igno + vcf_collision + udv_collision;
        printf ( "DONE:\t%lu\t%.2f\n", done, done * 100.0 / sum );
        printf ( "REF:\t%lu\t%.2f\n", igno, igno * 100.0 / sum );
        printf ( "UDVc:\t%lu\t%.2f\n", udv_collision, udv_collision * 100.0 / ( sum - igno ) );
        printf ( "VCFc:\t%lu\t%.2f\n", vcf_collision, vcf_collision * 100.0 / ( sum - igno ) );
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i++ ) {
        allele_destroy ( allele[i] );
    }
    for ( int i = 0; i < ( 2 * ploidy ); i++ ) {
        fclose ( output[i] );
    }
    free ( allele );
    free ( output );
    free ( str );
    kseq_destroy ( seq );
    gzclose ( fp );
    wr_destroy ( w );
    exit ( EXIT_SUCCESS );
}
