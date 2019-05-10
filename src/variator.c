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
#include "fileManager.h"
#include "allele.h"
#include "parse_frequency.h"
#include "wrapper.h"

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
    struct filemanager * fm;
    struct sequence_t * seq;
    // Wrapper
    wrapper_t * w;
    // Alleles
    allele_t ** allele;
    int gap;
    // Output
    FILE ** output;
    FILE ** alignment;
    char * str;
    // Statistics
    bool stats = false;
    unsigned long int done = 0;
    unsigned long int udv_collision = 0;
    unsigned long int vcf_collision = 0;

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
    output = malloc ( sizeof ( FILE * ) * ploidy );
    alignment = malloc ( sizeof ( FILE * ) * ploidy );


    // Initialize wrapper
    w = wr_init ( vcf_fn, udv_fn, ploidy );

    // FASTA file
    fm = filemanager_init ( fasta_fn );
    if ( fm == NULL ) {
        exit ( EXIT_FAILURE );
    }

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
        sprintf ( str, "%s_%d.fa", out_fn, i );
        output[i] = fopen ( str, "w+" );
        sprintf ( str, "%s_%d.fa.alg", out_fn, i );
        alignment[i] = fopen ( str, "w+" );
    }

    // Load of the first sequence
    seq = filemanager_next_seq ( fm, NULL );
    // Initialize alleles
    for ( int i = 0; i < ploidy; i++ ) {
        allele[i] = allele_init ( 0, NULL );
    }

    // While there are sequences to read in the FASTA file
    while ( seq != NULL ) {
        // Resize allele
        for ( int i = 0; i < ploidy; i++ ) {
            allele[i] = allele_init ( seq->sequence_size, allele[i] );
        }
        // Label separated by white space
        if ( stats )
            printf ( "%s\n", seq->label );
        // Seek to the desired region
        if ( wr_seek ( w, seq->label ) ) {
            // Up to the end of the region
            while ( wr_region ( w ) ) {
                if ( wr_update_wrapper ( w ) ) {
                    // Per allele
                    for ( int i = 0; i < ploidy; i++ ) {
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
                                &seq->sequence[allele[i]->ref],
                                gap
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
                &seq->sequence[allele[i]->ref],
                &seq->sequence[allele[i]->ref],
                allele[i]
            );
            // End of the sequence
            allele[i]->sequence[allele[i]->pos] = '\0';
            allele[i]->alignment[allele[i]->alg] = '\0';
        }
        // Write of the sequence on file
        for ( int i = 0; i < ploidy; i++ ) {
            fprintf ( output[i], ">%s\n", seq->label );
            fprintf ( output[i], "%s\n", allele[i]->sequence );
            fprintf ( alignment[i], ">%s\n", seq->label );
            fprintf ( alignment[i], "%s\n", allele[i]->alignment );
        }
        // Next sequence
        seq = filemanager_next_seq ( fm, seq );
    }
    if ( stats ) {
        unsigned long int sum = done + vcf_collision + udv_collision;
        printf ( "DONE:\t%lu\t%.2f\n", done, done * 100.0 / sum );
        printf ( "UDVc:\t%lu\t%.2f\n", udv_collision, udv_collision * 100.0 / sum );
        printf ( "VCFc:\t%lu\t%.2f\n", vcf_collision, vcf_collision * 100.0 / sum );
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i++ ) {
        allele_destroy ( allele[i] );
        fclose ( output[i] );
        fclose ( alignment[i] );
    }
    free ( allele );
    free ( output );
    free ( alignment );
    free ( str );
    wr_destroy ( w );
    filemanager_destroy ( fm );
    exit ( EXIT_SUCCESS );
}
