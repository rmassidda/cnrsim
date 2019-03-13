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

void usage ( char * name){
    fprintf(stderr, "Usage: %s [-n number of alleles] [-u udv_file] [-o output_name] fasta_file vcf_file\n", name );
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
    long int * ref_pos;
    // Wrapper
    wrapper_t * w;
    // Alleles
    allele_t ** allele;
    long int distance;
    char * subseq;
    int ref_len;
    int alt_len;
    int offset;
    // Output
    FILE ** output;
    FILE ** alignment;
    char * str;
    // Statistics
    char * all_check = NULL;
    unsigned long int done = 0;
    unsigned long int ignored = 0;
    unsigned long int udv_collision = 0;
    unsigned long int less_than_zero = 0;

    while ((opt = getopt(argc, argv, "n:u:o:")) != -1) {
        switch (opt) {
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
                if (optopt == 'u' || optopt == 'o')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                exit ( EXIT_FAILURE );
            default:
                usage ( argv[0] );
                exit(EXIT_FAILURE);
        }
    }
    // Non optional arguments
    if ( argc - optind < 2 ){
        usage ( argv[0] );
        exit(EXIT_FAILURE);
    }
    fasta_fn = argv[optind++];
    vcf_fn = argv[optind];

    // Allocate
    ref_pos = malloc ( sizeof ( long int ) * ploidy );
    allele = malloc ( sizeof ( allele_t* ) * ploidy );
    output = malloc ( sizeof ( FILE* ) * ploidy );
    alignment = malloc ( sizeof ( FILE* ) * ploidy );


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
    if ( out_fn == NULL ){
        out_fn = basename ( fasta_fn );
        out_fn = strtok ( out_fn, "." );
    }
    str = malloc ( sizeof ( char ) * ( strlen ( out_fn ) + 20 ) );
    for ( int i = 0; i < ploidy; i++ ) {
        sprintf ( str, "%s_%d.fa", out_fn, i );
        output[i] = fopen ( str, "w+" );
        sprintf ( str, "%s_%d.fa.alg",out_fn, i );
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
            ref_pos[i] = 0;
        }
        // Label separated by white space
        fprintf ( stdout, "%s\n", seq->label );
        // Seek to the desired region
        if ( wr_seek ( w, seq->label ) ) {
            // Up to the end of the region
            while ( wr_region ( w ) ) {
                if ( wr_update_wrapper ( w ) ) {
                    // Per allele
                    for ( int i = 0; i < ploidy; i++ ) {
                        // Distance between the reference and the variation pointers
                        distance = w->pos - ref_pos[i];
                        if ( distance > 0 ) {
                            /*
                             * The variation starts far from the current
                             * reference position, what is in between can
                             * be copied without any mutation.
                             */
                            memcpy (
                                &allele[i]->sequence[allele[i]->pos],
                                &seq->sequence[ref_pos[i]],
                                distance
                            );
                            memset ( 
                                &allele[i]->alignment[allele[i]->pos],
                                '=',
                                distance
                            );
                        }
                        /*
                         * The allele pointer points to the start of
                         * the variation, according to the offset.
                         */
                        allele[i]->pos += distance;
                        assert ( w->pos == allele[i]->pos + allele[i]->off );
                        // The position hasn't been written yet
                        if ( allele[i]->pos < 0 ){
                            allele[i]->pos -= distance;
                            less_than_zero ++;
                            continue;
                        }
                        /*
                        * If we want to applicate a certain variation,
                        * reference in the allele and VCF reference
                        * have to coincide.
                        */
                        all_check = & ( allele[i]->sequence[allele[i]->pos] );

                        if ( distance <= 0 ) {
                            // The variation describes something that is already written
                            if ( strncasecmp ( all_check, w->ref, strlen ( all_check ) ) != 0 ) {
                                // The reference and the allele doesn't match
                                // due to previous variations
                                allele[i]->pos -= distance;
                                ignored ++;
                                continue;
                            }
                        }
                        done ++;

                        subseq = w->alt[w->alt_index[i]];
                        // Application of the variation
                        ref_len = strlen ( w->ref );
                        alt_len = strlen ( subseq );
                        offset = ref_len - alt_len;
                        int min = ( ref_len < alt_len ) ? ref_len : alt_len;
                        char c = ( ref_len < alt_len ) ? 'd' : 'i';
                        memcpy (
                            &allele[i]->sequence[allele[i]->pos],
                            subseq,
                            alt_len
                        );
                        memset ( 
                            &allele[i]->alignment[allele[i]->pos],
                            '=',
                            min
                        );
                        memset ( 
                            &allele[i]->alignment[allele[i]->pos + min],
                            c,
                            abs ( offset )
                        );
                        // Update of the offset and the position of the allele
                        allele[i]->off += offset;
                        allele[i]->pos += alt_len;
                        // Reference position update
                        ref_pos[i] += ( distance + ref_len );
                    }
                } else {
                    udv_collision++;
                }
            }
        }
        // Copy of the remaining part of the sequence
        for ( int i = 0; i < ploidy; i++ ) {
            distance = seq->sequence_size - ref_pos[i];
            memcpy (
                &allele[i]->sequence[allele[i]->pos],
                &seq->sequence[ref_pos[i]],
                distance
            );
            memset ( 
                &allele[i]->alignment[allele[i]->pos],
                '=',
                distance
            );
            // Update position
            allele[i]->pos += distance;
            // End of the sequence
            allele[i]->pos = 0;
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
    unsigned long int sum = done + ignored + udv_collision + less_than_zero;
    printf ( "DONE:\t%lu\t%.2f\n", done, done * 100.0 / sum );
    printf ( "IGNO:\t%lu\t%.2f\n", ignored, ignored * 100.0 / sum );
    printf ( "UDVC:\t%lu\t%.2f\n", udv_collision, udv_collision * 100.0 / sum );
    printf ( "LESS:\t%lu\t%.2f\n", less_than_zero, less_than_zero * 100.0 / sum );

    // Cleanup
    for ( int i = 0; i < ploidy; i++ ) {
        allele_destroy ( allele[i] );
        fclose ( output[i] );
        fclose ( alignment[i] );
    }
    free ( ref_pos );
    free ( allele );
    free ( output );
    free ( alignment );
    free ( str );
    wr_destroy ( w );
    filemanager_destroy ( fm );
    exit ( EXIT_SUCCESS );
}
