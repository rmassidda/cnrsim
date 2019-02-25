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
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "fileManager.h"
#include "variator.h"
#include "parse_frequency.h"
#include "wrapper.h"

allele_t * allele_init(long int size, allele_t * allele){
    // First initialization
    if ( allele == NULL ){
        allele = malloc ( sizeof(allele_t) );
        allele->sequence = NULL;
    }
    // Update of internal values
    allele->buffer_size = floor( size * 1.5 );
    allele->sequence = realloc( allele->sequence, (sizeof(char)) * allele->buffer_size );
    // It's necessary to clean the memory
    memset ( allele->sequence, 0, sizeof(char)*allele->buffer_size );
    allele->pos = 0;
    allele->off = 0;
    return allele;
}

int main(int argc, char **argv){
    // FASTA
    struct filemanager *fm;
    struct sequence_t *seq;
    char * label;
    long int ref_pos[ALL_N];
    // Wrapper
    wrapper_t * w;
    // Alleles
    allele_t * allele[ ALL_N ];
    long int distance;
    char * subseq;
    // Output
    FILE * output[ ALL_N ];
    char * str;
    // Statistics
    char * all_check = NULL;
    int done = 0;
    int ignored = 0;
    int udv_collision = 0;

    // Parse arguments
    if ( argc != 5 ){
        printf("usage: variator vcf_filename udv_filename fasta_filename output_filename\n");
        exit (EXIT_FAILURE);
    }


    // Initialize wrapper
    w = wr_init ( argv[1], argv[2] );

    // FASTA file
    fm = filemanager_init (argv[3]);
    if ( fm == NULL ){
        exit ( EXIT_FAILURE );
    }
    
    /*
     * Output files, one per allele
     * [filename_0.fa, filename_N.fa)
     */
    str = malloc( sizeof(char) * (strlen(argv[4]) + 20 ));
    for ( int i = 0; i < ALL_N; i++ ){
        sprintf(str, "%s_%d.fa", argv[4], i);        
        output[i] = fopen( str, "w+" );
    }
    
    // Load of the first sequence
    seq = filemanager_next_seq (fm, NULL);
    // Initialize alleles
    for ( int i = 0; i < ALL_N; i++ ){
        allele[i] = allele_init (0, NULL);
    }

    // While there are sequences to read in the FASTA file
    while ( seq != NULL ){
        // Resize allele
        for ( int i = 0; i < ALL_N; i++ ){
            allele[i] = allele_init( seq->sequence_size, allele[i]);
            ref_pos[i] = 0;
        }
        // Label separated by white space
        label = strtok ( seq->label, " " );
        fprintf ( stdout, "%s\n", seq->label );
        // Seek to the desired region
        if ( wr_seek ( w, label ) ){
            // Up to the end of the region
            while ( wr_region ( w ) ){
                if ( wr_update_wrapper ( w ) ){
                    // Per allele
                    for ( int i = 0; i < ALL_N; i++ ){
                        // Distance between the reference and the variation pointers
                        distance = w->pos - ref_pos[i];
                        if ( distance > 0 ){
                            /*
                             * The variation starts far from the current
                             * reference position, what is in between can
                             * be copied without any mutation.
                             */
                            memcpy( 
                                    &allele[i]->sequence[allele[i]->pos], 
                                    &seq->sequence[ref_pos[i]] , 
                                    distance
                                    );
                        }
                        /*
                         * The allele pointer points to the start of
                         * the variation, according to the offset.
                         */
                        allele[i]->pos += distance;    
                        assert ( w->pos == allele[i]->pos + allele[i]->off );
                        /*
                        * If we want to applicate a certain variation,
                        * reference in the allele and VCF reference
                        * have to coincide.
                        */
                        all_check = &(allele[i]->sequence[allele[i]->pos]);

                        if ( distance <= 0 ){
                            // The variation describes something that is already written
                            if ( strncasecmp ( all_check, w->ref, strlen (all_check) ) != 0 ){
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
                        memcpy(
                                &allele[i]->sequence[allele[i]->pos],
                                subseq,
                                strlen(subseq)
                              );
                        // Update of the offset and the position of the allele
                        allele[i]->off += strlen( w->ref ) - strlen(subseq);
                        allele[i]->pos += strlen(subseq);
                        // Reference position update
                        ref_pos[i] += (distance + strlen( w->ref ));
                    }
                }
                else{
                    udv_collision++;
                }
            }
            // Copy of the remaining part of the sequence
            for ( int i = 0; i < ALL_N; i++ ){
                distance = seq->sequence_size - ref_pos[i];
                memcpy(
                        &allele[i]->sequence[allele[i]->pos],
                        &seq->sequence[ref_pos[i]], 
                        distance
                      );
                // Update position
                allele[i]->pos += distance;
                // End of the sequence
                allele[i]->pos = 0;
            }
            // Write of the sequence on file
            for ( int i = 0; i < ALL_N; i++ ){
                fprintf ( output[i], ">%s\n", label );
                fprintf ( output[i], "%s\n", allele[i]->sequence );
                printf ( "%s %d writed on file.\n", label, i );
            }
        }
        else{
            fprintf ( stderr, "Sequence %s not found in VCF\n", label );
        }
        int sum = done + ignored + udv_collision;
        printf ( "DONE:\t%d\t%.2f\n", done, done*100.0/sum );
        printf ( "IGNO:\t%d\t%.2f\n", ignored, ignored*100.0/sum );
        printf ( "UDVC:\t%d\t%.2f\n", udv_collision, udv_collision*100.0/sum );
        done = 0;
        ignored = 0;
        // Next sequence
        seq = filemanager_next_seq (fm, NULL);
    }

    // Cleanup
    for ( int i = 0; i < ALL_N; i++ ){
        free( allele[i]->sequence );
        free( allele[i] );
        fclose( output[i] );
    }
    // free ( p );
    // bcf_destroy( line );
    // bcf_hdr_destroy( hdr );
    filemanager_destroy( fm );
    exit(EXIT_SUCCESS);
}
