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

#define ALL_N 1


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
    // VCF
    bcf_srs_t * sr;
    bcf_hdr_t * hdr;
    bcf1_t * line; 
    int current_region = 0;
    bool first_line = false;
    // Alleles
    allele_t * allele[ ALL_N ];
    long int distance;
    char * subseq;
    // Output
    FILE * output[ ALL_N ];
    char * str;
    // Random number generator
    srand(time(NULL));
    double outcome;
    double threshold;
    //  Allelic frequency parser
    float *af = NULL;
    int af_size = 0;
    int af_ret;
    char *freq = NULL;
    int freq_size = 0;
    int freq_ret;
    double * p = NULL;
    // Test
    char * ref_check = NULL;
    char * all_check = NULL;
    int done = 0;
    int ignored = 0;

    // Parse arguments
    if ( argc != 4 ){
        printf("usage: variator vcf_filename fasta_filename output_filename\n");
        exit (EXIT_FAILURE);
    }


    // Initialize VCF readers
    sr = bcf_sr_init ();
    // VCF has to be indexed
    bcf_sr_set_opt ( sr, BCF_SR_REQUIRE_IDX );
    // Link the reader to the VCF file passed by argument
    if ( bcf_sr_add_reader ( sr, argv[1] ) != 1 ){
        exit ( EXIT_FAILURE );
    }

    /*
     * If the file is indexed the name
     * and the position of the regions
     * are loaded into the readers
     * structure.
     */
    if ( sr->regions == NULL ){
        // File not indexed, creation of the index
        // 14 is a suggested value by HTSLIB documentation
        int res = bcf_index_build ( argv[1], 14 );
        if ( res == 0 ){
            // File is now indexed
            // Relink of the file with the reader
            bcf_sr_remove_reader ( sr, 0 );
            if ( bcf_sr_add_reader ( sr, argv[1] ) != 1 ){
                exit ( EXIT_FAILURE );
            }
        }
        else{
            perror ( "File is not indexable.\nTry first compressing it with gzip.\n" );
            exit ( EXIT_FAILURE );
        }
    }

    /*
     * Load of the header, needed to parse
     * INFO values in the VCF.
     */
    hdr = bcf_sr_get_header ( sr, 0 );
    // Initialize structure that contains a record
    line = bcf_init();

    // FASTA file
    fm = filemanager_init (argv[2]);
    if ( fm == NULL ){
        exit ( EXIT_FAILURE );
    }
    
    /*
     * Output files, one per allele
     * [filename_0.fa, filename_N.fa)
     */
    str = malloc( sizeof(char) * (strlen(argv[3]) + 20 ));
    for ( int i = 0; i < ALL_N; i++ ){
        sprintf(str, "%s_%d.fa", argv[3], i);        
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
        // Seek on VCF file
        if ( bcf_sr_seek ( sr, label, 0 ) == 0 ){
            // Label found
            first_line = true; 
            // Up to the end of VCF file
            while ( bcf_sr_next_line ( sr ) ){
                line = bcf_sr_get_line ( sr, 0 );
                // Set the current region ID
                if ( first_line ){
                    current_region = line->rid;
                    first_line = false;
                }
                // If the region changed we need to stop
                if ( current_region != line->rid ){
                    printf ( "End of %s region\n", label );
                    break;
                }
                // Unpack up to INFO of the record
                if (bcf_unpack( line, BCF_UN_STR) != 0){
                    perror("Unpack error");
                    exit(EXIT_FAILURE);
                }
                // Per allele
                for ( int i = 0; i < ALL_N; i++ ){
                    // Distance between the reference and the variation pointers
                    distance = line->pos - ref_pos[i];
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
                    assert ( line->pos == allele[i]->pos + allele[i]->off );
                    /*
                     * Reference sequence and reference in the VCF
                     * MUST coincide.
                     */
                    ref_check = &seq->sequence[line->pos];
                    // assert ( strncasecmp ( ref_check, line->d.allele[0], strlen ( line->d.allele[0] ) ) == 0 );
                    // Questo assert viene momentaneamente tolto, in quanto strncasecmp non supporta il carattere
                    // N o il carattere n
                    /*
                    * If we want to applicate a certain variation,
                    * reference in the allele and VCF reference
                    * have to coincide.
                    */
                    all_check = &(allele[i]->sequence[allele[i]->pos]);

                    if ( distance <= 0 ){
                        // fprintf ( stderr, "%s\t%s\t", all_check, line->d.allele[0] );
                        // The variation describes something that is already written
                        if ( strncasecmp ( all_check, line->d.allele[0], strlen (all_check) ) != 0 ){
                            // Questo controllo può essere falsato dal fatto che strncasecmp non riconosca
                            // il carattere N o n come "matcha tutto"
                            // The reference and the allele doesn't match
                            // due to previous variations
                            allele[i]->pos -= distance;
                            ignored ++;
                            // fprintf ( stderr, "ERR\n" );
                            continue;
                        }
                        else{
                            // fprintf ( stderr, "OK\n" );
                        }
                    }
                    done ++;
                    
                    // Allelic frequency as defined by VCF
                    af_ret = bcf_get_info_float( hdr, line, "AF", af, &af_size );
                    // Allelic frequency as defined by dbSNP
                    freq_ret = bcf_get_info_string( hdr, line, "FREQ", &freq, &freq_size );
                    // Parse results
                    if ( af_ret >= 0  ){
                        p = parse_af( line->n_allele, af, p );
                    }
                    else if (freq_ret >= 0){
                        p = parse_db_snp_freq( line->n_allele, freq, p );
                    }
                    else{
                        p = linear( line->n_allele, p );
                    }
                    
                    // Random decision about the alternatives 
                    outcome = (double)rand() / RAND_MAX;
                    threshold = 0;
                    for (int i = 0; i < line->n_allele; i++){
                        if ( threshold <= outcome && outcome < threshold + p[i] ){
                            subseq = line->d.allele[i];
                            break;
                        }
                        else{
                            threshold += p[i];
                        }
                    }
                    // Application of the variation
                    memcpy(
                            &allele[i]->sequence[allele[i]->pos],
                            subseq,
                            strlen(subseq)
                          );
                    // Update of the offset and the position of the allele
                    allele[i]->off += strlen(line->d.allele[0]) - strlen(subseq);
                    allele[i]->pos += strlen(subseq);
                    // Reference position update
                    ref_pos[i] += (distance + strlen(line->d.allele[0]));
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
        printf ( "DONE:\t%d\t%.2f\n", done, done*100.0/(done+ignored ) );
        printf ( "IGNO:\t%d\t%.2f\n", ignored, ignored*100.0/(done+ignored ) );
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
    free ( p );
    bcf_destroy( line );
    bcf_hdr_destroy( hdr );
    filemanager_destroy( fm );
    exit(EXIT_SUCCESS);
}
