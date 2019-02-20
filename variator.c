#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "fileManager.h"
#include "variator.h"
#include "parse_frequency.h"

#define ALL_N 2

allele_t * allele_init(long int size, allele_t * allele){
    allele->buffer_size = floor( size * 1.5 );
    allele->sequence = realloc( allele->sequence, (sizeof(char)) * allele->buffer_size );
    allele->pos = 0;
    allele->off = 0;
    return allele;
}

int main(int argc, char **argv){
    // FASTA
    struct filemanager *fm;
    struct sequence_t *seq;
    char * label;
    long int ref_pos = 0;
    // VCF
    bcf_srs_t * sr;
    bcf_hdr_t * hdr;
    bcf1_t * line; 
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

    // Parse arguments
    if ( argc != 4 ){
        printf("usage: variator vcf_filename fasta_filename output_filename\n");
        exit (EXIT_FAILURE);
    }

    // VCF file
    sr = bcf_sr_init ();
    bcf_sr_set_opt ( sr, BCF_SR_REQUIRE_IDX );
    if ( bcf_sr_add_reader ( sr, argv[1] ) != 1 ){
        exit ( EXIT_FAILURE );
    }

    // File is not indexed
    if ( sr->regions == NULL ){
        int res = bcf_index_build ( argv[1], 4 );
        if ( res == 0 ){
            // File is now indexed
            bcf_sr_remove_reader ( sr, 0 );
            if ( bcf_sr_add_reader ( sr, argv[1] ) != 1 ){
                exit ( EXIT_FAILURE );
            }
        }
        else{
            perror ( "File is not indexable.\n" );
            exit ( EXIT_FAILURE );
        }
    }
    // VCF header
    hdr = bcf_sr_get_header ( sr, 0 );
    // VCF line
    line = bcf_init();

    // FASTA file
    fm = filemanager_init (argv[2]);
    if ( fm == NULL ){
        exit ( EXIT_FAILURE );
    }
    
    // OUTPUT
    str = malloc( sizeof(char) * (strlen(argv[3]) + 20 ));
    for ( int i = 0; i < ALL_N; i++ ){
        sprintf(str, "%s_%d.fa", argv[3], i);        
        output[i] = fopen( str, "w+" );
    }
    
    // FASTA sequence
    seq = filemanager_next_seq (fm, NULL);
    for ( int i = 0; i < ALL_N; i++ ){
        allele[i] = malloc ( sizeof(struct allele_t) );
        allele[i]->sequence = NULL;
    }
    while ( seq != NULL ){
        // Allocazione di una struttura per ogni allele
        for ( int i = 0; i < ALL_N; i++ ){
            allele[i] = allele_init( seq->sequence_size, allele[i]);
        }
        label = strtok ( seq->label, " " );
        // Seek on VCF file
        if ( bcf_sr_seek ( sr, label, 0 ) == 0 ){
            // Fino alla fine del VCF
            while ( bcf_sr_next_line ( sr ) ){
                line = bcf_sr_get_line ( sr, 0 );
                // Lettura della linea
                if (bcf_unpack( line, BCF_UN_STR) != 0){
                    perror("Unpack error");
                    exit(EXIT_FAILURE);
                }
                // Distanza tra la variazione ed il carattere ancora da leggere
                distance = line->pos - ref_pos;
                // Per ogni allele
                for ( int i = 0; i < ALL_N; i++ ){
                    if ( distance >= 0 ){
                        // Copia in blocco [ref_pos, var_pos)
                        // dal reference all'allele
                        memcpy( 
                                &allele[i]->sequence[allele[i]->pos], 
                                &seq->sequence[ref_pos] , 
                                distance
                                );
                    }
                    // La posizione sull'allele varia al netto del segno della distanza
                    allele[i]->pos += distance;    
                    // Test
                    // La posizione del puntatore sull'allele è coerente?
                    assert ( line->pos == allele[i]->pos + allele[i]->off );
                    // Il reference è coerente con quello descritto dalla variazione?
                    ref_check = realloc (ref_check, sizeof(char) * ( strlen( line->d.allele[0] ) + 1 ) );
                    sprintf ( ref_check, "%.*s", (int) strlen ( line->d.allele[0] ), &seq->sequence[ line->pos ] ); 
                    assert ( strcmp ( ref_check, line->d.allele[0] ) == 0 );
                    // Definito dallo standard
                    af_ret = bcf_get_info_float( hdr, line, "AF", af, &af_size );
                    // Usato da dbSNP
                    freq_ret = bcf_get_info_string( hdr, line, "FREQ", &freq, &freq_size );
                    // Controllo dei risultati
                    if ( af_ret >= 0  ){
                        p = parse_af( line->n_allele, af, p );
                    }
                    else if (freq_ret >= 0){
                        p = parse_db_snp_freq( line->n_allele, freq, p );
                    }
                    else{
                        p = linear( line->n_allele, p );
                    }
                    
                    // Estrazione dell'allele
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
                    // Inserimento dei caratteri
                    memcpy(
                            &allele[i]->sequence[allele[i]->pos],
                            subseq,
                            strlen(subseq)
                          );
                    // Aggiornamento dell'offset
                    // TODO: si può recuperare questa informazione dal VCF?
                    allele[i]->off += strlen(line->d.allele[0]) - strlen(subseq);
                    allele[i]->pos += strlen(subseq);
                }
                // Aggiornamento della posizione sul reference
                ref_pos += (distance + strlen(line->d.allele[0]));
            }
            // Inserimento della parte finale del reference
            for ( int i = 0; i < ALL_N; i++ ){
                distance = seq->sequence_size - ref_pos;
                memcpy(
                        &allele[i]->sequence[allele[i]->pos],
                        &seq->sequence[ref_pos], 
                        distance
                      );
                // Aggiornamento della posizione
                allele[i]->pos += distance;
                // Carattere di fine stringa
                allele[i]->pos = 0;
            }
            // Write of the sequence on file
            for ( int i = 0; i < ALL_N; i++ ){
                fprintf ( output[i], ">%s\n", label );
                fprintf ( output[i], allele[i]->sequence );
            }
        }
        else{
            fprintf ( stderr, "Sequence %s not found in VCF\n", label );
        }
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
