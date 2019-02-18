#include <htslib/vcf.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "fileManager.h"
#include "variator.h"
#include "parse_frequency.h"

#define ALL_N 2

allele_t * allele_init(long int size){
    allele_t * allele = malloc((sizeof(struct allele_t)));
    allele->buffer_size = floor( size * 1.5 );
    allele->sequence = malloc((sizeof(char)) * allele->buffer_size );
    allele->pos = 0;
    allele->off = 0;
    return allele;
}

int main(int argc, char **argv){
    // FASTA
    struct filemanager *fm;
    struct sequence_t *seq;
    long int ref_pos = 0;
    // VCF
    htsFile * inf;
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

    // Parse arguments
    if ( argc != 4 ){
        printf("usage: variator vcf_filename fasta_filename output_filename\n");
        exit (EXIT_FAILURE);
    }

    // VCF file
    inf = bcf_open (argv[1], "r" );
    // VCF header
    hdr = bcf_hdr_read(inf);
    // VCF line
    line = bcf_init();

    // FASTA file
    fm = filemanager_init (argv[2]);
    // FASTA sequence
    seq = filemanager_next_seq (fm, NULL);

    // OUTPUT
    str = malloc( sizeof(char) * (strlen(argv[3]) + 20 ));
    for ( int i = 0; i < ALL_N; i++ ){
        sprintf(str, "%s_%d.fa", argv[3], i);        
        output[i] = fopen( str, "w+" );
    }
    
    // Allocazione di una struttura per ogni allele
    for ( int i = 0; i < ALL_N; i++ ){
        allele[i] = allele_init( seq->sequence_size );
    }

    // Fino alla fine del VCF
    while (bcf_read(inf, hdr, line) == 0){
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
            
            assert ( line->pos == allele[i]->pos + allele[i]->off );
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

    // Scrittura su file
    for ( int i = 0; i < ALL_N; i++ ){
        fprintf( output[i], allele[i]->sequence );
        free( allele[i]->sequence );
        free( allele[i] );
        fclose( output[i] );
    }
    bcf_destroy( line );
    bcf_hdr_destroy( hdr );
    bcf_close( inf );
    filemanager_destroy( fm );
    exit(EXIT_SUCCESS);
}
