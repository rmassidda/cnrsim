/*
 * Second experiment with HTSLIB.
 * Try to read secondary informations from a VCF file.
 * This could be useful to parse genotypes, allel frequency and so on. 
 *
 * Achievement:
 * - Truly understand VCF file format
 * - Wider use of vcf.h and vcfutils.h 
 */

#include <string.h>
#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>

double * linear ( int n, double * p ){
    p = realloc(p, sizeof(double) * n );
    double val = 1.0 / n;
    for ( int i = 0; i < n; i++ ){
        p[i] = val;
    }
    return p;
}

double * parse_db_snp_freq( int n, char * freq, double * p ){
    p = realloc(p, sizeof(double) * n);
    char * research = strtok ( freq, "|" );
    char * substr = strtok ( research, ":" );
    int n_dots = 0;
    double sum = 0;
    double norm = 0;

    // Lettura dei dati noti
    for ( int i = 0; i < n; i++ ){
        substr = strtok ( NULL, "," );
        if ( substr[0] == '.' ){
            p[i] = -1;
            n_dots++;
        }
        else{
            sscanf( substr, "%le", &p[i] );
            sum += p[i];
        }
    }

    // Sostituzione dei punti
    for ( int i = 0; i < n; i++ ){
        if ( p[i] == -1 ){
            p[i] = (1 - sum) / n_dots;
        }
        norm += p[i];
    }

    // Normalizzazione dei risultati
    for ( int i = 0; i < n; i++ ){
        p[i] /= norm;
    }
    
    return p;
}

int main(int argc, char **argv){
    float *af = NULL;
    int af_size = 0;
    int af_ret;
    char *freq = NULL;
    int freq_size = 0;
    int freq_ret;
    char *substr;
    double * p = NULL;
    // Filename from argument
    if (argc != 2){
        printf("Filename required.\n");
        return -1;
    }
    // Open File
    htsFile * inf = bcf_open (argv[1], "r" );
    // Read of VCF header
    bcf_hdr_t * hdr = bcf_hdr_read(inf);
    if ( hdr == NULL ){
        return -1;
    } 
    // Struct for storing each line
    bcf1_t * line = bcf_init();
    if (line == NULL){
        return -1;
    }
    // Iteration on the line set
    while (bcf_read(inf, hdr, line) == 0){
        // Unpack 'till ALT 
        if (bcf_unpack( line, BCF_UN_INFO) != 0){
            printf("Unpack error");
            return -1;
        }
        printf("%d\t", line->pos);
        // Definito dallo standard
        af_ret = bcf_get_info_float( hdr, line, "AF", af, &af_size );
        // Usato da dbSNP
        freq_ret = bcf_get_info_string( hdr, line, "FREQ", &freq, &freq_size );
        if ( af_ret >= 0  ){
            p = realloc( p, sizeof(double) * af_size );
            memcpy( p, af, af_size );
            printf( "AF\t" );
        }
        else if (freq_ret >= 0){
            printf( "FR\t" );
            printf( "%s\t", freq );
            p = parse_db_snp_freq( line->n_allele, freq, p );
        }
        else{
            p = linear( line->n_allele, p );
            printf( "LI\t" );
        }
        // Estrazione
        for ( int i = 0; i < line->n_allele; i++ ){
            printf("%f\t", p[i]);
        }
        printf("\n");
    }

    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    return 0;
}
