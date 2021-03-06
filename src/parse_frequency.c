/*
 * CNRSIM
 * parse_frequency.c
 * Library that generates a discrete
 * distribution of alternative alleles
 * givent the allele frequency
 *
 * @author Riccardo Massidda
 */
#include "parse_frequency.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

double * linear ( int n, double * p ) {
    p = realloc ( p, sizeof ( double ) * n );
    double val = 1.0 / n;
    for ( int i = 0; i < n; i++ ) {
        p[i] = val;
    }
    return p;
}

double * parse_db_snp_freq ( int n, char * freq, double * p ) {
    p = realloc ( p, sizeof ( double ) * n );
    char * research = strtok ( freq, "|" );
    char * substr = strtok ( research, ":" );
    int n_dots = 0;
    double sum = 0;
    double norm = 0;

    // Read of data
    for ( int i = 0; i < n; i++ ) {
        substr = strtok ( NULL, "," );
        if ( substr[0] == '.' ) {
            p[i] = -1;
            n_dots++;
        } else {
            sscanf ( substr, "%le", &p[i] );
            sum += p[i];
        }
    }

    // Substitution of points (no info)
    for ( int i = 0; i < n; i++ ) {
        if ( p[i] == -1 ) {
            p[i] = ( 1 - sum ) / n_dots;
        }
        norm += p[i];
    }

    // Normalization
    for ( int i = 0; i < n; i++ ) {
        p[i] /= norm;
    }

    return p;
}

double * parse_af ( int n, float * freq, double * p ) {
    p = realloc ( p, sizeof ( double ) * n );
    for ( int i = 0; i < n; i++ ) {
        p[i] = ( double ) freq[i];
    }
    return p;
}
