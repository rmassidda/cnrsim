#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "align.h"

aligner_t * al_init ( char * reference, char * read, method_t method ){
    aligner_t * al = malloc ( sizeof ( aligner_t ) );
    // Sequences
    al->reference = reference;
    al->read = read;
    al->ref_len = strlen ( reference ) + 1;
    al->read_len = strlen ( read ) + 1;
    // Matrix
    al->nw = NULL;
    al->op = NULL;
    // Alignement result
    al->alignment = NULL;
    al->start = 0;
    // Init matrixes
    __initMatrix ( al );
    return al;
}

void __initMatrix ( aligner_t * al ) {
    int i, j;
    // Allocate matrix pointers
    al->nw = malloc ( al->read_len * sizeof ( int * ) );
    al->op = malloc ( al->read_len * sizeof ( char * ) );
    // Allocate rows
    for ( i = 0; i < al->read_len; i ++ ) {
        al->nw[i] = malloc ( al->ref_len * sizeof ( int ) );
        al->op[i] = malloc ( al->ref_len * sizeof ( char ) );
        // Initial conditions
        for ( j = 0; j < al->ref_len; j ++ ) {
            al->nw[i][j] = 0;
            al->op[i][j] = '!';
        }
        // GAP in the first column
        al->nw[i][0] = i * GAP;
        al->op[i][0] = 'T';
    }
}

void align ( aligner_t * al ) {
    int left, top;
    for ( int i = 1; i < al->read_len; i ++ ) {
        for ( int j = 1; j < al->ref_len; j ++ ) {
            left = al->nw[i][j-1] + GAP;
            top = al->nw[i-1][j] + GAP;
            // Search for max neighbors
            if ( al->reference[j-1] == al->read[i-1]){
                al->nw[i][j] = al->nw[i - 1][j - 1] + MATCH;
                al->op[i][j] = '=';
            }
            else {
                al->nw[i][j] = al->nw[i - 1][j - 1] + MISMATCH;
                al->op[i][j] = '!';
            }
            if ( top > al->nw[i][j] ){
                al->nw[i][j] = top;
                al->op[i][j] = 'T';
            }
            if ( left > al->nw[i][j] ){
                al->nw[i][j] = left;
                al->op[i][j] = 'L';
            }
        }
    }
}

char * build_alignment ( aligner_t * al ) {
    int i, j, k;
    char * s = malloc ( sizeof ( char ) * ( al->read_len + al->ref_len - 1 ) );

    // Search for the maximum value
    int max_val = al->nw[al->read_len - 1][0];
    int max_i = al->read_len - 1;
    int max_j = 0;
    for ( j = 0; j < al->ref_len; j ++ ){
        if ( al->nw[max_i][j] >= max_val ){
           max_j = j;
           max_val = al->nw[max_i][max_j];
        }
    }

    // Alignment
    i = max_i;
    j = max_j;
    k = 0;
    while ( j != 0 && i != 0 ) {
        s[k] = al->op[i][j];
        if ( s[k] == 'L' ){
            j --;
        }
        else if ( s[k] == 'T' ){
            i --;
        }
        else{
            i --;
            j --;
        }
        k ++;
    }
    while ( i > 0 ) {
        s[k] = 'T';
        i --;
        k ++;
    }
    s[k] = '\0';
    al->start = j;

    // Reverse string
    al->alignment = realloc ( al->alignment, sizeof ( char ) * ( k + 1 ) );
    for ( i = 0, k --; k >= 0; i++, k -- ) al->alignment[i] = s[k];
    al->alignment[i] = '\0';
    free ( s );

    return al->alignment;
}

char * alignment ( aligner_t * al ) {
    char * seq;
    int i = 0; // read index
    int j = 0; // reference index
    int k; // alignment index
    seq = malloc ( sizeof ( char ) * ( al->read_len + al->ref_len - 1 ) );
    for ( k = 0; k < al->start; k ++ ){
        seq[k] = ' ';
    }
    for ( k = 0; k < strlen ( al->alignment ); k ++ ) {
        if ( al->alignment[k] == '!' || al->alignment[k] == '=' ) {
            seq[k + al->start] = al->reference[j + al->start];
            i ++;
            j ++;
        }
        // Read presents an insertion
        else if ( al->alignment[k] == 'T' ) {
            seq[k + al->start] = al->read[i];
            i ++;
        }
        // Read presents a deletion
        else {
            seq[k + al->start] = al->reference[j];
            j ++;
        }
    }
    seq[k + al->start] = '\0';
    return seq;
}

void dump ( aligner_t * al ) {
    int i, j;
    int ** M = al->nw;
    for ( i = 0; i < al->read_len; i ++ ) {
        for ( j = 0; j < al->ref_len; j ++ )
            printf ( "%d\t", M[i][j] );
        printf ( "\n" );
    }
    for ( i = 0; i < al->read_len; i ++ ) {
        for ( j = 0; j < al->ref_len; j ++ )
            printf ( "%c ", al->op[i][j] );
        printf ( "\n" );
    }
}

void al_destroy ( aligner_t * al ){
    for ( int i = 0; i < al->read_len; i ++ ){
        free ( al->nw[i] );
        free ( al->op[i] );
    }
    free ( al->nw );
    free ( al->op );
    free ( al->alignment );
    free ( al );
}
