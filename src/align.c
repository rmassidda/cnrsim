#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "align.h"

aligner_t * al_init ( aligner_t * aligner, char * reference, int end, char * read ){
    aligner_t * al = realloc ( aligner , sizeof ( aligner_t ) );
    if ( al == NULL ){
        return NULL;
    }

    // Sequences
    al->reference = reference;
    al->read = read;
    al->ref_len = end + 1;
    al->read_len = strlen ( read ) + 1;
    
    // Matrix
    if ( aligner == NULL ){
        al->nw = NULL;
        al->op = NULL;
        al->alignment = NULL;
    }

    // Allocate matrix pointers
    int (*nw)[al->ref_len] = realloc( al->nw, sizeof ( int[al->read_len][al->ref_len] ) );
    char (*op)[al->ref_len] = realloc( al->op, sizeof ( char[al->read_len][al->ref_len] ) );
    al->nw = (void *)nw;
    al->op = (void *)op;

    if ( nw == NULL || op == NULL ){
        return NULL;
    }
    
    // Initial conditions
    for ( int i = 0; i < al->read_len; i ++ ) {
        // Initial conditions
        for ( int j = 0; j < al->ref_len; j ++ ) {
            nw[i][j] = 0;
            op[i][j] = '!';
        }
        // GAP in the first column
        nw[i][0] = i * GAP;
        op[i][0] = 'T';
    }
    
    // Alignement result
    al->start = 0;
    __align ( al );
    return al;
}

void __align( aligner_t * al ) {
    int (*nw)[al->ref_len] = al->nw;
    char (*op)[al->ref_len] = al->op;
    int left, top;
    for ( int i = 1; i < al->read_len; i ++ ) {
        for ( int j = 1; j < al->ref_len; j ++ ) {
            left = nw[i][j-1] + GAP;
            top = nw[i-1][j] + GAP;
            // Search for max neighbors
            if ( al->reference[j-1] == al->read[i-1]){
                nw[i][j] = nw[i - 1][j - 1] + MATCH;
                op[i][j] = '=';
            }
            else {
                nw[i][j] = nw[i - 1][j - 1] + MISMATCH;
                op[i][j] = '!';
            }
            if ( top > nw[i][j] ){
                nw[i][j] = top;
                op[i][j] = 'T';
            }
            if ( left > nw[i][j] ){
                nw[i][j] = left;
                op[i][j] = 'L';
            }
        }
    }
}

char * build_alignment ( aligner_t * al ) {
    int (*nw)[al->ref_len] = al->nw;
    char (*op)[al->ref_len] = al->op;
    int i, j, k;
    char * s = malloc ( sizeof ( char ) * ( al->read_len + al->ref_len - 1 ) );

    // Search for the maximum value
    int max_val = nw[al->read_len - 1][0];
    int max_i = al->read_len - 1;
    int max_j = 0;
    for ( j = 0; j < al->ref_len; j ++ ){
        if ( nw[max_i][j] >= max_val ){
           max_j = j;
           max_val = nw[max_i][max_j];
        }
    }

    // Alignment
    i = max_i;
    j = max_j;
    k = 0;
    while ( j != 0 && i != 0 ) {
        s[k] = op[i][j];
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
    int (*nw)[al->ref_len] = al->nw;
    char (*op)[al->ref_len] = al->op;
    for ( i = 0; i < al->read_len; i ++ ) {
        for ( j = 0; j < al->ref_len; j ++ )
            printf ( "%d\t", nw[i][j] );
        printf ( "\n" );
    }
    for ( i = 0; i < al->read_len; i ++ ) {
        for ( j = 0; j < al->ref_len; j ++ )
            printf ( "%c ", op[i][j] );
        printf ( "\n" );
    }
}

void al_destroy ( aligner_t * al ){
    int (*nw)[al->ref_len] = al->nw;
    char (*op)[al->ref_len] = al->op;
    free ( nw );
    free ( op );
    free ( al->alignment );
    free ( al );
}
