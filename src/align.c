#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "align.h"

aligner_t * al_init ( char * s1, char * s2, method_t method ){
    aligner_t * al = malloc ( sizeof ( aligner_t ) );
    if ( al == NULL ){
        return NULL;
    }
    // Sequences
    al->seq[0] = s1;
    al->seq[1] = s2;
    al->len[0] = strlen ( s1 );
    al->len[1] = strlen ( s2 );
    // Matrix
    al->M = NULL;
    al->op = NULL;
    // Method
    al->method = method;
    // Alignement result
    al->alignment = NULL;
    return al;
}

void __initMatrix ( aligner_t * al ) {
    int i, j;
    // Allocate matrix pointers
    al->M = malloc ( al->len[1] * sizeof ( int * ) );
    al->op = malloc ( al->len[1] * sizeof ( char * ) );
    // Allocate rows
    for ( i = 0; i < al->len[1]; i ++ ) {
        al->M[i] = malloc ( al->len[0] * sizeof ( int ) );
        al->op[i] = malloc ( al->len[0] * sizeof ( char ) );
        // Initial conditions
        for ( j = 0; j < al->len[0]; j ++ ) {
            al->M[i][j] = 0;
            al->op[i][j] = '!';
        }
    }
    // Global alignment
    if ( al->method == GLOBAL ) {
        for ( i = 0; i < al->len[0]; i ++ ) {
            al->M[0][i] = i * GAP;
            al->op[0][i] = '1';
        }
        for ( j = 0; j < al->len[1]; j ++ ) {
            al->M[j][0] = j * GAP;
            al->op[j][0] = '1';
        }
        al->op[0][0] = '=';
    }
}

int similarity ( char token1, char token2 ) {
    if ( token1 == token2 ){
        return MATCH;
    }
    return MISMATCH;
}

void align ( aligner_t * al ) {
    int i, j, g1, g2, s;
    int ** M;
    char simb;
    __initMatrix ( al );
    M = al->M;
    for ( j = 1; j < al->len[1]; j ++ ) {
        for ( i = 1; i < al->len[0]; i ++ ) {
            g1 = M[j - 1][i] + GAP;
            g2 = M[j][i - 1] + GAP;
            s = similarity ( al->seq[0][i - 1], al->seq[1][j - 1] );
            if ( s > 0 ) simb = '=';
            else simb = '!';
            /*  Sussume self.M[j][i] = max ( max (g1, g2), self.M[j-1][i-1] + s)  */
            M[j][i] = M[j - 1][i - 1] + s;
            if ( g1 > M[j][i] ) M[j][i] = g1;
            if ( g2 > M[j][i] ) M[j][i] = g2;
            if ( M[j][i] == g2 ) al->op[j][i] = '2';
            if ( M[j][i] == g1 ) al->op[j][i] = '1';
            if ( ( M[j][i] != g1 ) && ( M[j][i] != g2 ) ) al->op[j][i] = simb;
        }
    }
}

void __getLastScore ( aligner_t * al, int * lastC ) {
    int mj, mi, maxval, i, j;
    int ** M = al->M;

    if ( al->method == GLOBAL ) {
        lastC[0] = al->len[1] - 1;
        lastC[1] = al->len[0] - 1;
        return;
    }
    mj = al->len[1] - 1;
    mi = al->len[0] - 1;
    maxval = M[mj][mi];
    /* search the highest value in the right border of the matrix  */
    for ( j = 0; j < al->len[1]; j ++ ) {
        if ( M[j][al->len[0] - 1] > maxval ) {
            mj = j;
            maxval = M[j][al->len[0] - 1];
        }
    }
    /* search the highest value in the bottom of the matrix  */
    for ( i = 0; i < al->len[0]; i ++ ) {
        if ( M[al->len[1] - 1][i] > maxval ) {
            mi = i;
            mj = al->len[1] - 1;
            maxval = M[al->len[1] - 1][i];
        }
    }
    lastC[0] = mj;
    lastC[1] = mi;
}

char * build_alignment ( aligner_t * al ) {
    int lastC[2];
    int i, j, k;
    char * s;
    k = 0;
    s = malloc ( sizeof ( char ) * ( al->len[1] + al->len[0] + 1 ) );
    __getLastScore ( al, lastC );
    /* for local alignment only one loop will run  */
    i = al->len[0] - 1;
    while ( i > lastC[1] ) {
        s[k] = '2';
        i --;
        k ++;
    }
    j = al->len[1] - 1;
    while ( j > lastC[0] ) {
        s[k] = '1';
        j --;
        k ++;
    }
    while ( j != 0 && i != 0 ) {
        s[k] = al->op[j][i];
        if ( s[k] == '2' ) i --;
        else {
            if ( s[k] == '1' ) j --;
            else {
                j --;
                i --;
            }
        }
        k ++;
    }
    /* for local alignment only one loop will run  */
    while ( i > 0 ) {
        s[k] = '2';
        i --;
        k ++;
    }
    while ( j > 0 ) {
        s[k] = '1';
        j --;
        k ++;
    }
    s[k] = '\0';
    al->alignment = malloc ( sizeof ( char ) * ( k + 1 ) );
    for ( i = 0, k --; k >= 0; i++, k -- ) al->alignment[i] = s[k];
    al->alignment[i] = '\0';
    free ( s );
    return al->alignment;
}

char * alignment ( aligner_t * al, char * alignmentString ) {
    char * alignSeq;
    int c1, c2, k;
    alignSeq = malloc ( sizeof ( char ) * ( al->len[1] + al->len[0] + 1 ) );
    //alignSeq = (char *) malloc ((strlen (alignmentString) + 1) * sizeof (char));
    c1 = 0;
    c2 = 0;
    for ( k = 0; k < strlen ( alignmentString ); k ++ ) {
        if ( alignmentString[k] == '!' || alignmentString[k] == '=' ) {
            /* If one is N is ignered  XXX - only if !  */
            if ( al->seq[0][c1] == 'N' ) alignSeq[k] = al->seq[1][c2];
            else {
                if ( al->seq[1][c2] == 'N' ) alignSeq[k] = al->seq[0][c1];
                else alignSeq[k] = al->seq[0][c1];
            }
            c1 ++;
            c2 ++;
        }
        if ( alignmentString[k] == '1' ) {
            alignSeq[k] = al->seq[1][c2];
            c2 ++;
        }
        if ( alignmentString[k] == '2' ) {
            alignSeq[k] = al->seq[0][c1];
            c1 ++;
        }
    }
    alignSeq[k] = '\0';
    return alignSeq;
}

void dump ( aligner_t * al ) {
    int i, j;
    int ** M = al->M;
    for ( i = 0; i < al->len[1]; i ++ ) {
        for ( j = 0; j < al->len[0]; j ++ )
            printf ( "%d\t", M[i][j] );
        printf ( "\n" );
    }
    for ( i = 0; i < al->len[1]; i ++ ) {
        for ( j = 0; j < al->len[0]; j ++ )
            printf ( "%c ", al->op[i][j] );
        printf ( "\n" );
    }
}

void al_destroy ( aligner_t * al ){
    for ( int i = 0; i < al->len[1]; i ++ ){
        free ( al->M[i] );
        free ( al->op[i] );
    }
    free ( al->M );
    free ( al->op );
    free ( al->alignment );
    free ( al );
}
