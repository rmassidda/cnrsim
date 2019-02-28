/*
 * CNRSIM
 * error_profiler.c
 * Given a BAM containing real reads
 * analyzes the distribution and the
 * entity of the errors.
 *
 * @author Riccardo Massidda
 */
#include <stdlib.h>
#include <stdio.h>
#include "align.h"

int main ( int argc, char ** argv ) {
    char * s1 = argv[1];
    char * s2 = argv[2];
    aligner_t * al = NULL;
    printf ( "read\t%s\n", s2 );
    printf ( "seq\t%s\n", s1 );

    al = al_init ( s1, s2 );
    char * str = build_alignment ( al );
    printf ( "al\t%s\n", alignment ( al ) );
    printf ( "al\t%s\n", str );

    al_destroy ( al );
    return 0;
}
