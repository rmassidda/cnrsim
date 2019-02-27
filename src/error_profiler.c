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
    char * s1 = "GATTACA";
    char * s2 = "TGACTATA";
    aligner_t * al = NULL;
    printf ( "%s\n", s1 );
    printf ( "%s\n", s2 );

    al = al_init ( s1, s2, GLOBAL );
    align ( al );    
    char * t = buildAlignment ( al );
    printf ( "%s\n", t );
    al_destroy ( al );
    return 0;
}
