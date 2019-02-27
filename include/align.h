#ifndef ALIGN
#define ALIGN

typedef struct aligner_t aligner_t;
typedef enum alignment method_t;

enum alignment{
    GLOBAL,
    LOCAL
};

enum error{
    GAP = -1,
    MATCH = 1,
    MISMATCH = 0
}

struct aligner_t {
    char * seq[2]; // sequences to be aligned
    int len[2]; // length of the sequences
    int ** M; // matrix M
    char ** op; // matrix operations
    enum alignment method; // method of alignment
    char * aligned; // result of the alignment
};

/*
 * Initalize the structure used by the aligner
 *
 * @param s1 first sequence
 * @param s2 second sequence
 * @param method method of alignement
 * @ret initialized structure, NULL if error
 */
aligner_t * al_init ( char * s1, char * s2, method_t method );

void __initMatrix ( aligner_t * al ) ;

int similarity ( char token1, char token2 ) ;

void align ( aligner_t * al ) ;

void __getLastScore ( aligner_t * al, int * lastC ) ;

char * buildAlignment ( aligner_t * al ) ;

char * alignment ( aligner_t * al, char * alignmentString ) ;

void dump ( aligner_t * al );

/*
 * Frees the memory
 *
 * @param al pointer to the aligner
 */
void al_destroy ( aligner_t * al );

#endif
