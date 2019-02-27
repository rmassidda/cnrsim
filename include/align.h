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
};

struct aligner_t {
    char * seq[2]; // sequences to be aligned
    int len[2]; // length of the sequences (zero-included)
    int ** M; // matrix M
    char ** op; // matrix operations
    enum alignment method; // method of alignment
    char * alignment; // result of the alignment
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

/*
 * Checks if two char are equals.
 *
 * @param token1 first char
 * @param token2 second char
 * @ret MATCH if equals, MISMATCH otherwise.
 */
int similarity ( char token1, char token2 ) ;

void align ( aligner_t * al ) ;

void __getLastScore ( aligner_t * al, int * lastC ) ;

/*
 * Alignes the first string to the second.
 *
 * @param al pointer to the aligner
 * @ret alignment string
 */
char * build_alignment ( aligner_t * al ) ;

char * alignment ( aligner_t * al, char * alignmentString ) ;

/* 
 * Prints the matrixes.
 *
 * @param al pointer to the aligner
 */
void dump ( aligner_t * al );

/*
 * Frees the memory
 *
 * @param al pointer to the aligner
 */
void al_destroy ( aligner_t * al );

#endif
