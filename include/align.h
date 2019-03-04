#ifndef ALIGN
#define ALIGN

typedef struct aligner_t aligner_t;

enum error{
    GAP = -1,
    MATCH = 1,
    MISMATCH = 0
};

struct aligner_t {
    char * reference; // reference
    char * read; // read
    int ref_len; // length of the reference (zero-included)
    int read_len; // length of the read ( zero-included )
    void * nw; // Needleman–Wunsch matrix
    void * op; // matrix operations
    char * alignment; // result of the alignment
    int start; // start position of the alignment
};

/*
 * Initalize the structure used by the aligner
 *
 * @param s1 first sequence
 * @param s2 second sequence
 * @param aligner previously allocated structure to be reused
 * @ret initialized structure, NULL if error
 */
aligner_t * al_init ( aligner_t * aligner_t, char * reference, int end, char * read );

/*
 * Populates the matrix
 *
 * @param al pointer to the aligner
 */
void __align ( aligner_t * al ) ;

/*
 * Computes the alignment string
 *
 * @param al pointer to the aligner
 * @ret alignment string
 */
char * build_alignment ( aligner_t * al ) ;

char * alignment ( aligner_t * al ) ;

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
