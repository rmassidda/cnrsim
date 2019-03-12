/*
 * CNRSIM
 * wrapper.h
 * Library that wrappes the HTSLIB parser for
 * VCF file, and the CNRSIM parser for user
 * defined variants.
 *
 * @author Riccardo Massidda
 */

#ifndef WRAPPER
#define WRAPPER

#include <stdbool.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "user_variation.h"

typedef struct wrapper_t wrapper_t;

enum parser {
    NO = 0,
    VCF = 1,
    UDV = 2,
    BOTH = 3
};

struct wrapper_t {
    // VCF
    bcf_srs_t * sr;
    bcf_hdr_t * hdr;
    bcf1_t * vcf_line;
    // UDV
    variation_set_t * udv;
    variation_t * udv_line;
    /*
     * Status variables.
     * First bit for VCF
     * Second bit for UDV
     */
    int present; // which file has been correctly opened
    int seek; // which file contains the requested region
    int used; // which file needs to update is buffer
    // Public data
    char * region; // current region
    int pos; // position of the variation
    char * ref; // reference
    char ** alt; // alternatives
    int * alt_index; // alternative chosen for each allele
    int ploidy;
};

/*
 * Initialize a structure containing
 * the data used to arbitrate VCF
 * and UDV parsers.
 * At least one of the two file has to
 * be present.
 *
 * @param vcf_filename path to the vcf file
 * @param udv_filename path to the user defined variants file
 * @returns pointer to the initialized wrapper, NULL otherwise
 */
wrapper_t * wr_init ( char * udv_filename, char * vcf_filename, int ploidy );

/*
 * Sets the position of the readers
 * to the start point of the region.
 *
 * @param w pointer to the variation set
 * @param region label of the desired region
 * @returns 0 if no reader set
 */
int wr_seek ( wrapper_t * w, char * label );

/*
 * Checks if there are any lines available
 * in the current region.
 *
 * @param w pointer to the wrapper
 * @returns 0 if the region ended
 */
int wr_region ( wrapper_t * w );

/*
 * Select which record use, and updates
 * the structure.
 *
 * @param w pointer to the wrapper
 */
bool wr_update_wrapper ( wrapper_t * w );

/*
 * Free allocated memory
 *
 * @param w pointer to the reader
 */
void wr_destroy ( wrapper_t * w );

#endif
