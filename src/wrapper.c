/*
 * CNRSIM
 * wrapper.c
 * Library that wrappes the HTSLIB parser for
 * VCF file, and the CNRSIM parser for user
 * defined variants
 *
 * @author Riccardo Massidda
 */
#include <time.h>
#include "wrapper.h"
#include "parse_frequency.h"

// Global variables
double * p = NULL;

wrapper_t * wr_init ( char * vcf_filename, char * udv_filename ) {
    wrapper_t * w;
    // At least a filename is required
    if ( vcf_filename == NULL && udv_filename == NULL ) {
        return NULL;
    }
    // Allocate memory
    w = malloc ( sizeof ( wrapper_t ) );
    if ( w == NULL ) {
        return NULL;
    }

    w->present = 0;
    w->seek = 0;
    w->used = 0;
    w->udv = NULL;
    w->udv_line = NULL;

    // Random number generator
    srand ( time ( NULL ) );

    // VCF
    if ( vcf_filename != NULL ) {
        w->sr = bcf_sr_init ();
        // Index required
        bcf_sr_set_opt ( w->sr, BCF_SR_REQUIRE_IDX );
        // Reader file link
        if ( bcf_sr_add_reader ( w->sr, vcf_filename ) != 1 ) {
            return NULL;
        }

        /*
         * If the file is indexed the name
         * and the position of the regions
         * are loaded into the ws
         * structure.
         */
        if ( w->sr->regions == NULL ) {
            // File not indexed, creation of the index
            // 14 is a suggested value by HTSLIB documentation
            int res = bcf_index_build ( vcf_filename, 14 );
            if ( res == 0 ) {
                // File is now indexed
                // Relink of the file with the w
                bcf_sr_remove_reader ( w->sr, 0 );
                if ( bcf_sr_add_reader ( w->sr, vcf_filename ) != 1 ) {
                    return NULL;
                }
            } else {
                perror ( "File is not indexable.\nTry first compressing it with gzip.\n" );
                return NULL;
            }
        }

        /*
         * Load of the header, needed to parse
         * INFO values in the VCF.
         */
        w->hdr = bcf_sr_get_header ( w->sr, 0 );
        w->vcf_line = NULL;

        w->present += VCF;
    }

    // UDV
    if ( udv_filename != NULL ) {
        // Initalize UDV structure
        w->udv = udv_init ( udv_filename );
        if ( w->udv == NULL ) {
            return NULL;
        }
        w->present += UDV;
    }

    return w;
}

int wr_seek ( wrapper_t * w, char * label ) {
    // Initialize
    w->seek = 0;
    // New line must be read
    w->used = BOTH;
    // Update region
    w->region = label;
    // Update region
    w->region = label;
    if ( w->present & VCF ) {
        if ( bcf_sr_seek ( w->sr, label, 0 ) == 0 ) {
            w->seek += VCF;
        }
    }
    if ( w->present & UDV ) {
        if ( udv_seek ( w->udv, label ) ) {
            w->seek += UDV;
        }
    }
    return w->seek;
}

int wr_region ( wrapper_t * w ) {
    // Check if the buffered line has to be updated
    if ( w->seek & UDV && w->used & UDV ) {
        if ( udv_next_line ( w->udv ) ) {
            w->udv_line = udv_get_line ( w->udv );
            w->used -= UDV;
        }
        // Region ended
        else {
            w->seek -= UDV;
        }
    }

    if ( w->seek & VCF && w->used & VCF ) {
        if ( bcf_sr_next_line ( w->sr ) ) {
            w->vcf_line = bcf_sr_get_line ( w->sr, 0 );
            // bcf_sr_next_line doesn't return false on region change
            int ret = strncmp (
                          w->region,
                          bcf_hdr_id2name ( w->hdr, w->vcf_line->rid ),
                          strlen ( w->region )
                      );
            if ( ret != 0 ) {
                w->seek -= VCF;
            } else {
                w->used -= VCF;
            }
        } else {
            w->seek -= VCF;
        }
    }

    return w->seek;
}

bool _detect_collision ( wrapper_t * w ) {
    // Coordinates of the reference
    int vcf_start = w->vcf_line->pos;
    int udv_start = w->udv_line->pos;
    int vcf_end = w->vcf_line->pos + strlen ( w->vcf_line->d.allele[0] );
    int udv_end = w->udv_line->pos + strlen ( w->udv_line->ref );
    // Collision
    return ( ( vcf_end >= udv_start ) && ( vcf_start <= udv_end ) );
}

bool _vcf2wrapper ( wrapper_t * w ) {
    double outcome;
    double threshold;
    //  Allelic frequency parser
    float * af = NULL;
    int af_size = 0;
    int af_ret;
    char * freq = NULL;
    int freq_size = 0;
    int freq_ret;

    // Unpack line
    if ( bcf_unpack ( w->vcf_line, BCF_UN_STR ) != 0 ) {
        perror ( "Unpack error" );
        return false;
    }

    w->pos = w->vcf_line->pos;
    w->ref = w->vcf_line->d.allele[0];
    // Even the reference is a possible alternative!
    w->alt = w->vcf_line->d.allele;
    // Allelic frequency as defined by VCF
    af_ret = bcf_get_info_float ( w->hdr, w->vcf_line, "AF", af, &af_size );
    // Allelic frequency as defined by dbSNP
    freq_ret = bcf_get_info_string ( w->hdr, w->vcf_line, "FREQ", &freq, &freq_size );
    // Parse results
    if ( af_ret >= 0  ) {
        p = parse_af ( w->vcf_line->n_allele, af, p );
        free ( af );
    } else if ( freq_ret >= 0 ) {
        p = parse_db_snp_freq ( w->vcf_line->n_allele, freq, p );
        free ( freq );
    } else {
        p = linear ( w->vcf_line->n_allele, p );
    }

    for ( int i = 0; i < ALL_N; i++ ) {
        // Random decision about the alternatives
        outcome = ( double ) rand() / RAND_MAX;
        threshold = 0;
        for ( int j = 0; j < w->vcf_line->n_allele; j++ ) {
            if ( threshold <= outcome && outcome < threshold + p[j] ) {
                w->alt_index[i] = j;
                break;
            } else {
                threshold += p[j];
            }
        }
    }
    return true;
}

bool _udv2wrapper ( wrapper_t * w ) {
    w->pos = w->udv_line->pos;
    w->ref = w->udv_line->ref;
    w->alt = w->udv_line->all;
    // Allele alternatives already assigned
    for ( int i = 0; i < ALL_N; i++ ) {
        w->alt_index[i] = i;
    }
    return true;
}

bool wr_update_wrapper ( wrapper_t * w ) {
    /*
     * The function MUST be used
     * after wr_region, so at least
     * one line has to be usable.
     */
    assert ( w->used != BOTH );

    if ( w->used & VCF ) {
        // VCF isn't usable
        w->used += UDV;
        return _udv2wrapper ( w );
    } else if ( w->used & UDV ) {
        // UDV isn't usable
        w->used += VCF;
        return _vcf2wrapper ( w );
    } else {
        if ( w->vcf_line->pos < ( w->udv_line->pos + strlen ( w->udv_line->ref ) ) ) {
            w->used += VCF;
            // If there is a collision
            // line must be ignored
            return _vcf2wrapper ( w ) && ! ( _detect_collision ( w ) );
        } else {
            // Threshold reached
            w->used += UDV;
            return _udv2wrapper ( w );
        }
    }

    return true;
}

void wr_destroy ( wrapper_t * w ) {
    // VCF
    bcf_sr_destroy ( w->sr );
    // UDV
    if ( w->udv != NULL )
        udv_destroy ( w->udv );
    // Wrapper
    free ( w );
    free ( p );
    return;
}
