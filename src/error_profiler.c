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
#include <edlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <time.h>
#include "allele.h"
#include "model.h"
#include "stats.h"
#include "tandem.h"
#include "translate_notation.h"

#define MAX_MOTIF 6

// Init kseq structure
KSEQ_INIT ( gzFile, gzread );

void usage ( char * name ) {
    fprintf ( stderr, "Usage: %s [-d dictionary] [-t] [-v] [-s] bam_file fasta_file [allele_file ...]\n", name );
}

void dump_read ( char * ref, unsigned char * alignment, int alg_len, char * read, uint8_t * quality ) {
    // Quality
    for ( int z = 0; z < alg_len; z ++ ) {
        switch ( alignment[z] ) {
        case 2:
            printf ( " " );
            break;
        default:
            printf ( "%c", ( *quality + 33 ) );
            quality ++;
            break;
        }
    }
    printf ( "\n" );

    // Read
    for ( int z = 0; z < alg_len; z ++ ) {
        switch ( alignment[z] ) {
        case 2:
            printf ( " " );
            break;
        default:
            printf ( "%c", *read );
            read ++;
            break;
        }
    }
    printf ( "\n" );

    // Alignment
    for ( int z = 0; z < alg_len; z ++ ) {
        switch ( alignment[z] ) {
        case 0:
            printf ( "=" );
            break;
        case 1:
            printf ( "I" );
            break;
        case 2:
            printf ( "D" );
            break;
        case 3:
            printf ( "!" );
            break;
        }
    }
    printf ( "\n" );

    // Reference
    for ( int z = 0; z < alg_len; z ++ ) {
        switch ( alignment[z] ) {
        case 1:
            printf ( " " );
            break;
        default:
            printf ( "%c", *ref );
            ref ++;
            break;
        }
    }
    printf ( "\n" );
}

int main ( int argc, char ** argv ) {
    // Parser
    int opt;
    char * dictionary = NULL;
    char * bam_fn;
    bool verbose = false;
    bool silent = false;
    int ploidy;
    // FASTA
    gzFile * fasta;
    kseq_t ** seq;
    kseq_t * curr_seq;
    char * align;
    allele_t ** allele;
    bool last = false;
    // Tandem repeats
    int tandem = 16;
    tandem_set_t ** trs;
    // Insert size granularity
    int size_granularity = 1024;
    // BAM
    htsFile * fp;
    bam_hdr_t * hdr;
    bam1_t * line;
    hts_idx_t * index;
    hts_itr_t * itr;
    // Alias dictionary
    region_index_t * alias_index = NULL;
    char * alias = NULL;
    // Aligner
    EdlibAlignResult  * edlib_alg;
    EdlibAlignConfig config;
    EdlibEqualityPair additionalEqualities[4] = {
        {'A', 'a'},
        {'C', 'c'},
        {'G', 'g'},
        {'T', 't'}
    };
    char * read = NULL;
    int pos;
    int len;
    uint8_t * read_seq;
    uint8_t * qual;
    int insert_size;
    int max_insert_size = 4096;
    int orientation;
    int flank_1;
    int flank_2;
    int start;
    int end;
    // Statistics
    model_t * model;
    stats_t * curr_stats;
    int min_score;
    int min_start = 0;
    int min_index = 0;
    long long int read_counter = 0;
    long long int skipped = 0;
    int density = 1;

    // Init pseudorandom generator
    srand ( time ( NULL ) );

    while ( ( opt = getopt ( argc, argv, "svm:i:t:d:p:" ) ) != -1 ) {
        switch ( opt ) {
        case 's':
            silent = true;
            break;
        case 'v':
            verbose = true;
            break;
        case 't':
            tandem = atoi ( optarg );
            break;
        case 'i':
            size_granularity = atoi ( optarg );
            break;
        case 'm':
            max_insert_size = atoi ( optarg );
            break;
        case 'd':
            dictionary = optarg;
            break;
        case 'p':
            density = atoi ( optarg );
            break;
        case '?':
            if ( optopt == 'p' || optopt == 'd' || optopt == 'm' || optopt == 'i' || optopt == 't' )
                fprintf ( stderr, "Option -%c requires an argument.\n", optopt );
            else if ( isprint ( optopt ) )
                fprintf ( stderr, "Unknown option `-%c'.\n", optopt );
            else
                fprintf ( stderr, "Unknown option character `\\x%x'.\n", optopt );
            exit ( EXIT_FAILURE );
        default:
            usage ( argv[0] );
            exit ( EXIT_FAILURE );
        }
    }

    if ( dictionary != NULL ) {
        // Alias dictionary
        alias_index = tr_init ( dictionary );
        if ( alias_index == NULL ) {
            perror ( "Can't load alias dictionary" );
            exit ( EXIT_FAILURE );
        }
    }

    // Non optional arguments
    if ( argc - optind < 2 ) {
        usage ( argv[0] );
        exit ( EXIT_FAILURE );
    }

    // BAM file
    bam_fn = argv[optind++];
    fp = hts_open ( bam_fn, "r" );
    hdr = sam_hdr_read ( fp );
    line = bam_init1 ();
    // BAM index
    index = bam_index_load ( bam_fn );
    if ( index == NULL ) {
        // File not indexed
        perror ( "Can't load BAM index" );
        exit ( EXIT_FAILURE );
    }

    // FASTA files
    ploidy = argc - optind;

    // Malloc of the structures
    fasta = malloc ( sizeof ( gzFile ) * ploidy );
    seq = malloc ( sizeof ( kseq_t * ) * ploidy );
    allele = malloc ( sizeof ( allele_t * ) * ploidy );
    edlib_alg = malloc ( sizeof ( EdlibAlignResult ) * ploidy );
    trs = malloc ( sizeof ( tandem_set_t * ) * ploidy );

    // Init sequences
    for ( int i = 0; i < ploidy; i ++ ) {
        // Read of the allele
        fasta[i] = gzopen ( argv[optind], "r" );
        if ( fasta[i] == NULL ) {
            fprintf ( stderr, "File %s not found.\n", argv[optind] );
            exit ( EXIT_FAILURE );
        }
        seq[i] = kseq_init ( fasta[i] );
        trs[i] = NULL;
        allele[i] = NULL;
        optind ++;
    }

    // Edlib configuration
    config = edlibNewAlignConfig ( -1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4 );

    model = model_init ( MAX_MOTIF, tandem, max_insert_size, size_granularity );

    // While there are sequences to read in the FASTA file
    while ( ! last ) {
        // Next sequence
        for ( int i = 0; i < ploidy; i ++ ) {
            int x = kseq_read ( seq[i] );
            if ( x >= 0 || x == -2 ) {
                // Alignment of the sequence
                align = ( seq[i]->qual.l > 0 ) ? seq[i]->qual.s : NULL;
                // Update allele
                allele[i] = allele_point (
                        seq[i]->seq.l,
                        seq[i]->seq.s,
                        align,                        
                        allele[i]
                        );
                if ( tandem != 0 ) {
                    trs[i] = tandem_set_init ( seq[i]->seq.l, model->max_motif, tandem, trs[i] );
                    trs[i] = tandem_set_analyze ( seq[i]->seq.s, seq[i]->seq.l, trs[i] );
                }
            }
            else{
                last = true;
            }
        }
        
        if ( last ) { break; }
        // Seek on the BAM
        alias = NULL;
        if ( dictionary != NULL ) {
            alias = tr_translate ( alias_index, (*seq)->name.s );
        }
        if ( alias == NULL ){
            alias = (*seq)->name.s;
        }
        itr = bam_itr_querys ( index, hdr, alias );
        if ( itr != NULL ) {
            while ( bam_itr_next ( fp, itr, line ) > 0 ) {
                read_counter++;
                double done = 100.0 * ( read_counter - skipped ) / read_counter;
                fprintf ( stderr, "READ> %lld : %s ( %.3f done ) actually in %s  #\r", read_counter, bam_get_qname ( line ), done, alias );
                if ( read_counter % density != 0 ) {
                  skipped ++;
                  continue;
                }

                if ( (( line->core.flag & 1 ) && ( line->core.flag & 2 )) || line->core.flag == 0 || line->core.flag == 16 ){
                  // Read information
                  pos = line->core.pos;
                  len = line->core.l_qseq;
                  read_seq = bam_get_seq ( line ); // Read nucleotides
                  qual = bam_get_qual ( line ); // Quality score
                  curr_stats = ( line->core.flag & 64 ) ? model->single : model->pair;
                  curr_stats = ( line->core.flag & 128 ) ? model->pair : model->single;
                  if ( line->core.flag == 0 || line->core.flag == 16 ) {
                    curr_stats = model->single;
                  }
                }
                else{
                    skipped++;
                    continue;
                }

                if ( curr_stats == model->single && ( line->core.tid == line->core.mtid || line->core.mtid == -1 ) ){
                  insert_size = line->core.mpos - ( pos + len );
                  insert_size = ( insert_size < 0 ) ? -insert_size : insert_size;
                  insert_size = ( insert_size < max_insert_size ) ? insert_size : max_insert_size; 
                  insert_size = ( insert_size * size_granularity ) / max_insert_size;
                  insert_size = ( line->core.mtid == -1 ) ? 0 : insert_size;
                  source_update ( NULL, 0, 0, insert_size, model->insert_size );
                  orientation = ( line->core.flag & 48 ) >> 4;
                  source_update ( NULL, 0, 0, orientation, model->orientation );
                }

                // Interval of the reference
                flank_1 = floor ( log ( 2 * len ) / log ( 2 ) );
                flank_2 = flank_1;

                // Read string
                read = realloc ( read, sizeof ( char ) * ( len + 1 ) );
                int i;
                for ( i = 0; i < len; i++ ) {
                    read[i] = seq_nt16_str[ bam_seqi ( read_seq, i ) ];
                }
                read[i] = 0;

                for ( i = 0; i < ploidy; i ++ ) {
                    curr_seq = seq[i];
                    // Seek on the allele
                    allele_seek ( pos, true, allele[i] );

                    // Flanking regions
                    start = allele[i]->pos;
                    if ( start - flank_1 < 0 ) {
                        start = 0;
                    } else {
                        start -= flank_1;
                    }

                    end = allele[i]->pos + len;
                    if ( end + flank_2 >= curr_seq->seq.l ) {
                        end = curr_seq->seq.l - 1;
                    } else {
                        end += flank_2;
                    }

                    // Align
                    edlib_alg[i] = edlibAlign (
                                       read,
                                       len,
                                       &curr_seq->seq.s[start],
                                       end - start,
                                       config );

                    // Select best alignment
                    if ( i == 0 || edlib_alg[i].editDistance < min_score ) {
                        min_index = i;
                        min_start = start + edlib_alg[i].startLocations[0];
                        min_score = edlib_alg[i].editDistance;
                    }

                    // Note: edlib c(GAP) = c(MM)
                    if ( ( line->core.flag == 147 || line->core.flag == 83 ) && edlib_alg[i].alignment[edlib_alg[i].alignmentLength-1] == 1 )
                        edlib_alg[i].alignment[edlib_alg[i].alignmentLength-1] = 3;

                    if ( verbose ) {
                        printf ( "Sequence no.%d %d->%ld\n", i, pos, allele[i]->pos );
                        dump_read (
                            &curr_seq->seq.s[start + edlib_alg[i].startLocations[0]],
                            edlib_alg[i].alignment,
                            edlib_alg[i].alignmentLength,
                            read,
                            qual );
                    }
                }

                if ( tandem != 0 ){
                    // Tandem index
                    int tin = trs[min_index]->i;
                    // Number of tandems
                    int n = trs[min_index]->n;
                    // Tandems set
                    tandem_t * set = trs[min_index]->set;

                    // Tandem not in the reads
                    while ( tin < n && set[tin].pos < min_start ){
                        tin ++;
                    }
                    // Reads are sorted, so these tandems will never be used
                    trs[min_index]->i = tin;

                    // Possibile tandems
                    int read_end = min_start + strlen ( read );
                    while ( tin < n && set[tin].pos < read_end ){
                        // Learning example
                        int motif = set[tin].pat;
                        unsigned char in = set[tin].rep;
                        unsigned char out = 0;
                        // Find position
                        int read_pos = -1;
                        int ref_pos = min_start - 1;
                        int alg = 0;
                        unsigned char * alg_str = edlib_alg[min_index].alignment;
                        int alg_len = edlib_alg[min_index].alignmentLength;
                        // Adjust ref_pos and ref_pos
                        while ( alg < alg_len && ( ref_pos < 0 || ref_pos < set[tin].pos ) ) {
                            switch ( alg_str[alg] ){
                                case 1:
                                    read_pos ++;
                                    break;
                                case 2:
                                    ref_pos ++;
                                    break;
                                default:
                                    read_pos ++;
                                    ref_pos ++;
                            }
                            alg ++;
                        }
                        if ( read_pos + (out + 1) * motif < len ){
                            // Compare pattern
                            if ( strncasecmp ( &seq[min_index]->seq.s[set[tin].pos], &read[read_pos], motif ) == 0 ){
                                out = 1;
                                // Compare right patterns
                                while ( read_pos + (out + 1) * motif < len ){
                                    if ( strncasecmp ( &read[read_pos], &read[read_pos + out * motif], motif ) == 0 ){
                                        out ++;
                                    }
                                    else{
                                        break;
                                    }
                                }
                            }
                        }
                        if ( verbose ) {
                            alg = 0;
                            // Print of the reference tandem
                            int z = min_start - set[tin].pos;
                            while ( z < motif * in ){
                                if ( alg_str[alg] == 1 ){
                                    printf ( " " );
                                }
                                else if ( z < 0 ){
                                    printf ( " " );
                                    z ++;
                                }
                                else {
                                    printf ( "%c", seq[min_index]->seq.s[set[tin].pos+z] );
                                    z ++;
                                }
                                alg ++;
                            }
                            printf ( "\n" );
                            alg = 0;
                            // Print of the read tandem
                            z = - read_pos;
                            while ( z < motif * out ){
                                if ( alg_str[alg] == 2 ){
                                    printf ( " " );
                                }
                                else if ( z < 0 ){
                                    printf ( " " );
                                    z ++;
                                }
                                else {
                                    printf ( "%c", read[read_pos + z] );
                                    z ++;
                                }
                                alg ++;
                            }
                            printf ( "\n" );
                        }
                        // Update tandem statistics
                        if ( in < tandem && out < tandem ) {
                            source_update ( &in, 1, motif, out, model->amplification );
                        }
                        tin ++;
                    }
                }

                stats_update (
                    edlib_alg[min_index].alignment,
                    edlib_alg[min_index].alignmentLength,
                    read,
                    &allele[min_index]->sequence[min_start],
                    qual,
                    curr_stats                    
                );

                for ( int i = 0; i < ploidy; i ++ ) {
                    edlibFreeAlignResult ( edlib_alg[i] );
                }
            }
        }
        else{
            fprintf ( stderr, "%s not found.\n", (*seq)->name.s );
        }
    }

    // Dump statistics
    if ( !silent ) {
        model_dump ( stdout, model );
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i ++ ) {
        gzclose ( fasta[i] );
        kseq_destroy( seq[i] );
        tandem_set_destroy ( trs[i] );
        // Free and not destroy because
        // the sequence is externally allocated
        free ( allele[i] );
    }
    free ( fasta );
    free ( seq );
    free ( trs );
    free ( allele );
    free ( edlib_alg );
    bam_destroy1 ( line );
    bam_hdr_destroy ( hdr );
    bam_itr_destroy ( itr );
    hts_idx_destroy ( index );
    model_destroy ( model );
    sam_close ( fp );
    free ( read );
    tr_destroy ( alias_index );
    return 0;
}
