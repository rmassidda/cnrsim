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
#include <htslib/sam.h>
#include <edlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include "allele.h"
#include "fileManager.h"
#include "stats.h"
#include "translate_notation.h"

void usage ( char * name){
    fprintf(stderr, "Usage: %s [-d dictionary] [-a] [-v] [-e] bam_file fasta_file [allele_file ...]\n", name );
}

void dump_read ( char * ref, int ref_len, int start, unsigned char * alignment, int alg_len, char * read , uint8_t * quality ){
    int z;
    char c;
    // Reference
    printf ( "%.*s\n", ref_len, ref );
    // Aligned Read
    for ( z = 0; z < start; z ++ ) printf ( " " );
    int j = 0;
    for ( z = 0; z < alg_len; z++ ){
        if ( alignment[z] == 'X' || alignment[z] == 3 ){
            printf ( "|" );
            j++;
        }
        else if ( alignment[z] == 'I' || alignment[z] == 2 ){
            printf ( "|" );
        } 
        else if ( alignment[z] == 'D' || alignment [z] == 1){
            j++;
        }
        else{
            printf ( "%c", read[j] );
            j++;
        }
    }
    printf ("\n");
    for ( z = 0; z < start; z ++ ) printf ( " " );
    printf("%s\n", read );
    for ( z = 0; z < start; z ++ ) printf ( " " );
    for ( z = 0; z < strlen ( read ); z ++ ) printf ( "%c", quality[z] + 33 );
    printf ("\n");
    
    // Alignment
    for ( z = 0; z < start; z ++ ) printf ( " " );
    for ( z = 0; z < alg_len; z ++ ){
        switch ( alignment[z] ){
            case 0: c = '='; break;
            case 1: c = 'D'; break;
            case 2: c = 'I'; break;
            case 3: c = 'X'; break;
            default: c = alignment[z];
        }
        printf ( "%c", c );
    }
    printf ( "\n\n" );
}

int main ( int argc, char ** argv ) {
    // Parser
    int opt;
    char * dictionary = NULL;
    bool alternative = false;
    char * bam_fn;
    char * fasta_fn = NULL;
    bool verbose = false;
    bool silent = false;
    int ploidy;
    // FASTA
    struct filemanager ** fm;
    struct sequence_t ** seq;
    allele_t ** allele;
    // Pointers
    struct sequence_t * curr_seq;
    // BAM
    htsFile *fp;
    bam_hdr_t *hdr;
    bam1_t *line;
    hts_idx_t *index;
    hts_itr_t *itr;
    // Alias dictionary
    region_index_t * alias_index = NULL;
    char * alias = NULL;
    // Aligner
    EdlibAlignResult  * edlib_alg;
    EdlibAlignConfig config;
    char * read = NULL;
    int pos;
    int len;
    int flank_1;
    int flank_2;
    int start;
    int end;
    // Statistics
    stats_t * stats[2];
    int pair;
    int min_score;
    int min_index = 0;

    while ((opt = getopt(argc, argv, "svd:a")) != -1) {
        switch (opt) {
            case 's':
                silent = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'd':
                dictionary = optarg;
                break;
            case 'a':
                alternative = true;
                break;
            case '?':
                if (optopt == 'd')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                exit ( EXIT_FAILURE );
            default:
                usage ( argv[0] );
                exit(EXIT_FAILURE);
        }
    }

    if ( dictionary != NULL ){
        // Alias dictionary
        alias_index = tr_init ( dictionary );
        if ( alias_index == NULL ){
            perror ( "Can't load alias dictionary" );
            exit ( EXIT_FAILURE );
        }
    }

    // Non optional arguments
    if ( argc - optind < 2 ){
        usage ( argv[0] );
        exit(EXIT_FAILURE);
    }

    // BAM file
    bam_fn = argv[optind++];
    fp = hts_open ( bam_fn, "r" );
    hdr = sam_hdr_read ( fp );
    line = bam_init1 ();
    // BAM index
    index = bam_index_load ( bam_fn );
    if ( index == NULL ){
        // File not indexed
        perror ( "Can't load BAM index" );
        exit ( EXIT_FAILURE );
    }

    // FASTA files
    ploidy = alternative ? ( argc - optind ) : 1;

    // Malloc of the structures
    fm = malloc ( sizeof ( struct filemanager_t * ) * 2 * ploidy );
    seq = malloc ( sizeof ( struct sequence_t * ) * 2 * ploidy);
    allele = malloc ( sizeof ( struct allele_t * ) * ploidy);
    edlib_alg = malloc ( sizeof ( EdlibAlignResult ) * ploidy );

    // First read of all the sequences
    for ( int i = 0; i < ploidy; i ++ ){
        // Read of the allele
        fm[i] = filemanager_init ( argv[optind] );
        if ( fm[i] == NULL ){
            exit (EXIT_FAILURE);
        }
        seq[i] = filemanager_next_seq ( fm[i], NULL );
        if ( seq[i] == NULL ){
            exit ( EXIT_FAILURE );
        }
        // Alignment of the allele
        char * align_string = NULL;
        if ( i != 0 ){
            int j = i + ploidy - 1;
            // Name of the alignment file
            fasta_fn = realloc ( fasta_fn, sizeof(char) * ( strlen(argv[optind]) + 5 ) );
            strcpy ( fasta_fn, argv[optind] );
            strcat ( fasta_fn, ".alg" );
            // Read of the alignment
            fm [j] = filemanager_init ( fasta_fn );
            if ( fm[j] == NULL ){
                exit ( EXIT_FAILURE );
            }
            seq [j] = filemanager_next_seq ( fm[j], NULL );
            if ( seq [j] == NULL ){
                exit ( EXIT_FAILURE );
            }
            align_string = seq[j]->sequence;
        }
        // Allele creation
        allele[i] = allele_point (
                seq[i]->sequence_size,
                seq[i]->sequence,
                align_string,
                NULL
                );
        optind ++;
    }
    curr_seq = seq[0];

    free ( fasta_fn );

    // Edlib configuration
    config = edlibNewAlignConfig ( -1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0 );
    stats[0] = stats_init ();
    stats[1] = stats_init ();

    // While there are sequences to read in the FASTA file
    while ( curr_seq != NULL ) {
        // Seek on the BAM
        if ( dictionary != NULL ){
            alias = tr_translate ( alias_index, curr_seq->label );
        }
        else{
            alias = curr_seq->label;
        }
        itr = bam_itr_querys( index, hdr, alias);
        if ( itr != NULL ){
            while( bam_itr_next( fp, itr, line ) > 0){

                // Read information
                pos = line->core.pos;
                len = line->core.l_qseq;
                uint8_t *read_seq = bam_get_seq( line ); // Read nucleotides
                uint8_t *qual = bam_get_qual ( line ); // Quality score

                // Interval of the reference
                flank_1 = floor ( log ( 2 * len ) / log ( 2 ) );
                flank_2 = flank_1;

                // Read string
                read = realloc ( read, sizeof ( char ) * ( len + 1 ) );
                int i;
                for ( i = 0; i < len; i++ ){
                    read[i] = seq_nt16_str[ bam_seqi ( read_seq, i ) ];
                }
                read[i] = 0;

                // Paired end
                if ( line->core.flag & BAM_FREAD2 ){
                    pair = 1;
                }
                else{
                    pair = 0;
                }

                for ( i = 0; i < ploidy; i ++ ){
                    curr_seq = seq[i];
                    // Seek on the allele
                    allele_seek ( pos, allele[i] );
                    
                    // Flanking regions
                    start = allele[i]->pos;
                    if ( start - flank_1 < 0 ){
                        start = 0;
                    }
                    else {
                        start -= flank_1;
                    }

                    end = allele[i]->pos + len;
                    if ( end + flank_2 >= curr_seq->sequence_size ){
                        end = curr_seq->sequence_size - 1;
                    }
                    else{
                        end += flank_2;
                    }

                    // Align
                    edlib_alg[i] = edlibAlign (
                            read,
                            len,
                            &curr_seq->sequence[start],
                            end - start,
                            config );
                    
                    // Select best alignment
                    if ( i == 0 || edlib_alg[i].editDistance < min_score ){
                        min_index = i;
                        min_score = edlib_alg[i].editDistance;
                    }

                    if ( verbose ){
                        printf ( "%d\n", edlib_alg[i].editDistance );
                        dump_read (
                            &curr_seq->sequence[start],
                            end - start,
                            edlib_alg[i].startLocations[0],
                            edlib_alg[i].alignment,
                            edlib_alg[i].alignmentLength,
                            read,
                            qual);
                    }
                }

                stats_update (
                        edlib_alg[min_index].alignment,
                        edlib_alg[min_index].alignmentLength,
                        stats[pair]
                        );

                for ( int i = 0; i < ploidy; i ++ ){
                    edlibFreeAlignResult ( edlib_alg[i] );
                }
            }
        }		
        // Next sequence
        for ( int i = 0; i < ploidy; i ++ ){
            seq[i] = filemanager_next_seq ( fm[i], seq[i] );
            char * align_string = NULL;
            if ( i != 0 ){
                int j = i + ploidy - 1;
                seq[j] = filemanager_next_seq ( fm[j], seq[j] );
                if ( seq[j] != NULL ){
                    align_string = seq[j]->sequence;
                }
            }
            if ( seq[i] != NULL) {
                // Update allele
                allele[i] = allele_point (
                        seq[i]->sequence_size,
                        seq[i]->sequence,
                        align_string,
                        allele[i]
                        );
            }
        }
        curr_seq = seq[0];
    }

    // Dump statistics
    if ( !silent ){
        stats_dump ( stdout, stats[0] );
        stats_dump ( stdout, stats[1] );
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i ++ ){
        filemanager_destroy ( fm[i] );
        free ( allele[i] );
        if ( i != 0 ){
            int j = i + ploidy - 1;
            filemanager_destroy ( fm [j] );
        }
    }
    free ( fm );
    free ( seq );
    free ( allele );
    free ( edlib_alg );
    bam_destroy1( line );
    bam_hdr_destroy ( hdr );
    bam_itr_destroy ( itr );
    hts_idx_destroy ( index );
    stats_destroy ( stats[0] );
    stats_destroy ( stats[1] );
    sam_close( fp );
    free ( read );
    tr_destroy ( alias_index );
    return 0;
}
