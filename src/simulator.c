/*
 * CNRSIM
 * variator.c
 * Generatos two alleles from a FASTA reference
 * and a VCF file containing known variations.
 *
 * @author Riccardo Massidda
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <time.h>
#include "model.h"
#include "stats.h"
#include "source.h"
#include "tandem.h"

// Init kseq structure
KSEQ_INIT ( gzFile, gzread );
// Macro to set a nucleotide
#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )
#define MAX_TARGETS 2048
#define MAX_TARGET_LENGTH 256

void usage ( char * name ) {
    fprintf ( stderr, "Usage: %s coverage error_model fastq [fastq ...]\n", name );
}

int main ( int argc, char ** argv ) {
    // Input
    int ploidy;
    char * model_name;
    char * fastq;
    // Error
    FILE * model_fp;
    model_t * model;
    read_t * generated = NULL;
    // FASTA
    gzFile * fp;
    kseq_t ** seq;
    int kseq_return;
    // BAM
    samFile * bam_fp = NULL;
    bam_hdr_t * bam_hdr = NULL;
    bam1_t * bam_entry = NULL;
    int io_check;
    int current_region;
    char qname[4096];
    int read_counter;
    // TODO: parameter via commandline
    bool is_bam = true;

    // Amplification
    char * amplified_seq = NULL;
    tandem_set_t * tandem = NULL;
    // Coverage
    int coverage;
    int sequenced;
    // Generation
    int insert_size = 0;
    int orientation = 0;
    int length;
    int reverse;
    int pos;
    stats_t * curr_end;

    // Init pseudorandom generator
    srand ( time ( NULL ) );

    // Non optional arguments
    if ( argc - optind < 3 ) {
        usage ( argv[0] );
        exit ( EXIT_FAILURE );
    }

    // Coverage
    coverage = atoi ( argv[optind++] );

    // Error Model
    model_name = argv[optind++];
    // File open
    model_fp = fopen ( model_name, "r" );
    if ( model_fp == NULL ) {
        fprintf ( stderr, "File %s not found.\n", model_name );
        exit ( EXIT_FAILURE );
    }

    // Init model
    model = model_parse ( model_fp );
    fclose ( model_fp );

    // Input sequences
    ploidy = argc - optind;
    fp = malloc ( sizeof ( gzFile ) * ploidy );
    seq = malloc ( sizeof ( kseq_t * ) * ploidy );

    // Open input files
    for ( int i = 0; i < ploidy; i ++ ){
        fastq = argv[optind++];
        fp[i] = gzopen ( fastq, "r" );
        if ( fp[i] == NULL ) {
            fprintf ( stderr, "File %s not found.\n", fastq );
            exit ( EXIT_FAILURE );
        }
        seq[i] = kseq_init ( fp[i] );
    }

    if ( is_bam ) {
      bam_fp = sam_open ( "/dev/stdout", "wb" );
      if ( bam_fp == NULL ) { 
        fprintf ( stderr, "File /dev/stdout not found.\n" );
        exit ( EXIT_FAILURE );
      }

      bam_hdr = bam_hdr_init ();
      bam_hdr->text = strdup ( "@HD\tVN:1.4\tSO:unknown\n" );
      bam_hdr->l_text = strlen ( bam_hdr->text );
      bam_hdr->n_targets = 0;
      bam_hdr->sdict = 0;
      bam_hdr->target_len = malloc ( sizeof ( uint32_t ) * MAX_TARGETS );
      bam_hdr->target_name = malloc ( sizeof ( char * ) * MAX_TARGETS );
      for ( int i = 0; i < MAX_TARGETS; i ++ ) {
        bam_hdr->target_name[i] = malloc ( sizeof ( char ) * MAX_TARGET_LENGTH );
      }

      for ( int i = 0; i < ploidy; i ++ ){
        kseq_return = kseq_read ( seq[i] );
        while ( kseq_return >= 0 || kseq_return == -2 ) {
          // If the region isn't already been added
          // can't use bam_name2id here, because it
          // generates the dictionary on the first
          // call, so it would be empty
          bool already_in = false;
          for ( int j = 0; j < bam_hdr->n_targets; j ++ ) {
            if ( strncmp ( seq[i]->name.s, bam_hdr->target_name[j], seq[i]->name.l ) == 0 ) {
              already_in = true;
              break;
            }
          }
          if ( ! already_in ) {
            bam_hdr->target_len[bam_hdr->n_targets] = seq[i]->name.l;
            memcpy ( bam_hdr->target_name[bam_hdr->n_targets], seq[i]->name.s, seq[i]->name.l + 1 );
            bam_hdr->n_targets ++;
          }
          kseq_return = kseq_read ( seq[i] );
        }
        kseq_destroy ( seq[i] );
        gzrewind ( fp[i] );
        seq[i] = kseq_init ( fp[i] );
      }

      io_check = sam_hdr_write ( bam_fp, bam_hdr );
      if ( io_check != 0 ) {
        fprintf ( stderr, "IO error\tsam_hdr_write\n" );
        exit ( EXIT_FAILURE );
      }

      bam_entry = bam_init1 ();
      if ( bam_entry == NULL ) {
        fprintf ( stderr, "Unable to create BAM record.\n" );
        exit ( EXIT_FAILURE );
      }
    }

    // Simulated read generation
    for ( int i = 0; i < ploidy; i ++ ){
        kseq_return = kseq_read ( seq[i] );
        while ( kseq_return >= 0 || kseq_return == -2 ) {
            if ( is_bam ) {
              current_region = bam_name2id ( bam_hdr, seq[i]->name.s );
              read_counter = 0;
            }
            // Analysis of the repetitions in the original sequence
            tandem = tandem_set_init ( seq[i]->seq.l, model->max_motif, model->max_repetition, tandem );
            tandem = tandem_set_analyze ( seq[i]->seq.s, seq[i]->seq.l, tandem );
            amplified_seq = realloc ( amplified_seq, sizeof ( char ) * ( seq[i]->seq.l * 2 ) );

            // Index for the FASTA sequence
            int seq_p = 0;
            // Index for the amplified sequence
            int aseq_p = 0;
            // Statistics
            int amp = 0;
            int deamp = 0;
            for ( int t = 0; t < tandem->n; t ++ ) {
              unsigned char in = tandem->set[t].rep;
              if ( in >= model->max_repetition ) {
                continue;
              }
              int gap = tandem->set[t].pos - seq_p;
              // Copy of the nucleotides between different tandems
              if ( gap > 0 ) {
                memcpy (
                    &amplified_seq[aseq_p],
                    &seq[i]->seq.s[seq_p],
                    sizeof ( char ) * gap );
                aseq_p += gap; 
                seq_p += gap;
              }
              else if ( gap < 0 ) {
                fprintf ( stderr, "The tandem set isn't ordered.\n" );
                exit ( EXIT_FAILURE );
              }
              unsigned char out = source_generate (
                  &in,
                  1,
                  tandem->set[t].pat,
                  model->amplification );
              ( in < out ) ? ++amp : ++deamp;
              int rep = out;
              memcpy (
                  &amplified_seq[aseq_p],
                  &seq[i]->seq.s[seq_p],
                  sizeof ( char ) * rep * tandem->set[t].pat );
              seq_p += ( tandem->set[t].rep * tandem->set[t].pat );
              aseq_p += ( rep * tandem->set[t].pat );
            }
            // Copy of the remaining sequence
            memcpy (
                &amplified_seq[aseq_p],
                &seq[i]->seq.s[seq_p],
                sizeof ( char ) * ( seq[i]->seq.l - seq_p ) );
            aseq_p += ( seq[i]->seq.l - seq_p );
            seq_p = seq[i]->seq.l;
            amplified_seq[aseq_p] = '\0';
            fprintf ( stderr, "%s\n", seq[i]->name.s );
            fprintf ( stderr, "\t(amplified):\t%d\t%d\t%.3f\n", aseq_p, seq_p, (aseq_p*100.0/seq_p));

            // Initial conditions
            sequenced = 0;
            pos = 0;
            curr_end = model->single;

            // Reach the coverage
            while ( coverage > ( sequenced / aseq_p ) ) {
              if ( curr_end == model->single ) {
                // Two bits: ++,+-,-+,--
                orientation = source_generate ( NULL, 0, 0, model->orientation );
                // Not sequenced nucleotides between pairs
                insert_size = source_generate ( NULL, 0, 0, model->insert_size );
                int lo_bound = insert_size * ( model->max_insert_size / model->size_granularity );
                int up_bound = ( insert_size + 1 ) * ( model->max_insert_size / model->size_granularity );
                insert_size = rand () % ( up_bound - lo_bound + 1 );
                insert_size += lo_bound;
                if ( is_bam ) {
                  read_counter ++;
                }
              }
              else {
                // There is no insert size after mate pair
                insert_size = 0;
              }

              // Generate new read
              generated = stats_generate_read ( &amplified_seq[pos], generated, curr_end );
              length = strlen ( generated->read );
              reverse = ( curr_end == model->single ) ? ( orientation & 2 ) : ( orientation & 1 );

              // The read reached the limit of the sequence
              if ( !generated->cut ) {
                if ( is_bam ) {
									uint8_t *s;
                  sprintf ( qname, "%sr%x", seq[i]->name.s, read_counter );
                  // bam_entry->data: qname-cigar-seq-qual-aux
                  // qname
                  bam_entry->core.l_qname = strlen ( qname ) + 1;
                  kroundup32 ( bam_entry->core.l_qname );
                  bam_entry->core.l_extranul = bam_entry->core.l_qname - strlen ( qname ) - 1;
                  bam_entry->l_data = bam_entry->core.l_qname;
                  // seq
                  bam_entry->l_data += length>>1;
                  // quality
                  bam_entry->l_data += length + 1;

                  // m_data === memory
                  // l_data === needed memory
                  if ( bam_entry->m_data < bam_entry->l_data ) {
                    bam_entry->m_data = bam_entry->l_data;
                    kroundup32 ( bam_entry->m_data );
                    bam_entry->data = realloc ( bam_entry->data, bam_entry->m_data );
                    if ( bam_entry->data == NULL ) {
                      fprintf ( stderr, "Cannot reallocate memory for a BAM entry.\n" );
                      exit ( EXIT_FAILURE );
                    }
                  }

									// QNAME
                  memcpy ( bam_entry->data, qname, bam_entry->core.l_qname );
                  memset ( &bam_entry->data[bam_entry->core.l_qname], 0, bam_entry->core.l_extranul );
                  
									// FLAG
									bam_entry->core.flag = BAM_FMUNMAP;
									// RNAME  tid
                  bam_entry->core.tid = current_region;
									// POS    pos
                  bam_entry->core.pos = pos;
									// MAPQ
                  bam_entry->core.qual = 255;

                  // TODO: convert alignment to CIGAR and include it
									// CIGAR
									bam_entry->core.n_cigar = 0;

                  // LSEQ
                  bam_entry->core.l_qseq = length;

                  // TODO: mate reads
									// RNEXT  mtid
                  bam_entry->core.mtid = current_region;
									// PNEXT  mpos
                  bam_entry->core.mpos = -1;

                  // SEQ
									s = bam_get_seq ( bam_entry );
									for ( int p = 0; p < bam_entry->core.l_qseq; p++ ){
										bam1_seq_seti( s, p, seq_nt16_table[(unsigned char)generated->read[p]] );
									}

                  // QUAL
									s = bam_get_qual ( bam_entry );
                  memcpy ( s, generated->quality, length + 1 );

                  io_check = sam_write1 ( bam_fp, bam_hdr, bam_entry );
                  if ( io_check < 0 ) {
                    fprintf ( stderr, "Unable to write record to BAM file\n" );
                    exit ( EXIT_FAILURE );
                  }
                }
                else {
                  // Reverse the order of the nucleotides
                  if ( reverse ) {
                    char swap;
                    for ( int j = 0; j < (length/2); j ++ ) {
                      swap = generated->read[j];
                      generated->read[j] = generated->read[length-j-1];
                      generated->read[length-j-1] = swap;
                      swap = generated->quality[j];
                      generated->quality[j] = generated->quality[length-j-1];
                      generated->quality[length-j-1] = swap;
                    }
                  }

                  // Adjust quality score for visualization
                  for ( int j = 0; j < length; j ++ ) {
                    generated->quality[j] += 33;
                  }

                  // Print result to file
                  printf ( "@%s %d %c\n", seq[i]->name.s, pos, ( reverse ) ? '-' : '+' );
                  printf ( "%s\n", generated->read );
                  printf ( "+\n" );
                  printf ( "%s\n\n", generated->quality );
                }

                // Update sequenced bases
                sequenced += length;

                // Change read-end
                curr_end = ( curr_end == model->single ) ? model->pair : model->single;
              }

              // Update start position
              pos += ( insert_size + length );
              if ( pos >= aseq_p ) {
                // Start from the beginning
                pos = 0;
                curr_end = model->single;
              }
            }
            fprintf ( stderr, "\t(sequenced):\t%d\t%d\t%.3f\n", sequenced, aseq_p, (sequenced/aseq_p)*100.0);

            // Read next sequence
            kseq_return = kseq_read ( seq[i] );
        }
    }

    // Cleanup
    for ( int i = 0; i < ploidy; i++ ) {
        gzclose ( fp[i] );
        kseq_destroy ( seq[i] );
    }
    if ( generated != NULL ){
        free ( generated->read );
        free ( generated->align );
        free ( generated->quality );
    }
    free ( generated );
    sam_close ( bam_fp );
    for ( int i = bam_hdr->n_targets; i < MAX_TARGETS; i ++ ) {
      free ( bam_hdr->target_name[i] );
    }
    bam_hdr_destroy ( bam_hdr );
    bam_destroy1 ( bam_entry );
    free ( fp );
    free ( amplified_seq );
    free ( seq );
    tandem_set_destroy ( tandem );
    model_destroy ( model );
    exit ( EXIT_SUCCESS );
}
