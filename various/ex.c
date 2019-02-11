

struct filemanager *fm;  /*  fasta/fastq file maneger  */
struct sequence_t *seq;

fm = filemanager_init (cfg->svalue);

seq = NULL;
while (seq == NULL) {
  seq = filemanager_next_seq (fm, seq);
 }

filemanager_destroy (fm);  /*  Releasefile managemant resources  */
