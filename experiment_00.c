#include "fileManager.h"

struct filemanager *fm;  /*  fasta/fastq file maneger  */
struct sequence_t *seq;

int main(void){
    char *filename = "data/chrY.fa";
    fm = filemanager_init (filename);
    if (fm == NULL){
        exit(EXIT_FAILURE);
    }
    seq = NULL;
    while (seq == NULL) {
      seq = filemanager_next_seq (fm, seq);
      printf("%s",seq->sequence);
    }

    filemanager_destroy (fm);  /*  Releasefile managemant resources  */
    return 0;
}
