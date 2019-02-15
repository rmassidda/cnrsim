#ifndef VARIATOR_H
#define VARIATOR_H

typedef struct allele_t allele_t;

struct allele_t {
    char *sequence;
    long int buffer_size;
    long int pos;
    long int off;
};

#endif
