#ifndef QUANT_H
#define QUANT_H

#include "args.h"

#include <stdlib.h>
#include <complex.h>
#include <stdint.h>

struct quant_state_vector {
    uint64_t n;
    uint64_t len;
    double complex *elems;
};

void quant_state_vector_alloc(struct quant_state_vector *v, uint64_t n);

void quant_state_vector_init_random(struct quant_state_vector *v);

void quant_state_vector_print_info(struct quant_state_vector *v);

void quant_state_vector_print_elems(struct quant_state_vector *v);

void quant_state_vector_free(struct quant_state_vector *v);

typedef double complex quant_matrix[2][2];

void quant_matrix_init_hadamard(quant_matrix u);

int8_t quant_hadamard_transform(
    struct args *args,
    struct quant_state_vector *v,
    struct quant_state_vector *w,
    quant_matrix u);

#endif // QUANT_H