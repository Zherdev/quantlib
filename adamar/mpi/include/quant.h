#ifndef QUANT_H
#define QUANT_H

#include "args.h"

#include <stdlib.h>
#include <complex.h>
#include <stdint.h>
#include <mpi.h>

struct quant_state_vector {
    uint64_t len;
    double complex *elems;
};

void quant_state_vector_alloc(struct quant_state_vector *v, uint64_t qubits_num);

void quant_state_vector_init_random(struct quant_state_vector *v);

int8_t quant_state_vector_print_file(struct quant_state_vector *v, MPI_File fh);

int32_t quant_state_vector_snprint_info(
        char *buff, int32_t buff_size, struct quant_state_vector *v);

int32_t quant_state_vector_snprint_elems(
        char *buff, int32_t buff_size, struct quant_state_vector *v);

void quant_state_vector_free(struct quant_state_vector *v);

typedef double complex quant_matrix[2][2];

void quant_matrix_init_adamar(quant_matrix u);

int8_t quant_state_vector_print_file(struct quant_state_vector *v, MPI_File fh);

int8_t quant_state_vector_read_file(struct quant_state_vector *v, MPI_File fh);

int8_t quant_state_vector_cmp(struct quant_state_vector *a, struct quant_state_vector *b);

int8_t quant_adamar_transform(
    struct args *args,
    struct quant_state_vector *v,
    struct quant_state_vector *res,
    quant_matrix u);

#endif // QUANT_H