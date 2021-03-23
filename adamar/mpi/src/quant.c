#include "quant.h"
#include "args.h"

#include <complex.h>
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

static uint64_t
quant_state_vector_get_part_len(uint64_t qubits_num, uint64_t my_rank, uint64_t comm_size)
{
    uint64_t full_len = pow(2, qubits_num);
    uint64_t part_len = full_len / comm_size; // No reminder here.
    return part_len;
}

void quant_state_vector_alloc(struct quant_state_vector *v, uint64_t qubits_num)
{
    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    v->len = quant_state_vector_get_part_len(qubits_num, my_rank, comm_size);
    v->elems = calloc(v->len, sizeof(double complex));
}

static uint64_t quant_state_vector_size(struct quant_state_vector *v)
{
    return v->len * sizeof(double complex);
}

void quant_state_vector_init_random(struct quant_state_vector *v)
{
    const int max_val = 5;
    for (uint64_t i = 0; i < v->len; i++) {
        double re = rand() % max_val;
        double im = rand() % max_val;

        v->elems[i] = re + im * I;
    }

    double sum = 0;
    for (uint64_t i = 0; i < v->len; i++) {
        double complex x = v->elems[i] * conj(v->elems[i]);
        sum += crealf(x);
    }

    double norm = sqrt(sum);
    for (uint64_t i = 0; i < v->len; i++) {
        v->elems[i] = v->elems[i] / norm;
    }
}

void static quant_state_vector_seek_file_by_steps(
        struct quant_state_vector *v, MPI_File fh, int my_rank, int step)
{
    uint64_t offset = v->len * my_rank;
    uint64_t cur = 0;
    uint64_t rem = offset % (uint64_t) step;
    while (cur < offset - rem) {
        MPI_File_seek(fh, step, MPI_SEEK_CUR);
        cur += step;
    }
    if (offset - cur > 0) {
        MPI_File_seek(fh, offset - cur, MPI_SEEK_CUR);
    }
}

void static quant_state_vector_read_file_by_steps(
        struct quant_state_vector *v, MPI_File fh, int my_rank, int step)
{
    uint64_t cur = 0;
    uint64_t rem = v->len % (uint64_t) step;
    while (cur < v->len - rem) {
        MPI_File_read(fh, v->elems + cur, step, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        cur += step;
    }
    if (v->len - cur > 0) {
        MPI_File_read(fh, v->elems + cur, v->len - cur, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    }
}

int8_t quant_state_vector_read_file(struct quant_state_vector *v, MPI_File fh)
{
    int step = (INT_MAX / sizeof(double complex)) / 2;
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    uint64_t to_read_bytes = quant_state_vector_size(v);
    if (to_read_bytes < (uint64_t) INT_MAX) {
        return MPI_File_read_ordered(fh, v->elems, v->len, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    }

    quant_state_vector_seek_file_by_steps(v, fh, my_rank, step);
    quant_state_vector_read_file_by_steps(v, fh, my_rank, step);
    return 0;
}

void static quant_state_vector_print_file_by_steps(
        struct quant_state_vector *v, MPI_File fh, int my_rank, int step)
{
    uint64_t cur = 0;
    uint64_t rem = v->len % (uint64_t) step;
    while (cur < v->len - rem) {
        MPI_File_write(fh, v->elems + cur, step, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        cur += step;
    }
    if (v->len - cur > 0) {
        MPI_File_write(fh, v->elems + cur, v->len - cur, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    }
}

int8_t quant_state_vector_print_file(struct quant_state_vector *v, MPI_File fh)
{
    int step = (INT_MAX / sizeof(double complex)) / 2;
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    uint64_t to_write_bytes = quant_state_vector_size(v);
    if (to_write_bytes < (uint64_t) INT_MAX) {
        return MPI_File_write_ordered(fh, v->elems, v->len, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
    }

    quant_state_vector_seek_file_by_steps(v, fh, my_rank, step);
    quant_state_vector_print_file_by_steps(v, fh, my_rank, step);
    return 0;
}

int32_t quant_state_vector_snprint_info(
        char *buff, int32_t buff_size, struct quant_state_vector *v)
{
    if (!v) {
        return 0;
    }

    int32_t offset = 0;

    offset = snprintf(buff, buff_size,
            "\tQuantum state vector len: %lu\n"
            "\tDouble complex size: %lu\n"
            "\tQuantum state vector size: %lu\n",
            v->len, sizeof(double complex), v->len * sizeof(double complex));

    return offset;
}

int32_t quant_state_vector_snprint_elems(
        char *buff, int32_t buff_size, struct quant_state_vector *v)
{
    if (!v) {
        return 0;
    }

    int32_t offset = 0;

    for (uint64_t i = 0; i < v->len; i++) {
        offset += snprintf(buff + offset, buff_size - offset,
                "\t%5lu: (%.4f, %.4f)\n", i, creal(v->elems[i]), cimag(v->elems[i]));
    }

    return offset;
}

void quant_state_vector_free(struct quant_state_vector *v)
{
    if (!v) {
        return;
    }

    free(v->elems);
    v->elems = NULL;
    v->len = 0;
}

void quant_matrix_init_adamar(quant_matrix u)
{
    const double v = 1 / sqrt(2);

    u[0][0] = v; u[0][1] = v;
    u[1][0] = v; u[1][1] = -v;
}

static void quant_adamar_transform_local(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        quant_matrix u)
{
    uint64_t len = v->len;
    uint64_t n = args->qubits_num;
    uint64_t k = args->target_qubit_num;

    for (uint64_t i = 0; i < len; i++) {
        uint8_t ik = (i >> (n - k)) & 1LLU;

        uint64_t ik_one_mask = 1LLU << (n - k);
        uint64_t ik_zero_mask = ~ik_one_mask;
        uint64_t i_zero = i & ik_zero_mask;
        uint64_t i_one = i_zero | ik_one_mask;

        res->elems[i] = u[ik][0] * v->elems[i_zero] + u[ik][1] * v->elems[i_one];
    }
}

static void quant_state_vector_sendrecv_by_parts(
        int my_rank, int other_rank, double complex *to_send,
        double complex *to_recv, uint64_t len, int step)
{
    uint64_t cur = 0;
    uint64_t rem = len % (uint64_t) step;
    while (cur < len - rem) {
        MPI_Sendrecv(to_send + cur, step, MPI_DOUBLE_COMPLEX, other_rank, 0,
                to_recv + cur, step, MPI_DOUBLE_COMPLEX, other_rank, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cur += step;
    }
    if (len - cur > 0) {
        MPI_Sendrecv(to_send + cur, len - cur, MPI_DOUBLE_COMPLEX, other_rank, 0,
                to_recv + cur, len - cur, MPI_DOUBLE_COMPLEX, other_rank, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

static void quant_adamar_transform_shared(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        quant_matrix u,
        uint64_t log_size,
        uint64_t my_rank)
{
    uint64_t len = v->len;
    uint64_t n = args->qubits_num;
    uint64_t k = args->target_qubit_num;
    uint64_t size = quant_state_vector_size(v);

    struct quant_state_vector tmp = {0};
    quant_state_vector_alloc(&tmp, n);

    uint8_t other_rank = my_rank ^ (1LLU << (log_size - k));

    if (size < (uint64_t) INT_MAX) {
        MPI_Sendrecv(v->elems, len, MPI_DOUBLE_COMPLEX, other_rank, 0,
                tmp.elems, len, MPI_DOUBLE_COMPLEX, other_rank, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        int step = (INT_MAX / sizeof(double complex)) / 2;
        quant_state_vector_sendrecv_by_parts(
                my_rank, other_rank, v->elems, tmp.elems, len, step);
    }


    if (my_rank < other_rank) {
        // ik = 0, ~ik = 1
        for (uint64_t i = 0; i < len; i++) {
            res->elems[i] = u[0][0] * v->elems[i] + u[0][1] * tmp.elems[i];
        }
    } else {
        // ik = 1, ~ik = 0
        for (uint64_t i = 0; i < len; i++) {
            res->elems[i] = u[1][0] * tmp.elems[i] + u[1][1] * v->elems[i];
        }
    }

    quant_state_vector_free(&tmp);
}

int8_t quant_state_vector_cmp(struct quant_state_vector *a, struct quant_state_vector *b)
{
    if (a->len != b->len) {
        return -1;
    }

    for (uint64_t i = 0; i < a->len; i++) {
        if (a->elems[i] != b->elems[i]) {
            return -1;
        }
    }

    return 0;
}

int8_t quant_adamar_transform(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        quant_matrix u)
{
    if (!v || !u || !res || !args) {
        return -1;
    }

    if (v->len != res->len) {
        return -1;
    }

    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    uint64_t log_size = log2(comm_size);

    if (args->target_qubit_num > log_size) {
        // In this case all data are accessible for the current process.
        quant_adamar_transform_local(args, v, res, u);
    } else {
        // In this case we need data from another process.
        quant_adamar_transform_shared(args, v, res, u, log_size, my_rank);
    }

    return 0;
}