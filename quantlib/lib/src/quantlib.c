#include "quantlib.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <inttypes.h>

/*
 * See quantlib.h for details.
 * Zherdev, 2021.
 */

// quant util ==================================================================

/*
 * quant_util_sendrecv* is used to avoid internal MPI bugs
 * with int length buffer overflow.
 */

static void quant_util_sendrecv_by_parts(
        int other_rank, double complex *to_send,
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

static void quant_util_sendrecv(
        int other_rank, double complex *to_send,
        double complex *to_recv, uint64_t len, uint64_t size)
{
    if (size < (uint64_t) INT_MAX) {
        MPI_Sendrecv(
                to_send, len, MPI_DOUBLE_COMPLEX, other_rank, 0,
                to_recv, len, MPI_DOUBLE_COMPLEX, other_rank, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        int step = (INT_MAX / sizeof(double complex)) / 2;
        quant_util_sendrecv_by_parts(other_rank, to_send, to_recv, len, step);
    }
}

// quant util ^=================================================================

// quant state vector ==========================================================

static uint64_t
quant_state_vector_get_part_len(uint64_t qubits_num, uint64_t my_rank, uint64_t comm_size)
{
    uint64_t full_len = pow(2, qubits_num);
    uint64_t part_len = full_len / comm_size; // No reminder here.
    return part_len;
}

int8_t quant_state_vector_alloc(struct quant_state_vector *v, uint64_t qubits_num)
{
    if (!v || !qubits_num) {
        return -1;
    }

    int comm_size = 0;
    int8_t err = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    if (err != MPI_SUCCESS) {
        return -1;
    }

    int my_rank = 0;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (err != MPI_SUCCESS) {
        return -1;
    }

    v->qubits_num = qubits_num;
    v->len = quant_state_vector_get_part_len(qubits_num, my_rank, comm_size);
    v->elems = calloc(v->len, sizeof(double complex));
    if (!v->elems) {
        return -1;
    }

    return 0;
}

static uint64_t quant_state_vector_size(struct quant_state_vector *v)
{
    return v->len * sizeof(double complex);
}

int8_t quant_state_vector_fill_random(struct quant_state_vector *v)
{
    if (!v || !v->elems) {
        return -1;
    }

    const int max_val = 5;
    for (uint64_t i = 0; i < v->len; i++) {
        double re = rand() % max_val;
        double im = rand() % max_val;

        v->elems[i] = re + im * I;
    }

    double local_sum = 0;
    for (uint64_t i = 0; i < v->len; i++) {
        double complex x = v->elems[i] * conj(v->elems[i]);
        local_sum += creal(x);
    }

    double sum = 0;
    int8_t err = MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        return -1;
    }

    double norm = sqrt(sum);
    for (uint64_t i = 0; i < v->len; i++) {
        v->elems[i] = v->elems[i] / norm;
    }

    return 0;
}

int8_t static quant_state_vector_seek_file_by_steps(
        struct quant_state_vector *v, MPI_File fh, int my_rank, int step)
{
    uint64_t offset = v->len * my_rank;
    uint64_t cur = 0;
    uint64_t rem = offset % (uint64_t) step;

    while (cur < offset - rem) {
        int8_t err = MPI_File_seek(fh, step, MPI_SEEK_CUR);
        if (err != MPI_SUCCESS) {
            return -1;
        }
        cur += step;
    }

    if (offset - cur > 0) {
        int8_t err = MPI_File_seek(fh, offset - cur, MPI_SEEK_CUR);
        if (err != MPI_SUCCESS) {
            return -1;
        }
    }

    return 0;
}

int8_t static quant_state_vector_read_file_by_steps(
        struct quant_state_vector *v, MPI_File fh, int my_rank, int step)
{
    uint64_t cur = 0;
    uint64_t rem = v->len % (uint64_t) step;

    while (cur < v->len - rem) {
        int8_t err = MPI_File_read(
                fh, v->elems + cur, step, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        if (err != MPI_SUCCESS) {
            return -1;
        }
        cur += step;
    }

    if (v->len - cur > 0) {
        int8_t err = MPI_File_read(
                fh, v->elems + cur, v->len - cur, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        if (err != MPI_SUCCESS) {
            return -1;
        }
    }

    return 0;
}

int8_t quant_state_vector_read_file(struct quant_state_vector *v, MPI_File fh)
{
    if (!v) {
        return -1;
    }

    int my_rank = 0;
    int8_t err = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (err != MPI_SUCCESS) {
        return -1;
    }

    uint64_t to_read_bytes = quant_state_vector_size(v);
    if (to_read_bytes < (uint64_t) INT_MAX) {
        err = MPI_File_read_ordered(
                fh, v->elems, v->len, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        if (err != MPI_SUCCESS) {
            return -1;
        }
    }

    int step = (INT_MAX / sizeof(double complex)) / 2;
    err = quant_state_vector_seek_file_by_steps(v, fh, my_rank, step);
    if (err != MPI_SUCCESS) {
        return -1;
    }
    err = quant_state_vector_read_file_by_steps(v, fh, my_rank, step);
    if (err != MPI_SUCCESS) {
        return -1;
    }

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
            "\tQuantum state vector len: %" PRId64 "\n"
            "\tDouble complex size: %lu\n"
            "\tQuantum state vector size: %" PRId64 "\n",
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
                "\t%" PRId64 ": (%.4f, %.4f)\n", i, creal(v->elems[i]), cimag(v->elems[i]));
    }

    return offset;
}

int8_t quant_state_vector_free(struct quant_state_vector *v)
{
    if (!v) {
        return -1;
    }

    free(v->elems);
    v->elems = NULL;
    v->len = 0;
    v->qubits_num = 0;

    return 0;
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

// quant state vector ^=========================================================

// quant matrix ================================================================

int8_t quant_matrix_allocate(struct quant_matrix *u, uint64_t qubits_num)
{
    if (!u || !qubits_num) {
        return -1;
    }

    u->qubits_num = qubits_num;
    u->size = (uint64_t) round(pow(2, qubits_num));
    u->len = u->size * u->size;
    u->elems = calloc(u->len, sizeof(double complex));
    if (!u->elems) {
        return -1;
    }

    return 0;
}

int8_t quant_matrix_free(struct quant_matrix *u)
{
    if (!u) {
        return -1;
    }

    free(u->elems);
    u->elems = NULL;
    u->len = 0;
    u->qubits_num = 0;
    u->size = 0;

    return 0;
}

int8_t quant_matrix_set(struct quant_matrix *u, uint64_t i, uint64_t j, double complex val)
{
    if (!u || !u->elems) {
        return -1;
    }
    if (u->size <= i || u->size <= j) {
        return -1;
    }

    u->elems[i * u->size + j] = val;
    return 0;
}

double complex quant_matrix_get(struct quant_matrix *u, uint64_t i, uint64_t j)
{
    return u->elems[i * u->size + j];
}

int8_t quant_matrix_init_hadamard(struct quant_matrix *u)
{
    int8_t err = quant_matrix_allocate(u, 1);
    if (err) {
        return err;
    }

    const double v = 1 / sqrt(2);
    err = quant_matrix_set(u, 0, 0, v) || quant_matrix_set(u, 0, 1, v) ||
          quant_matrix_set(u, 1, 0, v) || quant_matrix_set(u, 1, 1, -v);
    return err;
}

int8_t quant_matrix_init_not(struct quant_matrix *u)
{
    int8_t err = quant_matrix_allocate(u, 1);
    if (err) {
        return err;
    }

    err = quant_matrix_set(u, 0, 0, 0) || quant_matrix_set(u, 0, 1, 1) ||
          quant_matrix_set(u, 1, 0, 1) || quant_matrix_set(u, 1, 1, 0);
    return err;
}

int8_t quant_matrix_init_cnot(struct quant_matrix *u)
{
    int8_t err = quant_matrix_allocate(u, 2);
    if (err) {
        return err;
    }

    err = quant_matrix_set(u, 0, 0, 1) || quant_matrix_set(u, 0, 1, 0) || quant_matrix_set(u, 0, 2, 0) || quant_matrix_set(u, 0, 3, 0) ||
          quant_matrix_set(u, 1, 0, 0) || quant_matrix_set(u, 1, 1, 1) || quant_matrix_set(u, 1, 2, 0) || quant_matrix_set(u, 1, 3, 0) ||
          quant_matrix_set(u, 2, 0, 0) || quant_matrix_set(u, 2, 1, 0) || quant_matrix_set(u, 2, 2, 0) || quant_matrix_set(u, 2, 3, 1) ||
          quant_matrix_set(u, 3, 0, 0) || quant_matrix_set(u, 3, 1, 0) || quant_matrix_set(u, 3, 2, 1) || quant_matrix_set(u, 3, 3, 0);
    return err;
}

int8_t quant_matrix_init_rot(struct quant_matrix *u)
{
    int8_t err = quant_matrix_allocate(u, 1);
    if (err) {
        return err;
    }

    err = quant_matrix_set(u, 0, 0, 1) || quant_matrix_set(u, 0, 1, 0) ||
          quant_matrix_set(u, 1, 0, 0) || quant_matrix_set(u, 1, 1, -1);
    return err;
}

int8_t quant_matrix_init_crot(struct quant_matrix *u)
{
    int8_t err = quant_matrix_allocate(u, 2);
    if (err) {
        return err;
    }

    err = quant_matrix_set(u, 0, 0, 1) || quant_matrix_set(u, 0, 1, 0) || quant_matrix_set(u, 0, 2, 0) || quant_matrix_set(u, 0, 3, 0) ||
          quant_matrix_set(u, 1, 0, 0) || quant_matrix_set(u, 1, 1, 1) || quant_matrix_set(u, 1, 2, 0) || quant_matrix_set(u, 1, 3, 0) ||
          quant_matrix_set(u, 2, 0, 0) || quant_matrix_set(u, 2, 1, 0) || quant_matrix_set(u, 2, 2, 1) || quant_matrix_set(u, 2, 3, 0) ||
          quant_matrix_set(u, 3, 0, 0) || quant_matrix_set(u, 3, 1, 0) || quant_matrix_set(u, 3, 2, 0) || quant_matrix_set(u, 3, 3, -1);
    return err;
}

// quant matrix ^===============================================================

// quant transfromations =======================================================

static uint64_t quant_transform_masks_num(uint64_t target_qubits_num)
{
    return 1llu << target_qubits_num;
}

static void quant_transform_masks_init(
        uint64_t *masks,         uint64_t  masks_num,
        uint64_t *target_qubits, uint64_t  target_qubits_num,
        uint64_t  qubits_num)
{
    for (uint64_t m = 0; m < masks_num; m++) {
		masks[m] = 0;

		for (uint64_t q = 0; q < target_qubits_num; q++) {
            if (!((m >> q) & 1)) {
                continue;
            }

            uint64_t offset = qubits_num - target_qubits[target_qubits_num - 1 - q];
            masks[m] ^= 1llu << offset;
        }
	}
}

static void quant_transform_ranks_init(
        uint64_t *ranks,         uint64_t *ranks_num,
        uint64_t *target_qubits, uint64_t  target_qubits_num,
        uint64_t log_size,       uint64_t  my_rank)
{
    	uint64_t shared_qubits[target_qubits_num];
        uint64_t j = 0;

        for (uint64_t q = 0; q < target_qubits_num; q++) {
            if (target_qubits[q] <= log_size) {
                shared_qubits[j] = target_qubits[q];
                j++;
            }
        }
        uint64_t shared_qubits_num = j;

        *ranks_num = quant_transform_masks_num(shared_qubits_num);
        quant_transform_masks_init(
            ranks, *ranks_num,
            shared_qubits, shared_qubits_num,
            log_size);

        ranks[0] = 1 << shared_qubits_num;
        for (uint64_t r = 1; r < ranks[0]; r++) {
            ranks[r] ^= my_rank;
        }
        --ranks[0];
}

int8_t quant_transform(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t *target_qubits,
        struct quant_matrix *u)
{
    if (!v || !u || !res || !target_qubits) {
        return -1;
    }
    if (!v->elems || !u->elems || !res->elems) {
        return -1;
    }
    if (v->len != res->len || v->qubits_num != res->qubits_num) {
        return -1;
    }

    int comm_size = 0;
    int8_t err = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    if (err != MPI_SUCCESS) {
        return -1;
    }

    int my_rank = 0;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (err != MPI_SUCCESS) {
        return -1;
    }

    uint64_t log_size = log2(comm_size);
    uint64_t vec_size = quant_state_vector_size(v);
    uint64_t qubits_num = v->qubits_num;
    uint64_t vec_len = v->len;

    // Number of the qubits that will be transformed.
    uint64_t target_qubits_num = u->qubits_num;

    uint64_t masks_num = quant_transform_masks_num(target_qubits_num);
	uint64_t masks[masks_num];
    quant_transform_masks_init(
            masks, masks_num,
            target_qubits, target_qubits_num,
            qubits_num);

    // Get MPI processes ranks for communications.
    uint64_t ranks_num = 0;
    uint64_t ranks[masks_num]; // ranks_num <= masks_num, so masks_num is always enough.
    quant_transform_ranks_init(
            ranks, &ranks_num,
            target_qubits, target_qubits_num,
            log_size, my_rank);

    int8_t need_share = ranks[0] > 0; // need sendrecv vector parts.

    // Local buffers for vector parts from other processes.
    complex double **buffs;
    if (need_share) {
        buffs = calloc(ranks[0], sizeof(complex double *));
        for (uint64_t r = 0; r < ranks[0]; r++) {
            buffs[r] = calloc(1, vec_size);
        }
    }

    uint64_t rank2id[comm_size];
	for (uint64_t r = 0; r < ranks[0]; r++) {
        int other_rank = ranks[r+1];
        double complex *to_send = v->elems;
        complex double *to_recv = buffs[r];

        quant_util_sendrecv(other_rank, to_send, to_recv, vec_len, vec_size);

		rank2id[ranks[r+1]] = r;
	}

    uint64_t l = 1 << target_qubits_num;
    double complex *vec_part = NULL;
    vec_part = calloc(l, sizeof(double complex));

    uint64_t x = 0;   // matrix line number.
    uint64_t id0 = 0; // rank id.
    uint64_t num = 0; // element number.
    uint64_t id1 = 0; // local id.
    uint64_t tmp = 0;
    const uint64_t tmp0 = (uint64_t) my_rank << (qubits_num - log_size);
    const uint64_t tmp1 = ~((~0llu) << (qubits_num - log_size));

    for (uint64_t j = 0; j < vec_len; j++) {
        tmp = (tmp0 ^ j) & ~masks[l-1];

        for (uint64_t z = 0; z < l; z++) {
            num = tmp ^ masks[z];
            id0 = num >> (qubits_num - log_size);
            id1 = num & tmp1;

            if (id0 == (uint64_t) my_rank) {
                vec_part[z] = v->elems[id1];
                if (id1 == j) {
                    x = z;
                }
            } else {
                uint64_t other_id = rank2id[id0];
                vec_part[z] = buffs[other_id][id1];
            }
        }

        res->elems[j] = 0;
        for (uint64_t z = 0; z < l; z++) {
            res->elems[j] += quant_matrix_get(u, x, z) * vec_part[z];
        }
    }

	if (need_share) {
        for (uint64_t r = 0; r < ranks[0]; r++) {
            free(buffs[r]);
        }
        free(buffs);
    }
    free(vec_part);

    return 0;
}

int8_t quant_transform_one(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_num,
        struct quant_matrix *u)
{
    if (target_qubit_num < 1 || target_qubit_num > v->qubits_num) {
        return -1;
    }
    if (u->qubits_num != 1) {
        return -1;
    }

    uint64_t target_qubits[1] = {target_qubit_num};

    int8_t err = quant_transform(v, res, target_qubits, u);
    return err;
}

int8_t quant_transform_two(
    struct quant_state_vector *v,
    struct quant_state_vector *res,
    uint64_t target_qubit_1,
    uint64_t target_qubit_2,
    struct quant_matrix *u)
{
    if (target_qubit_1 < 1 || target_qubit_1 > v->qubits_num ||
            target_qubit_2 < 1 || target_qubit_2 > v->qubits_num) {
        return -1;
    }
    if (u->qubits_num != 2) {
        return -1;
    }

    uint64_t target_qubits[2] = {target_qubit_2, target_qubit_1};

    int8_t err = quant_transform(v, res, target_qubits, u);
    return err;

    return 0;
}

int8_t quant_transform_hadamard(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_num)
{
    struct quant_matrix u = {0};
    int8_t err = quant_matrix_init_hadamard(&u);
    if (err) {
        return err;
    }

    err = quant_transform_one(v, res, target_qubit_num, &u);
    err |= quant_matrix_free(&u);
    return err;
}

int8_t quant_transform_hadamard_n(
        struct quant_state_vector *v,
        struct quant_state_vector *res)
{
    if (!v || !res) {
        return -1;
    }

    if (v->len != res->len || v->qubits_num != res->qubits_num) {
        return -1;
    }

    struct quant_matrix u = {0};
    int8_t err = quant_matrix_init_hadamard(&u);
    if (err) {
        return err;
    }

    // Temporary vectors for swapping during
    // the cycle of quant hadamard transformations.
    struct quant_state_vector next = *res;
    struct quant_state_vector tmp = {0};
    struct quant_state_vector prev = {0};

    err = quant_state_vector_alloc(&prev, v->qubits_num);
    if (err) {
        quant_matrix_free(&u);
        return err;
    }
    memcpy(prev.elems, v->elems, sizeof(double complex) * v->len);

    for (uint64_t i = 1; i <= v->qubits_num; i++) {
        err |= quant_transform_one(&prev, &next, i, &u);

        tmp = prev;
        prev = next;
        next = tmp;
    }
    *res = prev;

    quant_state_vector_free(&tmp);
    quant_matrix_free(&u);

    return err;
}

int8_t quant_transform_not(
    struct quant_state_vector *v,
    struct quant_state_vector *res,
    uint64_t target_qubit_num)
{
    struct quant_matrix u = {0};
    int8_t err = quant_matrix_init_not(&u);
    if (err) {
        return err;
    }

    err = quant_transform_one(v, res, target_qubit_num, &u);
    err |= quant_matrix_free(&u);
    return err;
}

int8_t quant_transform_cnot(
    struct quant_state_vector *v,
    struct quant_state_vector *res,
    uint64_t control_qubit,
    uint64_t target_qubit)
{
    struct quant_matrix u = {0};
    int8_t err = quant_matrix_init_cnot(&u);
    if (err) {
        return err;
    }

    err = quant_transform_two(v, res, control_qubit, target_qubit, &u);
    err |= quant_matrix_free(&u);
    return err;
}

int8_t quant_transform_rot(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_num)
{
    struct quant_matrix u = {0};
    int8_t err = quant_matrix_init_rot(&u);
    if (err) {
        return err;
    }

    err = quant_transform_one(v, res, target_qubit_num, &u);
    err |= quant_matrix_free(&u);
    return err;
}

int8_t quant_transform_crot(
    struct quant_state_vector *v,
    struct quant_state_vector *res,
    uint64_t control_qubit,
    uint64_t target_qubit)
{
    struct quant_matrix u = {0};
    int8_t err = quant_matrix_init_crot(&u);
    if (err) {
        return err;
    }

    err = quant_transform_two(v, res, control_qubit, target_qubit, &u);
    err |= quant_matrix_free(&u);
    return err;
}

// quant transfromations ^======================================================
