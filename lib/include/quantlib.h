#ifndef QUANT_H
#define QUANT_H

#include <complex.h>
#include <stdint.h>
#include <mpi.h>

/*
 * quantlib -- parallel library for quantum calculations with MPI support.
 * Provides tools for quantum state manipulations
 * and N-qubit quantum transformations.
 *
 * Such common gates as Hadamard, NOT, CNOT, ROT and etc are provided.
 * Transformations with custom matrices are available also.
 *
 * MPI-processes number are expected to be 2^N, N >= 0.
 * N may be equal to zero, so it is possible to run calculations
 * in the non-parallel mode.
 *
 * See test suites in quantlib/test/src/test.c source codes
 * for quantlib usage examples.
 *
 * Zherdev, 2021.
 */

// quant state vector ==========================================================

// Represents quantum state of the system of qubits_num qubits.
struct quant_state_vector {
    uint64_t qubits_num;
    uint64_t len;
    double complex *elems;
};

// Allocate and init with zeros n-qubits state vector.
// Vector is distributed over MPI-processes equally.
// It is always required to allocate vector before any usage.
// quant_state_vector_free must be called later for memory deallocation.
//
// Returns not zero if error.
int8_t quant_state_vector_alloc(struct quant_state_vector *v, uint64_t qubits_num);

// Fill vector v with random values. v will be normalized.
//
// Returns not zero if error.
int8_t quant_state_vector_fill_random(struct quant_state_vector *v);

// Read the state vector v from the file fh.
// MPI-processes read file in the ordered way. Raw binary fowmat is used.
// Returns not zero if error.
int8_t quant_state_vector_read_file(struct quant_state_vector *v, MPI_File fh);

// Write the state vector v to the file fh.
// MPI-processes write file in the ordered way. Raw binary fowmat is used.
// Returns not zero if error.
int8_t quant_state_vector_print_file(struct quant_state_vector *v, MPI_File fh);

// Write formatted infromation text about vector to the buffer.
// Returns the length of bytes written.
int32_t quant_state_vector_snprint_info(
        char *buff, int32_t buff_size, struct quant_state_vector *v);

// Write formatted text with vector elements to the local buffer.
// Only the vector related to current MPI-process is printed.
// Returns the length of bytes written.
int32_t quant_state_vector_snprint_elems(
        char *buff, int32_t buff_size, struct quant_state_vector *v);

// Returns not zero if error.
int8_t quant_state_vector_free(struct quant_state_vector *v);

// If vectors are equal returns zero, else non zero.
int8_t quant_state_vector_cmp(struct quant_state_vector *a, struct quant_state_vector *b);

// quant state vector ^=========================================================

// quant matrix ================================================================

// Represents matrix of qubits_num qubits transformation.
struct quant_matrix {
    uint64_t qubits_num;
    uint64_t size;
    uint64_t len; // size^2
    double complex *elems;
};

// Allocate and init with zeros matrix of N-qubit transformation.
// Returns not zero if error.
int8_t quant_matrix_allocate(struct quant_matrix *u, uint64_t qubits_num);

// Returns not zero if error.
int8_t quant_matrix_free(struct quant_matrix *u);

// Returns not zero if error.
// 0 <= i, j < u->size.
int8_t quant_matrix_set(struct quant_matrix *u, uint64_t i, uint64_t j, double complex val);

// Expects correct params.
// 0 <= i, j < u->size.
double complex quant_matrix_get(struct quant_matrix *u, uint64_t i, uint64_t j);

// Allocate and init matrix of NOT gate.
// Returns not zero if error.
int8_t quant_matrix_init_not(struct quant_matrix *u);

// Allocate and init matrix of CNOT gate.
// Returns not zero if error.
int8_t quant_matrix_init_cnot(struct quant_matrix *u);

// Allocate and init matrix of ROT gate.
// Returns not zero if error.
int8_t quant_matrix_init_rot(struct quant_matrix *u);

// Allocate and init matrix of CROT gate.
// Returns not zero if error.
int8_t quant_matrix_init_crot(struct quant_matrix *u);

// quant matrix ^===============================================================

// quant transfromations =======================================================

// N-qubit transfromation with matrix u.
// Returns not zero if error.
int8_t quant_transform(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t *target_qubits,
        struct quant_matrix *u);

// One-qubit transfromation with matrix u.
// Returns not zero if error.
int8_t quant_transform_one(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_num,
        struct quant_matrix *u);

// Two-qubit transfromation with matrix u.
// Returns not zero if error.
int8_t quant_transform_two(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_1,
        uint64_t target_qubit_2,
        struct quant_matrix *u);

// One-qubit Hadamard transfromation.
// Returns not zero if error.
int8_t quant_transform_hadamard(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_num);

// Hadamard^N transfromation.
// Returns not zero if error.
int8_t quant_transform_hadamard_n(
        struct quant_state_vector *v,
        struct quant_state_vector *res);

// Represents NOT gate.
// Returns not zero if error.
int8_t quant_transform_not(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_num);

// Represents CNOT gate.
// Returns not zero if error.
int8_t quant_transform_cnot(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t control_qubit,
        uint64_t target_qubit);

// Represents ROT gate.
// Returns not zero if error.
int8_t quant_transform_rot(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t target_qubit_num);

// Represents CROT gate.
// Returns not zero if error.
int8_t quant_transform_crot(
        struct quant_state_vector *v,
        struct quant_state_vector *res,
        uint64_t control_qubit,
        uint64_t target_qubit);

// quant transfromations ^======================================================

#endif // QUANT_H
