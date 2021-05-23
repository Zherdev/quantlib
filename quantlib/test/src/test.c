#include "quantlib.h"
#include "test.h"

#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>

struct test_runner {
    int32_t my_rank;
    int32_t world_size;

    struct quant_state_vector v1;
    struct quant_state_vector v2;
    struct quant_state_vector res1;
    struct quant_state_vector res2;
};

static void test_print(struct test_runner *tester, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    const int root_rank = 0;

    if (tester->my_rank == root_rank) {
        vprintf(fmt, args);
    }

    va_end(args);
}

static void test_split_vector(
        struct test_runner *tester,
        double complex *elems,
        struct quant_state_vector *v)
{
    uint64_t offset = v->len * tester->my_rank;
    for (uint64_t i = 0; i < v->len; i++) {
        v->elems[i] = elems[offset + i];
    }
}

// test suite ==================================================================

typedef int8_t (*test_func) (struct test_runner *tester);

struct test_suite {
    const char *name;
    test_func   test;
    uint32_t    qubits_num;
};

static void test_suite_before(struct test_runner *tester, struct test_suite *suite)
{
    quant_state_vector_alloc(&tester->v1, suite->qubits_num);
    quant_state_vector_alloc(&tester->v2, suite->qubits_num);
    quant_state_vector_alloc(&tester->res1, suite->qubits_num);
    quant_state_vector_alloc(&tester->res2, suite->qubits_num);

    quant_state_vector_fill_random(&tester->v1);
    quant_state_vector_fill_random(&tester->v2);
    quant_state_vector_fill_random(&tester->res1);
    quant_state_vector_fill_random(&tester->res2);
}

static int8_t test_suite_run(struct test_runner *tester, struct test_suite *suite)
{
    if (!suite) {
        return -1;
    }

    test_print(tester, "%s... ", suite->name);
    int8_t err = suite->test(tester);
    test_print(tester, err ? "\tERROR!\n" : "\tOK!\n");

    return err;
}

static void test_suite_after(struct test_runner *tester)
{
    quant_state_vector_free(&tester->v1);
    quant_state_vector_free(&tester->v2);
    quant_state_vector_free(&tester->res1);
    quant_state_vector_free(&tester->res2);
}

static int32_t test_run_suites(
        struct test_runner *tester, struct test_suite *suites, int32_t suites_num) {
    if (!suites) {
        return -1;
    }

    int32_t errs = 0;
    for (int32_t i = 0; i < suites_num; i++) {
        MPI_Barrier(MPI_COMM_WORLD);

        test_suite_before(tester, suites);

        int8_t err = test_suite_run(tester, suites);
        if (err) {
            errs++;
        }

        test_suite_after(tester);

        suites++;
    }

    return errs;
}

// test suite ^=================================================================

// test asserts ================================================================

static int8_t test_assert(
        struct test_runner *tester, double complex actual, double complex expected, const char *msg)
{
    const double precision = 1e-10;
    int8_t assertion = cabs(actual - expected) < precision;

    if (!assertion) {
        test_print(tester,
                "test_assert failed: got %lf+i%lf, expected %lf+i%lf, %s\n",
                creal(actual), cimag(actual), creal(expected), cimag(expected), msg);
    }

    return !assertion;
}

static int8_t test_assert_vector(
        struct test_runner *tester,
        struct quant_state_vector *actual,
        struct quant_state_vector *expected,
        const char *msg)
{
    if (expected->len != actual->len || expected->qubits_num != actual->qubits_num) {
        test_print(tester, "test_assert_vector failed: different length, %s\n", msg);
        return 1;
    }

    for (uint64_t i = 0; i < expected->len; i++) {
        int8_t err = test_assert(tester, actual->elems[i], expected->elems[i], msg);
        if (err) {
            return err;
        }
    }
    return 0;
}

static int8_t test_assert_is_normalized(struct test_runner *tester, struct quant_state_vector *v)
{
    double local_sum = 0;
    for (uint64_t i = 0; i < v->len; i++) {
        double complex x = v->elems[i] * conj(v->elems[i]);
        local_sum += creal(x);
    }

    double sum = 0;
    MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return test_assert(tester, sum, 1.0, "vector is not normalized!");
}

// Scalar multiply v1 by v2 vector.
static double complex test_scalar_product(
        struct quant_state_vector *v1,
        struct quant_state_vector *v2)
{
    double complex local_sum = 0;
    for (uint64_t i = 0; i < v1->len; i++) {
        local_sum += v1->elems[i] * conj(v2->elems[i]);
    }

    double complex prod = 0;
    MPI_Allreduce(&local_sum, &prod, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    return prod;
}

// Check that scalar product gives the same result for v and res vectors.
static int8_t test_assert_scalar_prod_equals(struct test_runner *tester)
{
    double complex prod_v = test_scalar_product(&tester->v1, &tester->v2);
    double complex prod_res = test_scalar_product(&tester->res1, &tester->res2);

    int8_t err = test_assert(tester, prod_v, prod_res,
            "scalar product is not preserved (real part)!");
    return err;
}

// test asserts ^===============================================================

// test case funcs =============================================================

// Checks that the unitary transform preserves the norm and the scalar product.
static int8_t test_blackbox_check(struct test_runner *tester)
{
    int8_t err = test_assert_is_normalized(tester, &tester->v1);
    if (err) {
        return err;
    }
    err = test_assert_is_normalized(tester, &tester->v2);
    if (err) {
        return err;
    }
    err = test_assert_is_normalized(tester, &tester->res1);
    if (err) {
        return err;
    }
    err = test_assert_is_normalized(tester, &tester->res2);
    if (err) {
        return err;
    }

    err = test_assert_scalar_prod_equals(tester);
    if (err) {
        return err;
    }

    return 0;
}

static int8_t test_hadamard_blackbox(struct test_runner *tester)
{
    uint64_t target_qubit_num = 2;

    int8_t err = 0;
    err = quant_transform_hadamard(&tester->v1, &tester->res1, target_qubit_num);
    err |= quant_transform_hadamard(&tester->v2, &tester->res2, target_qubit_num);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_blackbox_check(tester);
    return err;
}

static int8_t test_hadamard_whitebox(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            1.0 / 3.0,
            -2.0 * I / 3.0,
            2.0 / 3.0,
            0.0
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            1.0 / sqrt(2.0),
            -2.0 * I / (3 * sqrt(2.0)),
            -1.0 / (3 * sqrt(2.0)),
            -2.0 * I / (3 * sqrt(2.0))
    };
    test_split_vector(tester, expected, expected_vec);

    uint64_t target_qubit_num = 1;

    int8_t err = 0;
    err = quant_transform_hadamard(src_vec, res_vec, target_qubit_num);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

static int8_t test_hadamard_n_blackbox(struct test_runner *tester)
{
    int8_t err = 0;
    err = quant_transform_hadamard_n(&tester->v1, &tester->res1);
    err |= quant_transform_hadamard_n(&tester->v2, &tester->res2);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_blackbox_check(tester);
    return err;
}

static int8_t test_hadamard_n_whitebox(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            1.0 / 3.0,
            -2.0 * I / 3.0,
            2.0 / 3.0,
            0.0
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            0.5 - I / 3.0,
            0.5 + I / 3.0,
            -1.0 / 6.0 - I / 3.0,
            -1.0 / 6.0 + I / 3.0
    };
    test_split_vector(tester, expected, expected_vec);

    int8_t err = quant_transform_hadamard_n(src_vec, res_vec);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

static int8_t test_not_blackbox(struct test_runner *tester)
{
    uint64_t target_qubit_num = 2;

    int8_t err = 0;
    err = quant_transform_not(&tester->v1, &tester->res1, target_qubit_num);
    err |= quant_transform_not(&tester->v2, &tester->res2, target_qubit_num);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_blackbox_check(tester);
    return err;
}

static int8_t test_not_whitebox(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            1, 0, 0, 0
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            0, 0, 1, 0
    };
    test_split_vector(tester, expected, expected_vec);

    uint64_t target_qubit_num = 1;

    int8_t err = 0;
    err = quant_transform_not(src_vec, res_vec, target_qubit_num);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

static int8_t test_cnot_blackbox(struct test_runner *tester)
{
    uint64_t control_qubit = 1;
    uint64_t target_qubit = 4;

    int8_t err = 0;
    err = quant_transform_cnot(&tester->v1, &tester->res1, control_qubit, target_qubit);
    err |= quant_transform_cnot(&tester->v2, &tester->res2, control_qubit, target_qubit);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_blackbox_check(tester);
    return err;
}

static int8_t test_cnot_whitebox_no_invert(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            0, 0, 1, 0, 0, 0, 0, 0
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            0, 0, 1, 0, 0, 0, 0, 0
    };
    test_split_vector(tester, expected, expected_vec);

    uint64_t control_qubit = 1;
    uint64_t target_qubit = 3;

    int8_t err = 0;
    err = quant_transform_cnot(src_vec, res_vec, control_qubit, target_qubit);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

static int8_t test_cnot_whitebox_invert(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            0, 1, 0, 0, 0, 0, 0, 0
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            0, 0, 0, 0, 0, 1, 0, 0
    };
    test_split_vector(tester, expected, expected_vec);

    uint64_t control_qubit = 1;
    uint64_t target_qubit = 3;

    int8_t err = 0;
    err = quant_transform_cnot(src_vec, res_vec, control_qubit, target_qubit);

    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

static int8_t test_rot_blackbox(struct test_runner *tester)
{
    uint64_t target_qubit_num = 5;

    int8_t err = 0;
    err = quant_transform_rot(&tester->v1, &tester->res1, target_qubit_num);
    err |= quant_transform_rot(&tester->v2, &tester->res2, target_qubit_num);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_blackbox_check(tester);
    return err;
}

static int8_t test_rot_whitebox(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            0, 0, 1, 0, 0, 0, 0, 0
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            0, 0, -1, 0, 0, 0, 0, 0
    };
    test_split_vector(tester, expected, expected_vec);

    uint64_t target_qubit_num = 2;

    int8_t err = 0;
    err = quant_transform_rot(src_vec, res_vec, target_qubit_num);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

static int8_t test_crot_blackbox(struct test_runner *tester)
{
    uint64_t control_qubit = 2;
    uint64_t target_qubit = 2;

    int8_t err = 0;
    err = quant_transform_crot(&tester->v1, &tester->res1, control_qubit, target_qubit);
    err |= quant_transform_crot(&tester->v2, &tester->res2, control_qubit, target_qubit);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_blackbox_check(tester);
    return err;
}

static int8_t test_crot_whitebox_no_rotation(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            0, 0, 0, 0, 1, 0, 0, 0
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            0, 0, 0, 0, 1, 0, 0, 0
    };
    test_split_vector(tester, expected, expected_vec);

    uint64_t control_qubit = 1;
    uint64_t target_qubit = 2;

    int8_t err = 0;
    err = quant_transform_crot(src_vec, res_vec, control_qubit, target_qubit);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

static int8_t test_crot_whitebox_rotation(struct test_runner *tester)
{
    struct quant_state_vector *src_vec = &tester->v1;
    struct quant_state_vector *res_vec = &tester->res1;
    struct quant_state_vector *expected_vec = &tester->res2;

    double complex src[] = {
            0, 1, 0, 1
    };
    test_split_vector(tester, src, src_vec);

    double complex expected[] = {
            0, 1, 0, -1
    };
    test_split_vector(tester, expected, expected_vec);

    uint64_t control_qubit = 1;
    uint64_t target_qubit = 2;

    int8_t err = 0;
    err = quant_transform_crot(src_vec, res_vec, control_qubit, target_qubit);
    if (test_assert(tester, err, 0, "no error expected")) {
        return err;
    }

    err = test_assert_vector(tester, res_vec, expected_vec, "wrong result");
    return err;
}

// test case funcs =============================================================

void test_all(void)
{
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    struct test_runner tester = {0};
    tester.my_rank = my_rank;
    tester.world_size = world_size;

    int8_t err = test_assert(&tester, world_size <= 4, 1, "1 <= world size <= 4 expected");
    if (err) {
        return;
    }

    const int32_t suites_num = 14;
    struct test_suite suites[suites_num];
    suites[0] = (struct test_suite) {
            .name       = "Test hadamard transformation blackbox",
            .test       = test_hadamard_blackbox,
            .qubits_num = 5,
    };
    suites[1] = (struct test_suite) {
            .name       = "Test hadamard transformation whitebox",
            .test       = test_hadamard_whitebox,
            .qubits_num = 2,
    };
    suites[2] = (struct test_suite) {
            .name       = "Test hadamard^N transformation blackbox",
            .test       = test_hadamard_n_blackbox,
            .qubits_num = 6,
    };
    suites[3] = (struct test_suite) {
            .name       = "Test hadamard^N transformation whitebox",
            .test       = test_hadamard_n_whitebox,
            .qubits_num = 2,
    };
    suites[4] = (struct test_suite) {
            .name       = "Test NOT gate blackbox",
            .test       = test_not_blackbox,
            .qubits_num = 6,
    };
    suites[5] = (struct test_suite) {
            .name       = "Test NOT gate whitebox",
            .test       = test_not_whitebox,
            .qubits_num = 2,
    };
    suites[6] = (struct test_suite) {
            .name       = "Test CNOT gate blackbox",
            .test       = test_cnot_blackbox,
            .qubits_num = 5,
    };
    suites[7] = (struct test_suite) {
            .name       = "Test CNOT gate whitebox, no invert",
            .test       = test_cnot_whitebox_no_invert,
            .qubits_num = 3,
    };
    suites[8] = (struct test_suite) {
            .name       = "Test CNOT gate whitebox, with invert",
            .test       = test_cnot_whitebox_invert,
            .qubits_num = 3,
    };
    suites[9] = (struct test_suite) {
            .name       = "Test ROT gate blackbox",
            .test       = test_rot_blackbox,
            .qubits_num = 9,
    };
    suites[10] = (struct test_suite) {
            .name       = "Test ROT gate whitebox",
            .test       = test_rot_whitebox,
            .qubits_num = 3,
    };
    suites[11] = (struct test_suite) {
            .name       = "Test CROT gate blackbox",
            .test       = test_crot_blackbox,
            .qubits_num = 8,
    };
    suites[12] = (struct test_suite) {
            .name       = "Test CROT gate whitebox, no rotation",
            .test       = test_crot_whitebox_no_rotation,
            .qubits_num = 3,
    };
    suites[13] = (struct test_suite) {
            .name       = "Test CROT gate whitebox, with rotation",
            .test       = test_crot_whitebox_rotation,
            .qubits_num = 2,
    };

    int32_t errs = test_run_suites(&tester, suites, suites_num);
    int32_t passed = suites_num - errs;
    test_print(&tester, "%d/%d tests passed.\n", passed, suites_num);
}
