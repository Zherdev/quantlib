#include "quant.h"
#include "args.h"

#include <stdio.h>
#include <stdint.h>
#include <omp.h>

static void main_print_test_info(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *w);

// Usage:
// $ hadamard N K thread_num [-t]
int main(int argc, char *argv[])
{
    printf("Quantum hadamard transformation program started.\n");

    struct args args = {0};
    int8_t err = args_parse(&args, argc, argv);
    if (err) {
        printf("Bad args.\n");
        return -1;
    }

    omp_set_num_threads(args.thread_num);

    struct quant_state_vector v = {0};
    quant_state_vector_alloc(&v, args.n);
    quant_state_vector_init_random(&v);

    struct quant_state_vector w = {0};
    quant_state_vector_alloc(&w, args.n);

    quant_matrix u = {{0}};
    quant_matrix_init_hadamard(u);

    double start = omp_get_wtime();

    err = quant_hadamard_transform(&args, &v, &w, u);
    if (err) {
        printf("Transfrom failed.\n");
    }

    double end = omp_get_wtime();
    printf("Time elapsed: %lf\n", end - start);

    if (args.test) {
        main_print_test_info(&args, &v, &w);
    }

    quant_state_vector_free(&v);
    quant_state_vector_free(&w);

    return 0;
}

static void main_print_test_info(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *w)
{
    printf("\nArgs:\n");
    args_print(args);

    printf("\nState vector v:\n");
    quant_state_vector_print_info(v);
    quant_state_vector_print_elems(v);

    printf("\nState vector w:\n");
    quant_state_vector_print_info(w);
    quant_state_vector_print_elems(w);
}