#include "quant.h"
#include "args.h"

#include <complex.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

void quant_state_vector_alloc(struct quant_state_vector *v, uint64_t n)
{
    v->n = n;
    v->len = pow(2, n);
    v->elems = calloc(v->len, sizeof(float complex));
}

void quant_state_vector_init_random(struct quant_state_vector *v)
{
    printf("Random vector initialization started...\n");

    const int max_val = 5;
    float sum = 0;

    #pragma omp parallel
    {
        // random initialization
        unsigned int seed = omp_get_thread_num();
        #pragma omp for
        for (uint64_t i = 0; i < v->len; i++) {
            float re = rand_r(&seed) % max_val;
            float im = rand_r(&seed) % max_val;

            v->elems[i] = re + im * I;
        }

        // normalization
        #pragma omp for reduction(+ : sum)
        for (uint64_t i = 0; i < v->len; i++) {
            float complex x = v->elems[i] * conjf(v->elems[i]);
            sum += crealf(x);
        }
        float norm = sqrtf(sum);
        #pragma omp for
        for (uint64_t i = 0; i < v->len; i++) {
            v->elems[i] = v->elems[i] / norm;
        }
    }

    printf("Random vector initialization is done.\n");
}

void quant_state_vector_print_info(struct quant_state_vector *v)
{
    printf("\tQubits num: %lu\n"
           "\tQuantum state vector len: %lu\n"
           "\tFloat complex size: %lu\n"
           "\tQuantum state vector size: %lu\n",
           v->n, v->len, sizeof(float complex), v->len * sizeof(float complex));
}

void quant_state_vector_print_elems(struct quant_state_vector *v)
{
    for (uint64_t i = 0; i < v->len; i++) {
        printf("\t%5lu: (%.4f, %.4f)\n", i, crealf(v->elems[i]), cimagf(v->elems[i]));
    }
}

void quant_state_vector_free(struct quant_state_vector *v)
{
    if (!v) {
        return;
    }

    free(v->elems);
    v->elems = NULL;
    v->n = 0;
    v->len = 0;
}

void quant_matrix_init_adamar(quant_matrix u)
{
    const float v = 1 / sqrtf(2);
    u[0][0] = v;
    u[0][1] = v;
    u[1][0] = v;
    u[1][1] = -v;
}

int8_t quant_adamar_transform(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *w,
        quant_matrix u)
{
    if (!v || !u || !w || !args) {
        return -1;
    }

    if (v->len != w->len || v->n != w->n) {
        return -1;
    }

    uint64_t len = v->len;
    uint64_t n = args->n;
    uint64_t k = args->k;

    printf("Quant transform started...\n");

    #pragma omp parallel for default(none) shared(v, w, u, len, n, k)
    for (uint64_t i = 0; i < len; i++) {
        uint8_t ik = (i >> (n - k)) & 1LL;

        uint64_t ik_one_mask = 1LL << (n - k);
        uint64_t ik_zero_mask = ~ik_one_mask;
        uint64_t i_zero = i & ik_zero_mask;
        uint64_t i_one = i_zero | ik_one_mask;

        w->elems[i] = u[ik][0] * v->elems[i_zero] + u[ik][1] * v->elems[i_one];
    }

    printf("Quant transform is done.\n");

    return 0;
}