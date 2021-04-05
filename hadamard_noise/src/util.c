#include "util.h"

#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// N(0, 1)
double util_normal_distribution(void)
{
    double sum = 0;
    const int32_t iters = 100;
    const double max = 1000;
    const double mean = max / 2;
    const double var = (max * max) / 12;

    for (int32_t i = 0; i < iters; i++) {
        sum += rand() % (int) max;
    }

    double from_normal = (sum - mean * iters) / sqrt(var * iters);
    return from_normal;
}

void util_srand_mpi(void)
{
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    char proc_name[MPI_MAX_PROCESSOR_NAME];
    int name_len = 0;

    MPI_Get_processor_name(proc_name, &name_len);

    srand((unsigned int) (time(NULL) + my_rank * comm_size + name_len));
}

unsigned int util_seed_omp(int seed_mpi)
{
    const int threads = omp_get_num_threads();
    const int thread = omp_get_thread_num();
    unsigned int seed_omp = time(NULL) + thread * threads + seed_mpi;
    return seed_omp;
}
