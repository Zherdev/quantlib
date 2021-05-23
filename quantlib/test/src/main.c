#include "test.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

static void main_print(int my_rank, const char *fmt, ...);
static void main_init_random(void);
static int8_t main_check_world_size(void);

// Usage:
// $ test
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    main_init_random();

    int8_t err = main_check_world_size();
    if (err) {
        return -1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    test_all();
    MPI_Finalize();

    return err;
}

static void main_print(int my_rank, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    if (my_rank == 0) {
        vprintf(fmt, args);
    }

    va_end(args);
}

static void main_init_random(void)
{
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    srand((unsigned int) (time(NULL) * my_rank));
}

static int8_t main_check_world_size(void)
{
    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double log_size = log2(comm_size);
    if (log_size != floor(log_size)) {
        main_print(my_rank, "MPI comm world size is expected to be 2^m.\n");
        return -1;
    }

    return 0;
}
