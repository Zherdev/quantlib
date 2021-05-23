#include "quantlib.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <mpi.h>

/*
 * Util for running Hadamard^N transformation for quantlib benchmarking.
 *
 * Zherdev, 2021.
 */

static int8_t hadamard_n(uint64_t qubits_num, char *filename, int32_t my_rank);

// Usage:
// $ generate qubits_num input_filename
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc != 3) {
        return -1;
    }
    uint64_t qubits_num;
    if (sscanf(argv[1], "%" PRId64, &qubits_num) != 1) {
        return -1;
    }
    char *input_filename = argv[2];

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int8_t err = hadamard_n(qubits_num, input_filename, my_rank);

    MPI_Finalize();
    return err;
}

static int8_t hadamard_n(uint64_t qubits_num, char *filename, int32_t my_rank)
{
    MPI_File fh = MPI_FILE_NULL;
    int8_t err = MPI_File_open(MPI_COMM_WORLD, filename,
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err) {
        return err;
    }

    struct quant_state_vector v = {0};
    quant_state_vector_alloc(&v, qubits_num);

    struct quant_state_vector res = {0};
    quant_state_vector_alloc(&res, qubits_num);

    err = quant_state_vector_read_file(&v, fh);
    MPI_File_close(&fh);
    if (err) {
        quant_state_vector_free(&v);
        quant_state_vector_free(&res);
        return err;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    err = quant_transform_hadamard_n(&v, &res);
    if (!err) {
        double elapsed_local = MPI_Wtime() - start;
        double elapsed_max = 0;
        MPI_Reduce(&elapsed_local, &elapsed_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (my_rank == 0) {
            printf("Time elapsed:\n%lf\n", elapsed_max);
        }
    }

    quant_state_vector_free(&v);
    quant_state_vector_free(&res);

    return err;
}
