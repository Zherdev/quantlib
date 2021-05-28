#include "quantlib.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <mpi.h>

/*
 * Util for generating data for quantlib benchmarking.
 *
 * Zherdev, 2021.
 */

static int8_t generate(uint64_t qubits_num, char *filename);
static void util_srand(int32_t my_rank, int32_t comm_size);

// Usage:
// $ generate qubits_num output_filename
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
    char *output_filename = argv[2];

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    util_srand(my_rank, comm_size);

    int8_t err = generate(qubits_num, output_filename);

    MPI_Finalize();
    return err;
}

static int8_t generate(uint64_t qubits_num, char *filename)
{
    MPI_File fh = MPI_FILE_NULL;
    int8_t err = MPI_File_open(MPI_COMM_WORLD, filename,
            MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    if (err) {
        return err;
    }

    struct quant_state_vector v = {0};
    quant_state_vector_alloc(&v, qubits_num);
    quant_state_vector_fill_random(&v);

    err = quant_state_vector_print_file(&v, fh);

    MPI_File_close(&fh);
    quant_state_vector_free(&v);

    return err;
}

static void util_srand(int32_t my_rank, int32_t comm_size)
{
    srand((unsigned int) (time(NULL) + my_rank * comm_size));
}
