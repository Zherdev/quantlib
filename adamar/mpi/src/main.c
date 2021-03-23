#include "quant.h"
#include "args.h"

#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <string.h>

static void main_print_test_info(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *w);

static int8_t main_transform(struct args *args, int32_t my_rank, int32_t comm_size);

static int8_t main_generate(struct args *args, int32_t my_rank, int32_t comm_size);

static int8_t main_cmp(struct args *args, int32_t my_rank, int32_t comm_size);

#define ROOT 0

#define main_print(my_rank, ...) \
        if ((my_rank) == ROOT) { \
            printf(__VA_ARGS__); \
        }

// Usage:
// $ adamar transform qubits_num target_qubit_num input_filename output_filename [-t]
// $ adamar generate qubits_num output_filename [-t]
// $ adamar cmp qubits_num a_filename b_filename
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    srand((unsigned int) (time(NULL) * my_rank));

    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    double log_size = log2(comm_size);
    if (log_size != floor(log_size)) {
        main_print(my_rank, "MPI comm world size is expected to be 2^m.\n");
        return -1;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    struct args args = {0};
    int8_t err = args_parse(&args, argc, argv);
    if (err) {
        main_print(my_rank, "Bad args.\n");
        return -1;
    }

    int8_t res = 0;
    switch (args.mode) {
        case ARGS_MODE_TRANSFORM:
            res = main_transform(&args, my_rank, comm_size);
            break;

        case ARGS_MODE_GENERATE:
            res = main_generate(&args, my_rank, comm_size);
            break;

        case ARGS_MODE_CMP:
            res = main_cmp(&args, my_rank, comm_size);
            break;

        default:
            main_print(my_rank, "Unknown mode.\n");
            res = -1;
            break;
    }

    MPI_Finalize();
    return res;
}

static int8_t main_generate(struct args *args, int32_t my_rank, int32_t comm_size)
{
    MPI_File fh = MPI_FILE_NULL;
    int8_t err = MPI_File_open(MPI_COMM_WORLD, args->output_filename,
            MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    if (err) {
        main_print(my_rank, "Can't open the output file.\n");
        return err;
    }

    struct quant_state_vector v = {0};
    quant_state_vector_alloc(&v, args->qubits_num);

    main_print(my_rank, "Random vector initialization started...\n");
    quant_state_vector_init_random(&v);
    main_print(my_rank, "Random vector initialization is done.\n");

    main_print(my_rank, "Random vector output started...\n");
    err = quant_state_vector_print_file(&v, fh);
    if (err) {
        main_print(my_rank, "Can't write the vector to the file\n");
        return err;
    }
    main_print(my_rank, "Random vector output is done.\n");

    err = MPI_File_close(&fh);
    if (err) {
        main_print(my_rank, "Can't close the output file\n");
        return err;
    }

    if (args->test) {
        main_print_test_info(args, &v, NULL);
    }

    quant_state_vector_free(&v);

    return 0;
}

static int8_t main_transform(struct args *args, int32_t my_rank, int32_t comm_size)
{
    MPI_File fh_in = MPI_FILE_NULL;
    int8_t err = MPI_File_open(MPI_COMM_WORLD, args->input_filename,
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_in);
    if (err) {
        main_print(my_rank, "Can't open the input file.\n");
        return err;
    }

    MPI_File fh_out = MPI_FILE_NULL;
    err = MPI_File_open(MPI_COMM_WORLD, args->output_filename,
            MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh_out);
    if (err) {
        main_print(my_rank, "Can't open the output file.\n");
        return err;
    }

    struct quant_state_vector v = {0};
    quant_state_vector_alloc(&v, args->qubits_num);

    main_print(my_rank, "Vector read started...\n");
    quant_state_vector_read_file(&v, fh_in);
    main_print(my_rank, "Vector read is done.\n");

    struct quant_state_vector w = {0};
    quant_state_vector_alloc(&w, args->qubits_num);

    quant_matrix u = {{0}};
    quant_matrix_init_adamar(u);

    double start = MPI_Wtime();

    err = quant_adamar_transform(args, &v, &w, u);
    if (err) {
        main_print(my_rank, "Transfrom failed.\n");
    }

    double elapsed_local = MPI_Wtime() - start;
    double elapsed_max = 0;
    MPI_Reduce(&elapsed_local, &elapsed_max, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
    main_print(my_rank, "Time elapsed:\n%lf\n", elapsed_max);

    err = quant_state_vector_print_file(&w, fh_out);
    if (err) {
        main_print(my_rank, "Can't write the w vector to the file\n");
        return err;
    }

    err = MPI_File_close(&fh_out);
    if (err) {
        main_print(my_rank, "Can't close the output file\n");
        return err;
    }

    err = MPI_File_close(&fh_in);
    if (err) {
        main_print(my_rank, "Can't close the input file\n");
        return err;
    }

    if (args->test) {
        main_print_test_info(args, &v, &w);
    }

    quant_state_vector_free(&v);
    quant_state_vector_free(&w);

    return 0;
}

static int8_t main_cmp(struct args *args, int32_t my_rank, int32_t comm_size)
{
    MPI_File fh_a = MPI_FILE_NULL;
    int8_t err = MPI_File_open(MPI_COMM_WORLD, args->a_filename,
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_a);
    if (err) {
        main_print(my_rank, "Can't open the A file.\n");
        return err;
    }

    MPI_File fh_b = MPI_FILE_NULL;
    err = MPI_File_open(MPI_COMM_WORLD, args->b_filename,
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_b);
    if (err) {
        main_print(my_rank, "Can't open the B file.\n");
        return err;
    }

    struct quant_state_vector a = {0};
    quant_state_vector_alloc(&a, args->qubits_num);

    struct quant_state_vector b = {0};
    quant_state_vector_alloc(&b, args->qubits_num);

    main_print(my_rank, "Vector A read started...\n");
    quant_state_vector_read_file(&a, fh_a);
    main_print(my_rank, "Vector A read is done.\n");

    main_print(my_rank, "Vector B read started...\n");
    quant_state_vector_read_file(&b, fh_b);
    main_print(my_rank, "Vector B read is done.\n");

    main_print(my_rank, "Compare vectors...\n");
    err = quant_state_vector_cmp(&a, &b);
    if (err) {
        main_print(my_rank, "Vectors are NOT equal!\n");
        return err;
    }
    main_print(my_rank, "Vectors are equal, OK!\n");

    err = MPI_File_close(&fh_a);
    if (err) {
        main_print(my_rank, "Can't close the A file\n");
        return err;
    }

    err = MPI_File_close(&fh_b);
    if (err) {
        main_print(my_rank, "Can't close the B file\n");
        return err;
    }

    quant_state_vector_free(&a);
    quant_state_vector_free(&b);

    return 0;
}

static void main_print_test_info(
        struct args *args,
        struct quant_state_vector *v,
        struct quant_state_vector *w)
{
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int32_t buff_size = 1024;
    char buff[buff_size];
    memset(buff, ' ', buff_size);

    const int32_t root_buff_size = buff_size * comm_size;
    char root_buff[root_buff_size];
    memset(root_buff, ' ', root_buff_size);

    int32_t offset = 0;

    offset += snprintf(buff + offset, buff_size - offset,
            "\n********\nTest info. My rank: %d\n", my_rank);

    offset += snprintf(buff + offset, buff_size - offset, "\nArgs:\n");
    offset += args_snprint(buff + offset, buff_size - offset, args);

    offset += snprintf(buff + offset, buff_size - offset, "\nState vector v:\n");
    offset += quant_state_vector_snprint_info(buff + offset, buff_size - offset, v);
    //offset += quant_state_vector_snprint_elems(buff + offset, buff_size - offset, v);

    offset += snprintf(buff + offset, buff_size - offset, "\nState vector w:\n");
    offset += quant_state_vector_snprint_info(buff + offset, buff_size - offset, w);
    //offset += quant_state_vector_snprint_elems(buff + offset, buff_size - offset, w);

    offset += snprintf(buff + offset, buff_size - offset, "\n********\n");

    MPI_Gather(&(buff[0]), buff_size, MPI_CHAR, &(root_buff[0]), buff_size, MPI_CHAR, ROOT, MPI_COMM_WORLD);

    fwrite(root_buff, sizeof(root_buff), 1, stdout);
}
