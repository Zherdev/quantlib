#ifndef ARGS_H
#define ARGS_H

#include <stdint.h>

enum args_mode {
    ARGS_MODE_GENERATE,
    ARGS_MODE_TRANSFORM,
    ARGS_MODE_CMP
};

struct args {
    enum args_mode mode;
    uint64_t  qubits_num;
    uint64_t  target_qubit_num;
    int32_t   threads_num;
    double    accuracy;
    uint8_t   test;
    uint8_t   loss;
    char     *input_filename;
    char     *output_filename;
    char     *a_filename;
    char     *b_filename;
};

int8_t args_parse(struct args *res, int argc, char *argv[]);

int32_t args_snprint(char *buff, int32_t buff_size, struct args *args);

#endif // ARGS_H
