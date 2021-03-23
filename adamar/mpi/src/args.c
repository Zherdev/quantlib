#include "args.h"

#include <stdint.h>
#include <string.h>
#include <stdio.h>

static uint8_t args_process_params_transform(struct args *res, int argc, char *argv[])
{
    if (argc < 6) {
        return -1;
    }

    if (sscanf(argv[2], "%lu", &res->qubits_num) != 1) {
        return -1;
    }

    if (sscanf(argv[3], "%lu", &res->target_qubit_num) != 1) {
        return -1;
    }

    res->input_filename = argv[4];
    res->output_filename = argv[5];

    if (argc >= 7 && !strcmp(argv[6], "-t")) {
        res->test = 1;
    } else {
        res->test = 0;
    }

    if (res->qubits_num == 0 ||
            res->target_qubit_num == 0 || res->target_qubit_num > res->qubits_num ||
            strlen(res->input_filename) == 0 ||
            strlen(res->output_filename) == 0) {
        return -1;
    }

    return 0;
}

static uint8_t args_process_params_generate(struct args *res, int argc, char *argv[])
{
    if (argc < 4) {
        return -1;
    }

    if (sscanf(argv[2], "%lu", &res->qubits_num) != 1) {
        return -1;
    }

    res->output_filename = argv[3];

    if (res->qubits_num == 0 || strlen(res->output_filename) == 0) {
        return -1;
    }

    if (argc >= 5 && !strcmp(argv[4], "-t")) {
        res->test = 1;
    } else {
        res->test = 0;
    }

    return 0;
}

static uint8_t args_process_params_cmp(struct args *res, int argc, char *argv[])
{
    if (argc != 5) {
        return -1;
    }

    if (sscanf(argv[2], "%lu", &res->qubits_num) != 1) {
        return -1;
    }

    res->a_filename = argv[3];
    res->b_filename = argv[4];

    if (res->qubits_num == 0 || strlen(res->a_filename) == 0 || strlen(res->b_filename) == 0) {
        return -1;
    }

    return 0;
}

static int8_t args_process_params(struct args *res, int argc, char *argv[])
{
    if (!res || !argv || !argc) {
        return -1;
    }

    if (!strcmp(argv[1], "generate")) {
        res->mode = ARGS_MODE_GENERATE;
        return args_process_params_generate(res, argc, argv);
    }

    if (!strcmp(argv[1], "transform")) {
        res->mode = ARGS_MODE_TRANSFORM;
        return args_process_params_transform(res, argc, argv);
    }

    if (!strcmp(argv[1], "cmp")) {
        res->mode = ARGS_MODE_CMP;
        return args_process_params_cmp(res, argc, argv);
    }

    return -1;
}

int8_t args_parse(struct args *res, int argc, char *argv[])
{
    const char *usage = "Usage:\n"
                        "$ adamar transform qubits_num target_qubit_num input_filename output_filename [-t]\n"
                        "$ adamar generate qubits_num output_filename [-t]\n"
                        "$ adamar cmp qubits_num a_filename b_filename\n";

    int8_t err = args_process_params(res, argc, argv);
    if (err) {
        printf("%s", usage);
    }

    return err;
}

int32_t args_snprint(char *buff, int32_t buff_size, struct args *args)
{
    int32_t offset = 0;

    offset += snprintf(buff + offset, buff_size - offset, "\tMode: %d\n", args->mode);
    offset += snprintf(buff + offset, buff_size - offset, "\tQubits num: %lu\n", args->qubits_num);
    offset += snprintf(buff + offset, buff_size - offset, "\tTarget qubit num: %lu\n", args->target_qubit_num);
    offset += snprintf(buff + offset, buff_size - offset, "\tInput file: %s\n", args->input_filename ?: "");
    offset += snprintf(buff + offset, buff_size - offset, "\tA file: %s\n", args->a_filename ?: "");
    offset += snprintf(buff + offset, buff_size - offset, "\tB file: %s\n", args->b_filename ?: "");
    offset += snprintf(buff + offset, buff_size - offset, "\tOutput file: %s\n", args->output_filename ?: "");
    offset += snprintf(buff + offset, buff_size - offset, "\tIs test: %u\n", args->test);

    return offset;
}