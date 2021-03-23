#include "args.h"

#include <stdint.h>
#include <string.h>
#include <stdio.h>

static int8_t args_process_params(struct args *res, int argc, char *argv[])
{
    if (!res || !argv || argc < 4) {
        return -1;
    }

    if (sscanf(argv[1], "%lu", &res->n) != 1) {
        return -1;
    }

    if (sscanf(argv[2], "%lu", &res->k) != 1) {
        return -1;
    }

    if (sscanf(argv[3], "%u", &res->thread_num) != 1) {
        return -1;
    }

    if (argc >= 5 && !strcmp(argv[4], "-t")) {
        res->test = 1;
    } else {
        res->test = 0;
    }

    if (res->n == 0 || res->k == 0 ||
            res->k > res->n || res->thread_num == 0) {
        return -1;
    }

    return 0;
}

int8_t args_parse(struct args *res, int argc, char *argv[])
{
    const char *usage = "Usage:\n"
                        "$ adamar N K thread_num [-t]\n";

    int8_t err = args_process_params(res, argc, argv);
    if (err) {
        printf(usage);
    }

    return err;
}

void args_print(struct args *args)
{
    printf("\tN: %lu\n", args->n);
    printf("\tK: %lu\n", args->k);
    printf("\tThread num: %u\n", args->thread_num);
    printf("\tIs test: %u\n", args->test);
}