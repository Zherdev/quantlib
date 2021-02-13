#ifndef ARGS_H
#define ARGS_H

#include <stdint.h>

struct args {
    uint64_t n;
    uint64_t k;
    uint32_t thread_num;
    uint8_t  test : 1;
};

int8_t args_parse(struct args *res, int argc, char *argv[]);

void args_print(struct args *args);

#endif // ARGS_H