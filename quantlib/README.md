quantlib
========

Parallel C library for quantum calculations with MPI support
------------------------------------------------------------

* Provides tools for quantum state manipulations and N-qubit quantum transformations.
* Such common gates as Hadamard, NOT, CNOT, ROT and etc are provided.
* Transformations with custom matrices are available also.
* Quantum state vectors can be stored & loaded as binary files.
* Non-parallel mode is available with one MPI process.

Compilation & Testing
---------------------

```bash
$ cd quantlib
$ make
$ make test-run-mpiexec
```

* Static library in `quantlib/lib/bin/quantlib.a`.
* Headers in `quantlib/lib/include/quantlib.h`.
* See `quantlib/Makefile` for additional info.

Examples
--------

See test suites in `quantlib/test/src/test.c` source codes for quantlib usage examples.

```c
#include "quantlib.h"

...

// Hadamard transformation example.

const uint64_t qubits_num = 5;
const uint64_t target_qubit = 1;

struct quant_state_vector v = {0};
quant_state_vector_alloc(&v, qubits_num);
quant_state_vector_fill_random(&v);

struct quant_state_vector res = {0};
quant_state_vector_alloc(&res, qubits_num);

int8_t err = quant_transform_hadamard(&v, &res, target_qubit);
if (err) {
    return err;
}

quant_state_vector_free(&v);
quant_state_vector_free(&res);
```

Notes
-----

* MPI-processes number are expected to be 2^N, N >= 0.
* N may be equal to zero, so it is possible to run calculations in the non-parallel mode.

Performance
-----------

Performance was tested on CMC MSU Polus IBM Power 8 high performance cluster.

See `quantlib/benchmarking` for details and benchmarking utils.

![Hadamard^N gate performance](https://github.com/Zherdev/quantum/blob/master/quantlib/benchmarking/hadamard_n.png)

Author
------

Zherdev, 2021.
