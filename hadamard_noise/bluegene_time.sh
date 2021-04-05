#!/bin/bash
#
# Script for running Hadamard transformation on Bluegene HPC.
#

E=0.01
QUBITS=28
THREADS=(1 2 4)
PROCS=(16)

mkdir -p report/out
mkdir -p report/err

for P in ${PROCS[*]}
do
    for T in ${THREADS[*]}
    do
        echo "P=$P"
        echo "T=$T"
        echo "bin/hadamard transform $QUBITS $T $E data/input-$QUBITS.bin data/output-$QUBITS.bin"
        mpisubmit.bg -n "$P" --wtime 00:03:00 --stdout "report/out/$QUBITS-$P-$T.out" --stderr "report/err/$QUBITS-$P-$T.err" bin/hadamard -- transform $QUBITS $T $E data/input-$QUBITS.bin data/output-$QUBITS.bin
        while ! grep "^Done." "report/out/$QUBITS-$P-$T.out" > /dev/null 2>&1;
        do
            sleep 1
        done
        cat "report/out/$QUBITS-$P-$T.out"
    done
done
