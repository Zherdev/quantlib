#!/bin/bash
#
# Script for calculating precision loss for Hadamard transformation on Bluegene HPC.
#

E=0.01
QUBITS=(24 25 26 27 28)
THREADS=4
PROCS=256
RUNS=60

mkdir -p data
mkdir -p report/gen
mkdir -p report/out
mkdir -p report/err
mkdir -p report/loss

for Q in ${QUBITS[*]}
do
    for i in $(seq 1 $RUNS)
    do
        echo "Q=$Q"
        echo "i=$i"

        rm -f data/input-$Q.bin

        echo "bin/hadamard generate $Q $THREADS data/input-$Q.bin"
        mpisubmit.bg -n "$PROCS" --wtime 00:03:00 --stdout "report/gen/$Q-$PROCS-$THREADS-$i.out" --stderr "report/gen/$Q-$PROCS-$THREADS-$i.err" bin/hadamard -- generate $Q $THREADS data/input-$Q.bin
        while ! grep "^Done." "report/gen/$Q-$PROCS-$THREADS-$i.out" > /dev/null 2>&1;
        do
            sleep 1
        done

        echo "bin/hadamard transform $Q $THREADS $E data/input-$Q.bin data/output-$Q.bin -l"
        mpisubmit.bg -n "$PROCS" --wtime 00:03:00 --stdout "report/out/$Q-$PROCS-$THREADS-$i.out" --stderr "report/err/$Q-$PROCS-$THREADS-$i.err" bin/hadamard -- transform $Q $THREADS $E data/input-$Q.bin data/output-$Q.bin -l
        while ! grep "^Done." "report/out/$Q-$PROCS-$THREADS-$i.out" > /dev/null 2>&1;
        do
            sleep 1
        done

        cat "report/out/$Q-$PROCS-$THREADS-$i.out"
        loss=$(grep -oP '^Precision loss:\n[0-9]*\.[0-9]+$' report/out/$Q-$PROCS-$THREADS-$i.out | tail -n 1)
        echo $loss >> "report/loss/$Q-$PROCS-$THREADS.tmp"
    done

    cat "report/loss/$Q-$PROCS-$THREADS.tmp"
done
