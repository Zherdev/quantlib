#!/bin/bash
#
# Script for benchmarking on IBM Power8 Polus HPC.
#

module load SpectrumMPI

Ns=(23 24 25)
Ps=(1 2 4 8 16 32)

mkdir -p report
mkdir -p data

# Generate data

cd generate
make clean && make
cd ../

for i in ${!Ns[*]}
do
    N=${Ns[$i]}
    P=16
    REPORT_OUT="report/generate-$N.out"
    REPORT_ERR="report/generate-$N.err"

    echo "mpisubmit.pl -p $P --stdout $REPORT_OUT --stderr $REPORT_ERR generate/generate -- $N data/$N.bin"
    mpisubmit.pl -p $P --stdout $REPORT_OUT --stderr $REPORT_ERR generate/generate -- $N data/$N.bin
    while [ ! -f "$REPORT_OUT" ]
    do
        sleep 2
    done
done

# Hadamard^N

cd hadamard_n
make clean && make
cd ../

for i in ${!Ns[*]}
do
    N=${Ns[$i]}

    for j in ${!Ps[*]}
    do
        P=${Ps[$j]}
        REPORT_OUT="report/hadamard_n-$P-$N.out"
        REPORT_ERR="report/hadamard_n-$P-$N.err"

        echo "mpisubmit.pl -p "$P" --stdout "$REPORT_OUT" --stderr "$REPORT_ERR"  hadamard_n/hadamard_n -- $N data/$N.bin"
        mpisubmit.pl -p "$P" --stdout "$REPORT_OUT" --stderr "$REPORT_ERR"  hadamard_n/hadamard_n -- $N data/$N.bin
        while [ ! -f "$REPORT_OUT" ]
        do
            sleep 2
        done
    done
done

# CNOT

cd cnot
make clean && make
cd ../

for i in ${!Ns[*]}
do
    N=${Ns[$i]}

    for j in ${!Ps[*]}
    do
        P=${Ps[$j]}
        REPORT_OUT="report/cnot-$P-$N.out"
        REPORT_ERR="report/cnot-$P-$N.err"

        echo "mpisubmit.pl -p "$P" --stdout "$REPORT_OUT" --stderr "$REPORT_ERR"  cnot/cnot -- $N data/$N.bin"
        mpisubmit.pl -p "$P" --stdout "$REPORT_OUT" --stderr "$REPORT_ERR"  cnot/cnot -- $N data/$N.bin
        while [ ! -f "$REPORT_OUT" ]
        do
            sleep 2
        done
    done
done
