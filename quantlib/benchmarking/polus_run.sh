#!/bin/bash
#
# Script for benchmarking on IBM Power8 Polus HPC.
#

module load SpectrumMPI

Ns=(25 26 27)
Ps=(16 32 64)

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

cd hadamard
make clean && make
cd ../

for i in ${!Ns[*]}
do
    N=${Ns[$i]}

    for j in ${!Ps[*]}
    do
        P=${Ps[$j]}
        REPORT_OUT="report/hadamard-$P-$N.out"
        REPORT_ERR="report/hadamard-$P-$N.err"

        echo "mpisubmit.pl -p "$P" --stdout "$REPORT_OUT" --stderr "$REPORT_ERR"  hadamard/hadamard -- $N data/$N.bin"
        mpisubmit.pl -p "$P" --stdout "$REPORT_OUT" --stderr "$REPORT_ERR"  hadamard/hadamard -- $N data/$N.bin
        while [ ! -f "$REPORT_OUT" ]
        do
            sleep 2
        done
    done
done
