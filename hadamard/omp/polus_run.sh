#!/bin/bash
#
# Script for running hadamard transformation on Polus HPC.
#

Ns=(20 24 28)
Ks=(3 3 3)
Ts=(1 2 4 8 16)

mkdir -p report

for i in ${!Ns[*]}
do
    N=${Ns[$i]}
    K=${Ks[$i]}

    for T in ${Ts[*]}
    do
        echo "bin/hadamard $N $K $T"
        mpisubmit-new.pl --stdout "report/$K-$N-$T.out" --stderr "report/$K-$N-$T.err" bin/hadamard -- $N $K $T
    done
done

echo "Tasks submited"

watch "grep -r Time report | sort"