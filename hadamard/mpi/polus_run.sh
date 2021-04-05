#!/bin/bash
#
# Script for running hadamard transformation on Polus HPC with MPI.
#

Ns=(25 25 25  26 26 26  27 27 27)
Ks=(4  1  25  4  1  26  4  1  27)
Ps=(1 2 4 8 16)

mkdir -p report

for i in ${!Ns[*]}
do
    N=${Ns[$i]}
    K=${Ks[$i]}

    for P in ${Ps[*]}
    do
        echo "bin/hadamard transfrom $N $K data/input-$N.bin data/output-$N.bin"
        mpisubmit-new.pl -p "$P" --stdout "report/$K-$N-$P.out" --stderr "report/$K-$N-$P.err" bin/hadamard -- transform $N $K data/input-$N.bin data/output-$N.bin
        while [ ! -f "report/$K-$N-$P.out" ]
        do
            sleep 2
        done
        cat "report/$K-$N-$P.out"
    done
done
