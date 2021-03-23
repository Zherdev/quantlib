#!/bin/bash
#
# Script for testing Adamar transformation on Polus HPC with MPI.
#

N=16
K=10
Ps=(2 4 8)

mkdir -p report
mkdir -p data

rm -f data/output-$N.bin

echo "bin/adamar generate $N data/input-$N.bin"
mpisubmit-new.pl -p 16 "bin/adamar" --stdout "report/test-gen.out" -- generate $N data/input-$N.bin
while [ ! -f "report/test-gen.out" ]
do
    sleep 2
done
echo "Generation is done!"

echo "bin/adamar transfrom $N $K data/input-$N.bin data/output-$K-$N-1.bin"
mpisubmit-new.pl -p "1" --stdout "report/test-transform-$K-$N-1.out" --stderr "report/test-transform-$K-$N-1.err" bin/adamar -- transform $N $K data/input-$N.bin data/output-$K-$N-1.bin
while [ ! -f "report/test-transform-$K-$N-1.out" ]
do
    sleep 2
done
echo "First transformation is done!"

for P in ${Ps[*]}
do
    echo "bin/adamar transfrom $N $K data/input-$N.bin data/output-$K-$N-$P.bin"
    mpisubmit-new.pl -p "$P" --stdout "report/test-transform-$K-$N-$P.out" --stderr "report/test-transform-$K-$N-$P.err" bin/adamar -- transform $N $K data/input-$N.bin data/output-$K-$N-$P.bin
    while [ ! -f "report/test-transform-$K-$N-$P.out" ]
    do
        sleep 2
    done
    echo "Transformation on $P processes is done!"

    echo "bin/adamar cmp $N data/output-$K-$N-1.bin data/output-$K-$N-$P.bin"
    mpisubmit-new.pl -p "$P" --stdout "report/test-cmp-$K-$N-$P.out" --stderr "report/test-cmp-$K-$N-$P.err" bin/adamar -- cmp $N data/output-$K-$N-1.bin data/output-$K-$N-$P.bin
    while [ ! -f "report/test-cmp-$K-$N-$P.out" ]
    do
        sleep 2
    done
    cat "report/test-cmp-$K-$N-$P.out"
    echo "Comparison for $P processes is done!"
done

