#!/bin/bash

for nproc in {1,2,4,8,16}; do
# for nproc in {1,2}; do
    mkdir $nproc
    cp *.c Makefile $nproc/
    cd $nproc
    make NPES=$nproc
    for mi in {100,500,1000,2000,3000,4000}; do
    # for mi in {100,500}; do
        if (($nproc == 1)); then
            ./laplace_serial -q -m $mi -s summary_$nproc.txt
        fi
        mpirun -n $nproc laplace_mpi -q -m $mi -s summary_$nproc.txt
        upcrun -n $nproc laplace_upc -q -m $mi -s summary_$nproc.txt
    done
    cd ..
done