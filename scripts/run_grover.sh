#! /bin/bash

NWORKERS=16
BIN=$(dirname `pwd`)/build/test/benchmark

echo "Run with $NWORKERS workers"

for i in `seq 2 2 10`
do
    echo""
    echo""
    echo "$i qubits"    
    if ! [ -x "$(command -v taskset)" ]; then
        $BIN $i $NWORKERS
    else
        taskset -c 0-15 $BIN $i $NWORKERS
    fi
done
