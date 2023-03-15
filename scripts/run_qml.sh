#!/bin/bash -x
TIMEOUT=60000
THREADS=16
GATES=100

for BITS in $(seq 10 2 30)
do

	python3 scripts/gates_test.py $BITS $GATES

done
