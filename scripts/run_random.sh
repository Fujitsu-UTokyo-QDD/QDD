#!/bin/bash -x
TIMEOUT=60000
GATES=200

BITS=16

for THRESHOLD in $(seq 1 2 16)
do
#/home/swli0426/Code/QDD/build/test/main $BITS $THRESHOLD $GATES
echo " "
done


for BITS in $(seq 16 2 30)
do
/home/swli0426/Code/QDD/build/test/main $BITS 11 $GATES
echo " "
done

