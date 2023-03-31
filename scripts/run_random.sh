#!/bin/bash -x
TIMEOUT=60000
GATES=200

BITS=16




for BITS in $(seq 16 2 30)
do
/home/swli0426/Code/QDD/build/test/main $BITS  $GATES
echo " "
done

