#!/bin/bash -x
TIMEOUT=60000
THREADS=16
GATES=200

#echo "ddsim"

for BITS in $(seq 10 2 8)
do
/home/swli0426/Code/ddsim/build/apps/ddsim_random $BITS $GATES
done

for BITS in $(seq 24 2 14)
do

# NO CACHE
#    SINGLE THREAD
echo " "
echo "NO CACHE"
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_RULE_MESSAGES=OFF &> /dev/null; cmake --build build/ -j &> /dev/null ; time  build/test/fiber_test 1 $BITS $GATES
#    THREAD POOL (fiber)
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time build/test/fiber_test $THREADS $BITS $GATES
#    OPENMP
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  build/test/fiber_test $THREADS $BITS $GATES

echo " "
# LOCAL CACHE
echo "LOCAL CACHE"
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on &> /dev/null ; cmake --build build/ -j &> /dev/null; time ./build/test/fiber_test 1 $BITS $GATES
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GATES
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GATES

echo " "
# GLOBAL CACHE
echo "GLOBAL CACHE"
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on &> /dev/null ; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GATES
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GATES
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GATES

done

echo " "
echo "qiskit "
for BITS in $(seq 10 2 34)
do

	python3 scripts/gates_test.py $BITS $GATES

done

echo " "

