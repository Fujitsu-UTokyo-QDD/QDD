#!/bin/bash -x
TIMEOUT=60000
THREADS=16
BITS=33
GCFREQ=1000000
GATES=200

for BITS in $(seq 10 3 16)
do

# NO CACHE
#    SINGLE THREAD
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_RULE_MESSAGES=OFF; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES
#    THREAD POOL (fiber)
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_threadpool=on; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES
#    OPENMP
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_openmp=on; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES

# LOCAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on &> /dev/null ; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES

# GLOBAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_clifford $THREADS $BITS $GCFREQ $GATES

done
