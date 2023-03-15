#!/bin/bash -x
TIMEOUT=60000
THREADS=16
BITS=33
GCFREQ=1000000

for BITS in $(seq 31 1 34)
do

# NO CACHE
#    SINGLE THREAD
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_RULE_MESSAGES=OFF &> /dev/null; cmake --build build/ -j &> /dev/null ; time  build/test/fiber_test $THREADS $BITS $GCFREQ
#    THREAD POOL (fiber)
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time build/test/fiber_test $THREADS $BITS $GCFREQ
#    OPENMP
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  build/test/fiber_test $THREADS $BITS $GCFREQ

# LOCAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on &> /dev/null ; cmake --build build/ -j &> /dev/null; time ./build/test/fiber_test $THREADS $BITS 1000000
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GCFREQ
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GCFREQ

# GLOBAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on &> /dev/null ; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS 1000000
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GCFREQ
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time  ./build/test/fiber_test $THREADS $BITS $GCFREQ

done
