#!/bin/bash -x
TIMEOUT=60000
THREADS=16
BITS=33
GCFREQ=1000000

for BITS in $(seq 36 3 40)
do

# NO CACHE
#    SINGLE THREAD
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_RULE_MESSAGES=OFF; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_test $THREADS $BITS $GCFREQ
#    THREAD POOL (fiber)
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_threadpool=on; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_test $THREADS $BITS $GCFREQ
#    OPENMP
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_openmp=on; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_test $THREADS $BITS $GCFREQ

# LOCAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on &> /dev/null ; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_test $THREADS $BITS 1000000
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_test $THREADS $BITS $GCFREQ
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_test $THREADS $BITS $GCFREQ

# GLOBAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_test $THREADS $BITS $GCFREQ
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_test $THREADS $BITS $GCFREQ

done
