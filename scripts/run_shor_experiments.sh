#!/bin/bash -x
TIMEOUT=600000
THREADS=16
GCFREQ=1000000

for NUM in 12
do

# NO CACHE
#    SINGLE THREAD
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_RULE_MESSAGES=OFF; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_shor $NUM  $GCFREQ
#    THREAD POOL (shor)
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_threadpool=on; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_shor $NUM  $GCFREQ
#    OPENMP
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_openmp=on; cmake --build build/ -j; time timeout $TIMEOUT build/test/fiber_shor $NUM  $GCFREQ

# LOCAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on &> /dev/null ; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_shor $NUM  $GCFREQ
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_shor $NUM  $GCFREQ
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_shor $NUM  $GCFREQ

# GLOBAL CACHE
rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_threadpool=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_shor $NUM  $GCFREQ
#rm -rf build/; cmake -B build/ -DCMAKE_BUILD_TYPE=RELEASE -Dwith_cache_global=on -Dwith_openmp=on &> /dev/null; cmake --build build/ -j &> /dev/null; time timeout $TIMEOUT ./build/test/fiber_shor $NUM  $GCFREQ

done
