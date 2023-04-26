# QDD

## Build and Test
You can build the repository as follows.
```
$ Boost_DIR=/usr/lib/x86_64-linux-gnu/cmake cmake -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build -j
```
BUILD_TYPE can be either `Release`, `Debug` or `RelWithDebInfo`.

You can run the test as follows.
```
$ mpirun -np 4 ./build/test/mpi_test
```