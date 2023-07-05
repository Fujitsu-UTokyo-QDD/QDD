# QDD

## Build and Test
You can build the repository as follows.
```
$ CC=mpicc CXX=mpicxx Boost_DIR=/usr/lib/x86_64-linux-gnu/cmake cmake -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build -j
```
BUILD_TYPE can be either `Release`, `Debug` or `RelWithDebInfo`.

You can run the test as follows.
```
$ ./build/test/qdd_test
```

## Python Bindings
After the above commands, try the following commands.
```
$ cp build/qdd/pyQDD.*.so .
$ poetry install
$ poetry run pytest
$ poetry build
```
You can find installable files in 'dist' directory.

## MPI Test
You can run the MPI programs as follows.
```
$ mpirun -np 4 ./build/test/mpt_test
$ mpirun -np 4 ./build/test/mpi_test_grover 20
```
Currently, python bindings does NOT support MPI.
