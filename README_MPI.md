## MPI
You need to add some options.
```sh
$ CC=mpicc CXX=mpicxx Boost_DIR=/usr/lib/x86_64-linux-gnu/cmake cmake -B build -DCMAKE_BUILD_TYPE=Release -DisMPI=ON
$ cmake --build build -j
```

You can run the MPI programs as follows.
```sh
$ mpirun -np 4 ./build/test/mpt_test
$ mpirun -np 4 ./build/test/mpi_test_grover 20
```
Currently, python bindings does NOT support MPI.
