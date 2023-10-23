# QDD

## Build and Test
You can build the repository as follows.
```
$ cmake -B build -DCMAKE_BUILD_TYPE=Release
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
$ cp build/qdd/pyQDD.*.so qdd
$ poetry install
$ poetry run pytest
$ poetry build
```
You can find installable files in 'dist' directory.

## MPI
You need to add some options.
```
$ CC=mpicc CXX=mpicxx Boost_DIR=/usr/lib/x86_64-linux-gnu/cmake cmake -B build -DCMAKE_BUILD_TYPE=Release -DisMPI=ON
$ cmake --build build -j
```

You can run the MPI programs as follows.
```
$ mpirun -np 4 ./build/test/mpt_test
$ mpirun -np 4 ./build/test/mpi_test_grover 20
```
Currently, python bindings does NOT support MPI.

# Limitation of Liability
In no event and under no legal theory, whether in tort (including negligence), contract, or otherwise, unless required by applicable law (such as deliberate and grossly negligent acts) or agreed to in writing, shall any Contributor be liable to You for damages, including any direct, indirect, special, incidental, or consequential damages of any character arising as a result of this License or out of the use or inability to use the Work (including but not limited to damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses), even if such Contributor has been advised of the possibility of such damages.
