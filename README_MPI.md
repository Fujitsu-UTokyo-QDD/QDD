## Installing qdd-mpi from PyPI

qdd-mpi is distributed as an sdist on PyPI. The package name is different from the standard qdd package, but both provide the same `import qdd` interface. They are mutually exclusive -- do not install both into the same environment.

### Prerequisites

- C++17 compatible compiler
- MPI implementation (OpenMPI or MPICH)
- Boost with MPI and serialization components

On Debian/Ubuntu:

```sh
sudo apt-get install -y libopenmpi-dev openmpi-bin libboost-mpi-dev libboost-serialization-dev
```

### Install

```sh
CC=mpicc CXX=mpicxx pip install qdd-mpi --no-binary qdd-mpi
```

The `--no-binary qdd-mpi` flag is defensive (only sdist is published) but makes the intent explicit.

### Mutual exclusion with qdd

qdd and qdd-mpi install to the same `qdd` Python package namespace. If you want to switch between them, uninstall the other first:

```sh
pip uninstall qdd        # before installing qdd-mpi
pip uninstall qdd-mpi    # before installing qdd
```

## Manual build
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

Python bindings with MPI:
```sh
$ CMAKE_ARGS="-DisMPI=ON" CC=mpicc CXX=mpicxx pip install .
$ mpirun -np 2 python -m pytest test/python/test_mpi.py
```
