# MPI Support

QDD provides an MPI-enabled distribution named `qdd-mpi`. It installs the same
Python package namespace as the standard package, so the two distributions are
mutually exclusive in one environment.

## Prerequisites

Install an MPI implementation and Boost components for MPI and serialization.
On Debian or Ubuntu:

```sh
sudo apt-get install -y libopenmpi-dev openmpi-bin \
    libboost-mpi-dev libboost-serialization-dev
```

## Install From PyPI

```sh
CC=mpicc CXX=mpicxx pip install qdd-mpi --no-binary qdd-mpi
```

The `--no-binary qdd-mpi` option makes the source-build requirement explicit.

## Manual Build

```sh
CC=mpicc CXX=mpicxx Boost_DIR=/usr/lib/x86_64-linux-gnu/cmake \
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DisMPI=ON
cmake --build build -j
```

Python bindings with MPI:

```sh
CMAKE_ARGS="-DisMPI=ON" CC=mpicc CXX=mpicxx pip install .
mpirun -np 2 python -m pytest test/python/test_mpi.py
```
