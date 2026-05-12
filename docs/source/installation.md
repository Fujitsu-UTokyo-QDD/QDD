# Installation

## Standard Package

Install QDD from PyPI:

```sh
pip install qdd
```

The package exposes the Python package name `qdd`.

## MPI Package

Install the MPI-enabled package from PyPI:

```sh
CC=mpicc CXX=mpicxx pip install qdd-mpi --no-binary qdd-mpi
```

The standard `qdd` package and the MPI-enabled `qdd-mpi` package both provide
the same `import qdd` namespace. Do not install both in the same environment.
Uninstall the existing package before switching variants.

```sh
pip uninstall qdd
CC=mpicc CXX=mpicxx pip install qdd-mpi --no-binary qdd-mpi
```

See {doc}`mpi` for MPI prerequisites and build details.

## Supported Environment

QDD currently targets Linux on x86_64 and aarch64 with Python 3.9 through 3.13.
Source builds require a C++17 compiler and CMake 3.25 or newer.
