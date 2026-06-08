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

Python runs that set `use_mpi=True` also need `mpi4py` installed for the same
MPI implementation. See {doc}`mpi` for prerequisites, build details, and a
complete `mpirun -n 2 python ...` example.

## Supported Environment

| Platform | Architecture | Python |
|----------|-------------|--------|
| Linux    | x86_64, aarch64 | 3.9 – 3.14 |
| macOS 10.12 or newer | x86_64 (Intel) | 3.9 – 3.14 |
| macOS 11.0 or newer | arm64 (Apple Silicon) | 3.9 – 3.14 |
| Windows  | x86_64 (AMD64), arm64 (ARM64) | 3.9 – 3.14 |

Source builds require a C++17 compiler and CMake 3.25 or newer.

> **Note:** MPI support is available on Linux only. The `qdd-mpi` package is not
> available for macOS or Windows.
