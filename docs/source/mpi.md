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

Python MPI runs with `use_mpi=True` also need `mpi4py` installed against the
same MPI implementation used to build `qdd-mpi`.

```sh
pip install mpi4py
```

## Install From PyPI

```sh
CC=mpicc CXX=mpicxx pip install qdd-mpi --no-binary qdd-mpi
```

Use the MPI wrapper compilers so the source build uses the MPI implementation's
compiler and linker flags.
The `--no-binary qdd-mpi` option makes the source-build requirement explicit.

## Run Qiskit Code With MPI

The Python import name is still `qdd`. Enable MPI in QDD with `use_mpi=True`
and launch the script with `mpirun` or `mpiexec`.

```python
# mpi_bell.py
from mpi4py import MPI
from qiskit import QuantumCircuit

from qdd import QddProvider

backend = QddProvider().get_backend("qasm_simulator")
backend.set_options(use_mpi=True)

circuit = QuantumCircuit(3)
circuit.h(0)
circuit.cx(0, 1)
circuit.cx(1, 2)
circuit.measure_all()

result = backend.run(circuit, shots=1024).result()

if MPI.COMM_WORLD.Get_rank() == 0:
    print(result.get_counts())
```

Run it with multiple MPI ranks:

```sh
mpirun -n 2 python mpi_bell.py
```

You can set `use_mpi=True` per run instead of storing it on the backend:

```python
result = backend.run(circuit, shots=1024, use_mpi=True).result()
```

QDD primitives use the same backend option:

```python
from qdd.qdd_estimator import Estimator
from qdd.qdd_sampler import Sampler

sampler = Sampler(backend_options={"use_mpi": True})
estimator = Estimator(backend_options={"use_mpi": True})
```

Every MPI rank executes the Python script, so guard printing, file writes, and
other side effects with the MPI rank when they should happen only once.

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
