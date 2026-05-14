[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://fujitsu-utokyo-qdd.github.io/QDD/)

# QDD
[Documentation](https://fujitsu-utokyo-qdd.github.io/QDD/) is available on GitHub Pages.

QDD is a decision diagram based quantum computing simulator for Qiskit. It can
reduce memory usage compared with typical state-vector simulators when the
decision diagram representation remains compact.

## Installation

Install QDD from PyPI:

```sh
pip install qdd
```

Supported environment:

- Linux x86_64 and aarch64
- Python 3.9 through 3.13

## Quick Start

QDD works as a Qiskit backend.

```python
from qiskit import QuantumCircuit

from qdd import QddProvider

backend = QddProvider().get_backend("qasm_simulator")

circuit = QuantumCircuit(3)
circuit.h(0)
circuit.cx(0, 1)
circuit.measure_all()

job = backend.run(circuit, shots=1024, seed_simulator=1234)
print(job.result().get_counts())
```

For statevector simulation:

```python
from qiskit import QuantumCircuit

from qdd import QddProvider

backend = QddProvider().get_backend("statevector_simulator")

circuit = QuantumCircuit(2)
circuit.h(0)
circuit.cx(0, 1)

statevector = backend.run(circuit).result().get_statevector()
print(statevector)
```

## MPI Usage

MPI support is optional. Use the MPI-enabled `qdd-mpi` distribution only when
you want to run QDD across MPI ranks.

The `qdd` and `qdd-mpi` distributions provide the same `import qdd` Python
package namespace. Do not install both in the same environment.

Install the MPI-enabled package from PyPI:

```sh
pip install mpi4py
CC=mpicc CXX=mpicxx pip install qdd-mpi --no-binary qdd-mpi
```

Minimal MPI run with Qiskit:

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

Run the script under MPI:

```sh
mpirun -n 2 python mpi_bell.py
```

You can also pass `use_mpi=True` per run, or through QDD primitive backend
options such as `Sampler(backend_options={"use_mpi": True})` and
`Estimator(backend_options={"use_mpi": True})`.

## Documentation

Full documentation, including source builds, testing, MPI usage, and Python/C++
API references, is published with GitHub Pages:

https://fujitsu-utokyo-qdd.github.io/QDD/

The documentation source is in `docs/source`.

## License

BSD 3-Clause Clear License

## Limitation of Liability

In no event and under no legal theory, whether in tort (including negligence), contract, or otherwise, unless required by applicable law (such as deliberate and grossly negligent acts) or agreed to in writing, shall any Contributor be liable to You for damages, including any direct, indirect, special, incidental, or consequential damages of any character arising as a result of this License or out of the use or inability to use the Work (including but not limited to damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses), even if such Contributor has been advised of the possibility of such damages.
