[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://fujitsu-utokyo-qdd.github.io/QDD/)

# QDD
[Documentation](https://fujitsu-utokyo-qdd.github.io/QDD/) is available on GitHub Pages.

QDD is a decision diagram based quantum computing simulator for Qiskit. It can
reduce memory usage compared with typical state-vector simulators when the
decision diagram representation remains compact.

## Installation

Install the standard package from PyPI:

```sh
pip install qdd
```

Install the MPI-enabled package from PyPI:

```sh
CC=mpicc CXX=mpicxx pip install qdd-mpi --no-binary qdd-mpi
```

The `qdd` and `qdd-mpi` distributions provide the same `import qdd` Python
package namespace. Do not install both in the same environment.

Supported environment:

- Linux x86_64 and aarch64
- Python 3.9 through 3.13

## Quick Start

QDD works as a Qiskit backend.

```python
from qiskit import QuantumCircuit
from qiskit.primitives import BackendSampler

from qdd import QddProvider

backend = QddProvider().get_backend()

circuit = QuantumCircuit(3)
circuit.h(0)
circuit.cx(0, 1)
circuit.measure_all()

sampler = BackendSampler(backend=backend)
job = sampler.run(circuits=circuit)
print(job.result())
```

For statevector simulation:

```python
backend = QddProvider().get_backend("statevector_simulator")
```

## Documentation

Full documentation, including source builds, testing, MPI usage, and Python/C++
API references, is published with GitHub Pages:

https://fujitsu-utokyo-qdd.github.io/QDD/

The documentation source is in `docs/source`.

## License

BSD 3-Clause Clear License

## Limitation of Liability

In no event and under no legal theory, whether in tort (including negligence), contract, or otherwise, unless required by applicable law (such as deliberate and grossly negligent acts) or agreed to in writing, shall any Contributor be liable to You for damages, including any direct, indirect, special, incidental, or consequential damages of any character arising as a result of this License or out of the use or inability to use the Work (including but not limited to damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses), even if such Contributor has been advised of the possibility of such damages.
