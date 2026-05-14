# Quick Start

QDD can be used as a Qiskit backend with normal Qiskit circuits.

```python
from qiskit import QuantumCircuit

from qdd import QddProvider

backend = QddProvider().get_backend("qasm_simulator")

circuit = QuantumCircuit(3)
circuit.h(0)
circuit.cx(0, 1)
circuit.cx(1, 2)
circuit.measure_all()

job = backend.run(circuit, shots=1024, seed_simulator=1234)
print(job.result().get_counts())
```

To retrieve statevectors, use the statevector backend:

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

For MPI execution, install `qdd-mpi`, set `use_mpi=True`, and launch the script
with `mpirun`; see {doc}`mpi`.
