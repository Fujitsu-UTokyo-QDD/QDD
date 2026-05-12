# Quick Start

QDD can be used as a Qiskit backend.

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

To retrieve statevectors, use the statevector backend:

```python
from qdd import QddProvider

backend = QddProvider().get_backend("statevector_simulator")
```
