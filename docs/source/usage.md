# Usage

## Backends

`QddProvider` exposes the available QDD backends.

```python
from qdd import QddProvider

provider = QddProvider()
backend = provider.get_backend("qasm_simulator")
statevector_backend = provider.get_backend("statevector_simulator")
```

The default backend name is `qasm_simulator`.

## QDD Sampler

QDD also provides a Sampler implementation.

```python
from qiskit import QuantumCircuit
from qdd.qdd_sampler import Sampler

circuit = QuantumCircuit(2)
circuit.h(0)
circuit.cx(0, 1)
circuit.measure_all()

sampler = Sampler()
job = sampler.run([circuit])
result = job.result()
```

Pass backend options, transpile options, or skip Qiskit transpilation when the
input circuit is already compatible with QDD.

```python
sampler = Sampler(
    backend_options={"shots": 4096},
    transpile_options={"optimization_level": 1},
)
```

## QDD Estimator

The Estimator computes expectation values for observables represented by
Qiskit primitive inputs.

```python
from qdd.qdd_estimator import Estimator

estimator = Estimator(default_precision=0.0)
job = estimator.run([(circuit, observable)])
result = job.result()
```

`default_precision=0.0` requests exact probabilities. Positive precision values
are translated into a shot count.
