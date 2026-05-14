# Usage

## Backends

`QddProvider` exposes QDD as Qiskit `BackendV2` objects. The default backend
is `qasm_simulator`, which returns shot counts or exact probabilities. The
`statevector_simulator` backend returns statevectors.

```python
from qdd import QddProvider

provider = QddProvider()
backend = provider.get_backend("qasm_simulator")
statevector_backend = provider.get_backend("statevector_simulator")
```

If the backend name is omitted, `QddProvider().get_backend()` returns
`qasm_simulator`.

## Run a Qiskit Circuit

Use the backend with normal Qiskit `QuantumCircuit` objects. Transpiling for the
QDD backend is recommended when the circuit may contain gates outside QDD's
native basis.

```python
from qiskit import QuantumCircuit, transpile

from qdd import QddProvider

backend = QddProvider().get_backend("qasm_simulator")

circuit = QuantumCircuit(3)
circuit.h(0)
circuit.cx(0, 1)
circuit.cx(1, 2)
circuit.measure_all()

compiled = transpile(circuit, backend, seed_transpiler=1234)
job = backend.run(
    compiled,
    shots=2048,
    memory=True,
    seed_simulator=1234,
)
result = job.result()

print(result.get_counts())
print(result.get_memory()[:5])
```

Set `shots=None` to compute exact probabilities instead of sampling:

```python
probability_circuit = circuit.remove_final_measurements(inplace=False)
compiled = transpile(probability_circuit, backend, seed_transpiler=1234)
result = backend.run(compiled, shots=None).result()
print(result.data(0)["probabilities"])
```

Runtime options can be passed to `run` or stored on the backend:

```python
backend.set_options(shots=4096, seed_simulator=1234)
result = backend.run(compiled).result()

result = backend.run(compiled, shots=1024, seed_simulator=5678).result()
```

Common runtime options are `shots`, `memory`, `seed_simulator`, `use_mpi`,
`show_progress`, and `show_progress_frequency`.

## Statevector Simulation

Use `statevector_simulator` for final statevectors. The circuit should not end
with measurements when you want the unmeasured quantum state.

```python
from qiskit import QuantumCircuit, transpile

from qdd import QddProvider

backend = QddProvider().get_backend("statevector_simulator")

circuit = QuantumCircuit(2)
circuit.h(0)
circuit.cx(0, 1)

compiled = transpile(circuit, backend, seed_transpiler=1234)
statevector = backend.run(compiled).result().get_statevector()
print(statevector)
```

## QDD Sampler

QDD also provides a Sampler implementation for Qiskit primitive-style inputs.
Each pub is a tuple of `(circuit, parameter_values)`.

```python
from qiskit import QuantumCircuit

from qdd.qdd_sampler import Sampler

circuit = QuantumCircuit(2)
circuit.h(0)
circuit.cx(0, 1)
circuit.measure_all()

sampler = Sampler()
job = sampler.run(pubs=[(circuit, [])], shots=4096)
result = job.result()[0]

print(result.data.quasi_dist.binary_probabilities())
print(result.metadata["shots"])
```

Use `is_exact=True` to request exact probabilities with `shots=None`:

```python
exact_job = sampler.run(pubs=[(circuit, [])], is_exact=True)
exact_result = exact_job.result()[0]
print(exact_result.data.quasi_dist.binary_probabilities())
```

Pass backend options, transpile options, or skip Qiskit transpilation in the
constructor. For example, this sets QDD backend shots and asks Qiskit
transpilation to use optimization level 1:

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
from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import SparsePauliOp

from qdd.qdd_estimator import Estimator

ansatz = RealAmplitudes(num_qubits=2, reps=1)
observable = SparsePauliOp.from_list([("ZZ", 1.0), ("XI", 0.5)])
theta = [0.1, 0.2, 0.3, 0.4]

estimator = Estimator(default_precision=0.0)
job = estimator.run([(ansatz, [observable], theta)])
result = job.result()[0]

print(result.data.evs)
```

`default_precision=0.0` requests exact probabilities. Positive precision values
are translated into a shot count. You can also override precision on each
`run` call:

```python
job = estimator.run([(ansatz, [observable], theta)], precision=0.01)
result = job.result()[0]
print(result.data.evs)
print(result.data.stds)
```

Backend options are accepted by both primitives. With the `qdd-mpi` package,
pass `{"use_mpi": True}` and launch the Python process with `mpirun`; see
{doc}`mpi` for a complete MPI example.
