# Performance Characteristics

QDD is a decision-diagram-based quantum circuit simulator. Its performance is
best when the simulated state and operators have enough repeated structure for
the decision diagram to share nodes. It can be slower when the circuit quickly
creates many distinct amplitudes or many distinct edge weights, because the
decision diagram then loses sharing and approaches the cost of an explicit
state-vector representation with additional bookkeeping overhead.

This page gives practical guidance for choosing workloads where QDD is likely
to work well.

## When QDD Tends To Be Fast

QDD is well suited to circuits whose state evolution remains compact as a
decision diagram:

- Grover-style search circuits with repeated oracle and diffusion blocks.
- Shor-style arithmetic circuits made from reversible modular arithmetic,
  controlled operations, and repeated structured subcircuits.
- Circuits with many Clifford, permutation, controlled-NOT, Toffoli, and other
  operations that preserve regular structure.
- Problems where a large Hilbert space is explored, but the amplitudes still
  fall into repeated patterns.

Grover and Shor circuits often have these properties. They may use many qubits
and many gates, but the circuit structure is highly regular, and intermediate
states often contain repeated subgraphs that a decision diagram can represent
compactly.

## When QDD Can Be Slow

QDD is less favorable for circuits that create many unrelated amplitudes:

- VQE circuits with hardware-efficient ansatze such as `EfficientSU2` or
  `RealAmplitudes` at large depth.
- Circuits with many independent small-angle rotations.
- Random or problem-instance-specific parameterized circuits with little
  repeated structure.
- Optimization loops that evaluate thousands of slightly different parameter
  points.

In these cases, each parameter value can introduce distinct edge weights. When
those weights are all different, the decision diagram cannot merge as many
nodes. A state-vector simulator may then be faster because it performs dense
linear algebra directly without maintaining the decision-diagram structure.

This does not mean VQE is unsupported. QDD can run VQE-style circuits, but it is
usually not the workload where QDD has the strongest advantage.

## Example: Grover With Qiskit's Grover Class

The examples in this section use optional Qiskit ecosystem packages such as
`qiskit-algorithms` and `qiskit-nature` in addition to QDD.

Qiskit Algorithms provides `Grover` and `AmplificationProblem`, so you do not
need to write the Grover iterator by hand. The example below uses Qiskit's
`PhaseOracleGate` to turn a 3-SAT instance into a phase oracle, then runs the
algorithm with QDD's Sampler.

```python
from qiskit.circuit.library.phase_oracle import PhaseOracleGate
from qiskit.synthesis.boolean.boolean_expression import BooleanExpression
from qiskit_algorithms import AmplificationProblem, Grover

from qdd.qdd_sampler import Sampler


input_3sat_instance = """
c example DIMACS-CNF 3-SAT
p cnf 3 5
-1 -2 -3 0
1 -2 3 0
1 2 -3 0
1 -2 -3 0
-1 2 3 0
"""

expression = BooleanExpression.from_dimacs(input_3sat_instance).expression
oracle = PhaseOracleGate(expression)
problem = AmplificationProblem(oracle)

grover = Grover(sampler=Sampler())
result = grover.amplify(problem)

print(result.assignment)
```

This is a favorable pattern for QDD because the algorithm repeatedly applies a
structured oracle and diffusion operator. The search space grows with the
number of variables, but the circuit can still contain repeated subgraphs that
a decision diagram can share.

The same example can use MPI by constructing the sampler with backend options
and launching the script with `mpirun`:

```python
grover = Grover(sampler=Sampler(backend_options={"use_mpi": True}))
```

```sh
mpirun -n 4 python grover_3sat.py
```

This example follows the public Qiskit Algorithms Grover and 3-SAT examples.

## Example: Shor From IBM Quantum's Tutorial

Shor's algorithm is another good match for a decision diagram simulator because
its core is order finding through reversible modular arithmetic. The circuit is
not necessarily small, but it is made from repeated arithmetic blocks,
controlled modular multiplication, and an inverse QFT.

Unlike Grover and VQE, Qiskit documents Shor as a tutorial circuit rather than
as a single main algorithm API class. For a generally available reference
implementation, use IBM Quantum's "Shor's algorithm" tutorial. It constructs
the order-finding circuit for `N = 15` and `a = 2` using Qiskit's standard
`QFT` and `UnitaryGate` building blocks, then runs the circuit with a Sampler
and classically post-processes the measured phase.

The QDD-facing part is the same as for any Qiskit circuit:

```python
from qiskit import transpile
from qiskit.circuit.library import QFT, UnitaryGate

from qdd import QddProvider

# Build the order-finding circuit as in IBM Quantum's Shor tutorial:
# https://qiskit.qotlabs.org/docs/tutorials/shors-algorithm
circuit = build_order_finding_circuit_for_N15_a2()

backend = QddProvider().get_backend("qasm_simulator")
compiled = transpile(circuit, backend, optimization_level=1)
counts = backend.run(compiled, shots=1024).result().get_counts()
print(counts)
```

The placeholder `build_order_finding_circuit_for_N15_a2()` is intentionally not
a QDD-specific helper. It should be replaced by the circuit construction from
the public IBM Quantum tutorial. That keeps the example independent of this
repository's benchmark code while documenting how to run the resulting circuit
with QDD.

## Example: H2 VQE With Qiskit Nature

Qiskit Nature provides chemistry-specific components for VQE. The example
below uses its standard H2 setup: `PySCFDriver` builds the electronic
structure problem, `JordanWignerMapper` maps it to qubits, and `UCCSD` with a
`HartreeFock` initial state provides a chemically motivated ansatz.

```python
import numpy as np
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import SLSQP
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper

from qdd.qdd_estimator import Estimator


driver = PySCFDriver(atom="H 0 0 0; H 0 0 0.735", basis="sto-3g")
problem = driver.run()

mapper = JordanWignerMapper()

ansatz = UCCSD(
    problem.num_spatial_orbitals,
    problem.num_particles,
    mapper,
    initial_state=HartreeFock(
        problem.num_spatial_orbitals,
        problem.num_particles,
        mapper,
    ),
)

vqe = VQE(Estimator(default_precision=0.0), ansatz, SLSQP())
vqe.initial_point = np.zeros(ansatz.num_parameters)

solver = GroundStateEigensolver(mapper, vqe)
result = solver.solve(problem)

print(result.total_energies[0])
```

This is a realistic and convenient way to express a molecular VQE problem, but
it is also representative of workloads where QDD may be slower than a dense
state-vector simulator. UCC-style and hardware-efficient ansatze contain many
parameterized rotations. During an optimizer loop, each evaluation uses a
different parameter vector, which often produces many distinct
decision-diagram edge weights and reduces node sharing.

For this class of workloads, use a state-vector simulator as a baseline and
choose QDD only when measurements show that the specific ansatz and Hamiltonian
remain compact.

## Practical Selection Guide

Use QDD first when:

- the circuit contains repeated structured blocks;
- the algorithm is search, arithmetic, or reversible logic heavy;
- you expect many amplitudes or suboperators to be identical or closely
  reusable;
- memory use is the limiting factor for a state-vector simulator.

Benchmark against a state-vector simulator when:

- the circuit is dominated by arbitrary `rx`, `ry`, `rz`, `u`, or controlled
  rotation gates with many distinct angles;
- the circuit resembles a random parameterized ansatz;
- the workload is an optimizer loop that repeatedly evaluates slightly
  different circuits;
- the decision diagram grows quickly during early test runs.

For a new workload, start with a smaller instance, compare QDD with a
state-vector backend, and scale only if the QDD run shows the expected compact
representation advantage.

## Benchmark Script

The repository includes `bench/compare_simulators.py` to compare QDD with
Qiskit Aer statevector and matrix-product-state simulation on the example
workloads above:

```sh
python bench/compare_simulators.py --benchmark grover --grover-qubits 20
python bench/compare_simulators.py --benchmark shor --counting-qubits 22
python bench/compare_simulators.py --benchmark vqe --molecule h4 --vqe-reps 8
python bench/compare_simulators.py --benchmark all
```

The VQE benchmark builds molecular Hamiltonians with Qiskit Nature. By default it
uses a deeper `EfficientSU2` ansatz to make the small-angle-rotation workload
large enough to time; pass `--vqe-ansatz uccsd` to use the chemistry-specific
UCCSD ansatz instead. Install the project with the `test` extra before running
the benchmark:

```sh
pip install .[test]
```

The defaults are intentionally larger than smoke tests so the timing
differences are visible. For a quick functional check, reduce the scale, for
example `--grover-qubits 12`, `--counting-qubits 8`, or `--vqe-steps 1`.
The VQE benchmark is shot based: each Pauli term is measured with
`--vqe-shots` shots for each optimizer evaluation.
The output also includes simple correctness checks: Grover reports the measured
probability of the marked prefix, Shor reports the measured probability of
recovering order 4 for `N = 15`, and VQE reports the sampled energy difference
from Aer statevector.

### Example Timing Results

The following measurements were taken on one development environment with
`.venv_qdd/bin/python bench/compare_simulators.py`. They are intended as a
representative comparison, not as portable absolute performance numbers.

Grover prefix search used 20 qubits, 256 shots, marked prefix `101010101010`,
and the automatically selected Grover iteration count. The correctness metric
is the measured probability that the sampled bitstring starts with the marked
prefix.

| Simulator | Transpile [s] | Run [s] | Total [s] | Correctness |
| --- | ---: | ---: | ---: | ---: |
| Aer statevector | 0.080887 | 7.784454 | 7.865341 | prefix success 1.000 |
| Aer MPS | 1.398227 | 12.111329 | 13.509557 | prefix success 1.000 |
| QDD | 0.022978 | 0.109171 | 0.132149 | prefix success 1.000 |

This result matches the expected QDD behavior: the Grover oracle and diffusion
operator are highly structured, and the state remains compact enough for the
decision diagram representation to share nodes effectively.

Shor order finding used `N = 15`, base `a = 2`, 22 counting qubits, and 256
shots. The correctness metric is the probability that classical
post-processing recovers the expected order 4.

| Simulator | Transpile [s] | Run [s] | Total [s] | Correctness |
| --- | ---: | ---: | ---: | ---: |
| Aer statevector | 0.087232 | 9.524342 | 9.611574 | order-4 success 0.512 |
| Aer MPS | 0.064690 | 0.100180 | 0.164869 | order-4 success 0.508 |
| QDD | 0.039742 | 0.197238 | 0.236979 | order-4 success 0.469 |

This also broadly matches the expected picture for arithmetic circuits, but it
shows an important nuance: MPS can be very competitive when the particular
arithmetic circuit has low effective entanglement in the chosen qubit order.
QDD remains much faster than Aer statevector in this run, while MPS is fastest
for this small `N = 15` order-finding instance.

VQE used the H4 Hamiltonian from Qiskit Nature, `EfficientSU2(reps=8)`, one
optimizer evaluation, and 128 shots per Pauli term. The correctness metric is
reported as the sampled energy difference from Aer statevector, using the same
initial parameters and shot count.

| Simulator | Transpile [s] | Run [s] | Total [s] | Energy Delta vs Aer Statevector |
| --- | ---: | ---: | ---: | ---: |
| Aer statevector | 2.635268 | 1.146089 | 3.781358 | +0.000000 |
| Aer MPS | 2.710340 | 5.776522 | 8.486862 | +0.019080 |
| QDD | 2.356933 | 1.914973 | 4.271906 | -0.023605 |

This VQE result needs careful interpretation. It is consistent with the
guidance above when QDD is compared with dense statevector simulation: QDD is
slower than Aer statevector on this small-angle-rotation workload. However, QDD
is faster than Aer MPS in this specific run. That does not mean VQE is
generally favorable for QDD. The `EfficientSU2` ansatz uses full entanglement,
which is often unfavorable for MPS because MPS performance depends strongly on
low entanglement and a good one-dimensional qubit ordering. In contrast, this
particular VQE instance is still small enough that QDD's overhead does not
dominate as much as it would for larger or more random parameterized circuits.
For VQE-style workloads, treat Aer statevector as the primary baseline and
measure the specific ansatz before choosing QDD.

## References

- Qiskit Algorithms: Grover's algorithm examples,
  https://qiskit-community.github.io/qiskit-algorithms/tutorials/07_grover_examples.html
- Qiskit Nature: getting started with H2 VQE,
  https://qiskit-community.github.io/qiskit-nature/getting_started.html
- IBM Quantum Documentation: Shor's algorithm,
  https://qiskit.qotlabs.org/docs/tutorials/shors-algorithm
