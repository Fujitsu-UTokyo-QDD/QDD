#!/usr/bin/env python3
"""Compare QDD with Qiskit Aer statevector and MPS simulators.

The script builds representative Grover, Shor-style order-finding, and H2 VQE
circuits, then runs the same workload with:

* Qiskit Aer, method="statevector"
* Qiskit Aer, method="matrix_product_state"
* QDD

Examples:
    python bench/compare_simulators.py --benchmark grover --grover-qubits 26
    python bench/compare_simulators.py --benchmark shor --counting-qubits 22
    python bench/compare_simulators.py --benchmark vqe --molecule h4 --vqe-reps 8
    python bench/compare_simulators.py --benchmark all
"""

from __future__ import annotations

import argparse
import math
import os
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Iterable

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qdd-matplotlib")

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFTGate, UnitaryGate
from qiskit.circuit.library.grover_operator import grover_operator
from qiskit.quantum_info import Pauli
from qiskit_aer import AerSimulator
from scipy.optimize import minimize

from qdd import QddProvider, pyQDD


@dataclass
class Timing:
    simulator: str
    transpile_seconds: float
    run_seconds: float
    detail: str

    @property
    def total_seconds(self) -> float:
        return self.transpile_seconds + self.run_seconds


def format_timing_rows(title: str, rows: Iterable[Timing]) -> None:
    print(f"\n== {title} ==")
    print(
        f"{'simulator':<24} {'transpile[s]':>12} {'run[s]':>12} "
        f"{'total[s]':>12} detail"
    )
    for row in rows:
        print(
            f"{row.simulator:<24} {row.transpile_seconds:12.6f} "
            f"{row.run_seconds:12.6f} {row.total_seconds:12.6f} {row.detail}"
        )


def run_counts(
    name: str,
    backend,
    circuit: QuantumCircuit,
    *,
    shots: int,
    seed_simulator: int,
) -> tuple[Timing, dict[str, int]]:
    start = time.perf_counter()
    compiled = transpile(circuit, backend=backend, optimization_level=1)
    transpile_seconds = time.perf_counter() - start

    start = time.perf_counter()
    result = backend.run(compiled, shots=shots, seed_simulator=seed_simulator).result()
    run_seconds = time.perf_counter() - start

    counts = result.get_counts()
    top = max(counts.items(), key=lambda item: item[1])
    return (
        Timing(name, transpile_seconds, run_seconds, f"top={top[0]}:{top[1]}"),
        counts,
    )


def time_measured_batch_run(
    name: str,
    backend,
    circuits: list[QuantumCircuit],
    *,
    shots: int,
    seed_simulator: int,
) -> tuple[Timing, list[dict[str, int]]]:
    start = time.perf_counter()
    compiled = transpile(circuits, backend=backend, optimization_level=1)
    transpile_seconds = time.perf_counter() - start

    start = time.perf_counter()
    result = backend.run(compiled, shots=shots, seed_simulator=seed_simulator).result()
    run_seconds = time.perf_counter() - start

    return (
        Timing(name, transpile_seconds, run_seconds, f"shots={shots}"),
        [result.get_counts(index) for index in range(len(circuits))],
    )


def aer_backends() -> list[tuple[str, AerSimulator]]:
    return [
        ("Aer statevector", AerSimulator(method="statevector")),
        ("Aer MPS", AerSimulator(method="matrix_product_state")),
    ]


def qdd_qasm_backend():
    return "QDD", QddProvider().get_backend("qasm_simulator")


def build_prefix_phase_oracle(num_qubits: int, prefix: str) -> QuantumCircuit:
    """Build a phase oracle that marks bitstrings with a given leading prefix."""
    if len(prefix) > num_qubits:
        raise ValueError("prefix length must be <= num_qubits")

    oracle = QuantumCircuit(num_qubits, name=f"mark_prefix_{prefix}")
    prefix_qubits = list(range(num_qubits - len(prefix), num_qubits))
    for qubit, bit in zip(prefix_qubits, reversed(prefix)):
        if bit == "0":
            oracle.x(qubit)
    oracle.h(prefix_qubits[-1])
    oracle.mcx(prefix_qubits[:-1], prefix_qubits[-1])
    oracle.h(prefix_qubits[-1])
    for qubit, bit in zip(prefix_qubits, reversed(prefix)):
        if bit == "0":
            oracle.x(qubit)
    return oracle


def build_grover_bitstring(
    num_qubits: int,
    iterations: int | None = None,
    prefix: str = "101010101010",
) -> QuantumCircuit:
    """Build a scalable Grover circuit using Qiskit's Grover operator."""
    oracle = build_prefix_phase_oracle(num_qubits, prefix)
    grover_op = grover_operator(oracle)
    if iterations is None:
        iterations = max(1, round(math.pi / 4 * math.sqrt(2 ** len(prefix))))
    circuit = QuantumCircuit(num_qubits)
    circuit.h(range(num_qubits))
    for _ in range(iterations):
        circuit.compose(grover_op, inplace=True)
    circuit.measure_all()
    return circuit


def prefix_success_rate(counts: dict[str, int], prefix: str) -> float:
    shots = sum(counts.values())
    if shots == 0:
        return 0.0
    success = sum(count for bitstring, count in counts.items() if bitstring.startswith(prefix))
    return success / shots


def multiplication_mod_15_gate(a: int) -> UnitaryGate:
    """Return a permutation gate for |x> -> |a*x mod 15> on 4 qubits.

    This is the small N=15 modular-multiplication block commonly used in public
    Shor tutorials. Values outside 0..14 are mapped to themselves so the matrix
    remains a valid permutation over the full 4-qubit space.
    """
    if math.gcd(a, 15) != 1:
        raise ValueError("a must be coprime to 15")

    matrix = np.zeros((16, 16), dtype=complex)
    for x in range(16):
        y = (a * x) % 15 if x < 15 else x
        matrix[y, x] = 1.0
    return UnitaryGate(matrix, label=f"{a}x mod 15")


def build_shor_n15_order_finding(a: int = 2, counting_qubits: int = 8) -> QuantumCircuit:
    """Build the order-finding circuit for factoring 15 with base ``a``."""
    work_qubits = 4
    circuit = QuantumCircuit(counting_qubits + work_qubits, counting_qubits)
    counting = list(range(counting_qubits))
    work = list(range(counting_qubits, counting_qubits + work_qubits))

    circuit.h(counting)
    circuit.x(work[0])

    for control_index, control in enumerate(counting):
        multiplier = pow(a, 2**control_index, 15)
        gate = multiplication_mod_15_gate(multiplier).control()
        circuit.append(gate, [control, *work])

    inverse_qft = QFTGate(counting_qubits).inverse()
    circuit.append(inverse_qft, counting)
    circuit.measure(counting, counting)
    return circuit


def estimate_order_from_bitstring(bitstring: str, a: int, n: int) -> int | None:
    measured = bitstring.replace(" ", "")
    phase = int(measured, 2) / (2 ** len(measured))
    denominator = Fraction(phase).limit_denominator(n).denominator
    if denominator > 0 and pow(a, denominator, n) == 1:
        return denominator
    return None


def shor_order_success_rate(counts: dict[str, int], a: int, n: int, order: int) -> float:
    shots = sum(counts.values())
    if shots == 0:
        return 0.0
    success = sum(
        count
        for bitstring, count in counts.items()
        if estimate_order_from_bitstring(bitstring, a, n) == order
    )
    return success / shots


def build_vqe_problem(args) -> tuple[QuantumCircuit, object, str]:
    """Build a molecular VQE ansatz and Hamiltonian with Qiskit Nature."""
    from qiskit.circuit.library import efficient_su2
    from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper

    molecules = {
        "h2": ("H 0 0 0; H 0 0 0.735", "H2"),
        "h4": (
            "H 0 0 0; H 0 0 0.735; H 0 0 1.47; H 0 0 2.205",
            "H4",
        ),
        "lih": ("Li 0 0 0; H 0 0 1.6", "LiH"),
    }
    atom, molecule_label = molecules[args.molecule]
    driver = PySCFDriver(atom=atom, basis="sto-3g")
    problem = driver.run()
    mapper = JordanWignerMapper()
    hamiltonian = mapper.map(problem.hamiltonian.second_q_op())
    if args.vqe_ansatz == "uccsd":
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
        return ansatz, hamiltonian, f"Qiskit Nature {molecule_label}/UCCSD"

    ansatz = efficient_su2(
        hamiltonian.num_qubits,
        reps=args.vqe_reps,
        entanglement="full",
    )
    return (
        ansatz,
        hamiltonian,
        f"Qiskit Nature {molecule_label}/EfficientSU2(reps={args.vqe_reps})",
    )


def bind_parameters(circuit: QuantumCircuit, values: np.ndarray) -> QuantumCircuit:
    parameters = list(circuit.parameters)
    return circuit.assign_parameters(dict(zip(parameters, values)))


def measurement_circuit_for_pauli(
    ansatz: QuantumCircuit,
    pauli: Pauli,
) -> QuantumCircuit:
    """Return a measured copy of ``ansatz`` in the basis of ``pauli``."""
    circuit = ansatz.copy()
    for qubit, basis in enumerate(reversed(pauli.to_label())):
        if basis == "X":
            circuit.h(qubit)
        elif basis == "Y":
            circuit.sdg(qubit)
            circuit.h(qubit)
    circuit.measure_all()
    return circuit


def sampled_pauli_expectation(pauli: Pauli, counts: dict[str, int]) -> float:
    """Estimate a Pauli expectation value from computational-basis counts."""
    label = pauli.to_label()
    active_positions = [index for index, basis in enumerate(label) if basis != "I"]
    shots = sum(counts.values())
    if shots == 0:
        return 0.0

    total = 0
    for bitstring, count in counts.items():
        compact = bitstring.replace(" ", "")
        parity = 0
        for position in active_positions:
            if compact[position] == "1":
                parity ^= 1
        total += (-1 if parity else 1) * count
    return total / shots


def run_grover(args) -> None:
    circuit = build_grover_bitstring(
        num_qubits=args.grover_qubits,
        iterations=args.grover_iterations,
        prefix=args.grover_prefix,
    )
    rows = []
    for name, backend in [*aer_backends(), qdd_qasm_backend()]:
        timing, counts = run_counts(
            name,
            backend,
            circuit,
            shots=args.shots,
            seed_simulator=args.seed,
        )
        success_rate = prefix_success_rate(counts, args.grover_prefix)
        timing.detail += f", correctness=prefix_success:{success_rate:.3f}"
        rows.append(timing)
    format_timing_rows(
        f"Grover prefix search, {args.grover_qubits} qubits",
        rows,
    )


def run_shor(args) -> None:
    circuit = build_shor_n15_order_finding(a=args.shor_a, counting_qubits=args.counting_qubits)
    rows = []
    for name, backend in [*aer_backends(), qdd_qasm_backend()]:
        start = time.perf_counter()
        compiled = transpile(circuit, backend=backend, optimization_level=1)
        transpile_seconds = time.perf_counter() - start

        start = time.perf_counter()
        result = backend.run(compiled, shots=args.shots, seed_simulator=args.seed).result()
        run_seconds = time.perf_counter() - start

        counts = result.get_counts()
        top = max(counts.items(), key=lambda item: item[1])[0]
        success_rate = shor_order_success_rate(counts, args.shor_a, 15, order=4)
        detail = f"top={top}, correctness=order4_success:{success_rate:.3f}"
        rows.append(Timing(name, transpile_seconds, run_seconds, detail))

    format_timing_rows("Shor order finding, N=15", rows)


def run_vqe(args) -> None:
    ansatz, hamiltonian, source = build_vqe_problem(args)
    rng = np.random.default_rng(args.seed)
    x0 = rng.uniform(-0.05, 0.05, len(ansatz.parameters))

    eval_rows: list[Timing] = []
    summary_rows: list[Timing] = []

    pauli_terms = list(zip(hamiltonian.paulis, hamiltonian.coeffs))

    def run_backend_energy(name, backend) -> Callable[[np.ndarray], float]:
        calls = {"count": 0}

        def energy(parameters: np.ndarray) -> float:
            bound_ansatz = bind_parameters(ansatz, parameters)
            value = 0.0
            circuits = [
                measurement_circuit_for_pauli(bound_ansatz, pauli)
                for pauli, _ in pauli_terms
            ]
            timing, counts_list = time_measured_batch_run(
                name,
                backend,
                circuits,
                shots=args.vqe_shots,
                seed_simulator=args.seed,
            )
            for (pauli, coeff), counts in zip(pauli_terms, counts_list):
                value += float(np.real(coeff)) * sampled_pauli_expectation(pauli, counts)

            eval_rows.append(
                Timing(
                    f"{name} eval {calls['count']}",
                    timing.transpile_seconds,
                    timing.run_seconds,
                    f"energy={value:.12f}",
                )
            )
            calls["count"] += 1
            return value

        return energy

    print(f"\nVQE source: {source}")
    aer_statevector_backend, aer_mps_backend = aer_backends()
    backends = [
        aer_statevector_backend,
        aer_mps_backend,
        qdd_qasm_backend(),
    ]
    for name, backend in backends:
        backend_rows_before = len(eval_rows)
        objective = run_backend_energy(name, backend)
        start = time.perf_counter()
        result = minimize(
            objective,
            x0,
            method="COBYLA",
            options={"maxiter": args.vqe_steps, "rhobeg": 0.05},
        )
        optimizer_seconds = time.perf_counter() - start
        backend_eval_rows = eval_rows[backend_rows_before:]
        transpile_seconds = sum(row.transpile_seconds for row in backend_eval_rows)
        run_seconds = sum(row.run_seconds for row in backend_eval_rows)
        summary_rows.append(
            Timing(
                name,
                transpile_seconds,
                run_seconds,
                f"best_energy={result.fun:.12f}, optimizer_total={optimizer_seconds:.6f}s",
            )
        )

    reference_energy = summary_rows[0].detail.split("best_energy=")[1].split(",")[0]
    reference_energy_float = float(reference_energy)
    for row in summary_rows:
        energy = float(row.detail.split("best_energy=")[1].split(",")[0])
        row.detail += f", delta_vs_aer_sv={energy - reference_energy_float:+.6f}"

    format_timing_rows(
        f"{args.molecule.upper()} VQE, {len(ansatz.parameters)} params, "
        f"{args.vqe_steps} evaluations",
        summary_rows,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark",
        choices=["grover", "shor", "vqe", "all"],
        default="all",
        help="benchmark to run",
    )
    parser.add_argument("--shots", type=int, default=512, help="shots for sampled runs")
    parser.add_argument("--seed", type=int, default=1234, help="random simulator seed")
    parser.add_argument(
        "--grover-qubits",
        type=int,
        default=20,
        help="number of qubits for the Grover bitstring-search circuit",
    )
    parser.add_argument(
        "--grover-iterations",
        type=int,
        default=None,
        help="number of Grover iterations",
    )
    parser.add_argument(
        "--grover-prefix",
        default="101010101010",
        help="leading bit prefix marked by the Grover oracle",
    )
    parser.add_argument(
        "--counting-qubits",
        type=int,
        default=22,
        help="counting-register size for Shor order finding",
    )
    parser.add_argument("--shor-a", type=int, default=2, help="base for N=15 order finding")
    parser.add_argument(
        "--vqe-steps",
        type=int,
        default=2,
        help="maximum optimizer evaluations per simulator for VQE",
    )
    parser.add_argument(
        "--vqe-shots",
        type=int,
        default=256,
        help="shots per Pauli term and optimizer evaluation for VQE",
    )
    parser.add_argument(
        "--molecule",
        choices=["h2", "h4", "lih"],
        default="h4",
        help="molecule for the VQE benchmark",
    )
    parser.add_argument(
        "--vqe-ansatz",
        choices=["efficient-su2", "uccsd"],
        default="efficient-su2",
        help="ansatz for the H2 VQE benchmark",
    )
    parser.add_argument(
        "--vqe-reps",
        type=int,
        default=8,
        help="EfficientSU2 repetitions when --vqe-ansatz=efficient-su2",
    )
    parser.add_argument(
        "--qdd-gc-threshold",
        type=int,
        default=1 << 26,
        help="QDD garbage-collection threshold used during benchmarks",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pyQDD.set_gc_thr(args.qdd_gc_threshold, args.qdd_gc_threshold)
    if args.benchmark in {"grover", "all"}:
        run_grover(args)
    if args.benchmark in {"shor", "all"}:
        run_shor(args)
    if args.benchmark in {"vqe", "all"}:
        run_vqe(args)


if __name__ == "__main__":
    main()
