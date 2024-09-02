"""Note: This module tests Qiskit gates each of which has 1-to-many mapping to a Qdd gates.
The following gates are tested in tests/test_conversion.py.
- Qiskit gates each of which has 1-to-1 mapping to a Qdd gate
- Instructions (reset, measure, and initialize)
- Conditional gates
"""

from qiskit import QuantumCircuit
from qiskit_aer import Aer

from qdd import QddProvider
from test.python.helpers.circuit_helper import (
    assert_probabilities_are_close,
    get_counts,
)


def test_gates_1_to_many_mapping():
    gates = [
        ("cy", [], [1, 3]),
        ("cz", [], [1, 3]),
        ("csx", [], [1, 3]),
        ("ch", [], [1, 3]),
        ("crz", [1], [1, 3]),
        ("cp", [1], [1, 3]),
        ("ccx", [], [1, 2, 5]),
        ("mcx", [], [[1, 2, 3], 5]),
        ("mcx", [], [[1, 2, 5], 4]),
        ("cswap", [], [1, 2, 5]),
    ]

    aer_backend = Aer.get_backend("aer_simulator")
    qdd_backend = QddProvider().get_backend()
    for qiskit_gate_name, gate_params, target_qubits in gates:
        circ = QuantumCircuit(6)
        circ.h(1)
        circ.h(2)
        qiskit_gate_method = getattr(circ, qiskit_gate_name)
        qiskit_gate_method(*gate_params, *target_qubits)
        circ.measure_all()
        _, aer_counts = get_counts(
            circ, aer_backend, n_shots=5000, optimization_level=0
        )
        _, qdd_counts = get_counts(
            circ, qdd_backend, n_shots=5000, optimization_level=0
        )
        print(f"Gate {qiskit_gate_name}:")
        print(aer_counts)
        print(qdd_counts)
        assert_probabilities_are_close(aer_counts, qdd_counts)


def test_ccx():
    aer_backend = Aer.get_backend("aer_simulator")
    qdd_backend = QddProvider().get_backend()
    circ = QuantumCircuit(3)
    circ.h(0)
    circ.h(1)
    circ.ccx(0, 1, 2)
    circ.measure_all()
    print(circ)
    _, aer_counts = get_counts(circ, aer_backend, n_shots=5000, optimization_level=0)
    _, qdd_counts = get_counts(circ, qdd_backend, n_shots=5000, optimization_level=0)
    print(aer_counts)
    print(qdd_counts)
    assert_probabilities_are_close(aer_counts, qdd_counts)
