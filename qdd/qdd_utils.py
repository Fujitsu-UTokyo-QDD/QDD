from qiskit.circuit import QuantumCircuit
from qiskit.qasm2 import dumps as qasm2_dumps, QASM2ExportError
from qiskit.quantum_info import SparsePauliOp


def _circuit_key(circuit: QuantumCircuit):
    """Convert a circuit to a hashable key.
    If a circuit is serializable, a string is returned.
    Otherwise, the id of the circuit is returned.
    """
    if circuit.num_parameters == 0 and "if_else" not in circuit.count_ops():
        try:
            return qasm2_dumps(circuit)
        except QASM2ExportError:
            return id(circuit)
    return id(circuit)


def _observable_key(observable: SparsePauliOp) -> tuple:
    """Convert an observable to a tuple for key of cache."""
    return tuple(sorted((pauli.to_label(), coeff) for pauli, coeff in zip(observable.paulis, observable.coeffs)))