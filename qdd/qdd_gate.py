from . import pyQDD
from qiskit.circuit.gate import Gate
from qiskit.circuit.exceptions import CircuitError


class QDDGate(Gate):
    """Qiskit gate wrapper around a QDD matrix edge."""

    def __init__(self, nQubits: int, medge, basis: str = "QDD"):
        """Create a QDD-backed gate.

        Args:
            nQubits: Number of qubits the gate acts on.
            medge: QDD matrix edge representing the operation.
            basis: Reserved basis label for compatibility.
        """
        super().__init__("qdd", nQubits, [medge])

    def validate_parameter(self, parameter):
        """Validate that the gate parameter is a QDD matrix edge."""
        if not isinstance(parameter, pyQDD.mEdge):
            raise CircuitError(f"Invalid number of params for gate {self.name}.")
        return parameter
