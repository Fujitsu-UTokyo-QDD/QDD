from qdd import pyQDD
from qiskit.circuit.gate import Gate
from qiskit.circuit.exceptions import CircuitError


class QDDGate(Gate):

    def __init__(self, nQubits: int, medge, basis: str = "QDD"):
        super().__init__("qdd", nQubits, [medge])

    def validate_parameter(self, parameter):
        if not isinstance(parameter, pyQDD.mEdge):
            raise CircuitError(f"Invalid number of params for gate {self.name}.")
        return parameter