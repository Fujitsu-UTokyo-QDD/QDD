from qiskit.providers import Provider
from qiskit.transpiler import Target

from qdd.qdd_backend import QddBackend


class QddSVBackend(QddBackend):
    def __init__(self, provider: Provider):
        self._save_SV = True
        self._NUM_QUBITS = 20
        super().__init__(provider, name="statevector_simulator")
