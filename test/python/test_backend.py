from qiskit import QuantumCircuit, execute
from qdd import QddProvider

def test_simple_circuit():
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.cx(1,0)

    backend = QddProvider().get_backend()
    job = execute(circ, backend=backend, shots=20)