import pytest
from qiskit.utils import QuantumInstance
from qiskit.algorithms import Shor
from qiskit.exceptions import QiskitError
from qdd import QddProvider

def test_shor_circuit():
    num=5
    backend = QddProvider().get_backend()
    
    quantum_instance = QuantumInstance(backend, shots=1024)
    shor = Shor(quantum_instance=quantum_instance)
    print(f'Actual number of qubits of circuit: {shor.construct_circuit(num).num_qubits}')
    with pytest.raises(QiskitError):
        shor.factor(num)