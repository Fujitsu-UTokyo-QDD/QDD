"""
See https://qiskit.org/documentation/locale/ja_JP/tutorials/algorithms/08_factorizers.html
"""

from qiskit import Aer
import sys
from mqt import ddsim
from qiskit.utils import QuantumInstance
from qiskit.algorithms import Shor
import time

def main(num,back):
    if back == 'Aer':
        backend = Aer.get_backend('qasm_simulator')
    elif back == 'ddsim':
        backend = ddsim.DDSIMProvider().get_backend('qasm_simulator')
    else:
        raise(NotImplementedError)
    quantum_instance = QuantumInstance(backend, shots=1024)
    shor = Shor(quantum_instance=quantum_instance)
    print(f'Actual number of qubits of circuit: {shor.construct_circuit(num).num_qubits}')
    start = time.time()
    shor.factor(num)
    print('Time: {:.2f} sec'.format( time.time()-start))

if __name__ == "__main__":
    args = sys.argv
    main(int(args[1]),args[2])
