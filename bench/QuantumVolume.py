from qiskit.circuit.library import QuantumVolume
import qiskit
import sys
from mqt import ddsim
from qiskit import Aer

def main(num, back):
    if back == 'Aer':
        backend = Aer.get_backend('qasm_simulator')
    elif back == 'ddsim':
        backend = ddsim.DDSIMProvider().get_backend('qasm_simulator')
    else:
        raise(NotImplementedError)
    c = QuantumVolume(num,10)
    job = qiskit.execute(c, backend, shots=1)
    print('Time: {:.2f} sec'.format( job.result().time_taken ))
    return

if __name__ == "__main__":
    args = sys.argv
    main(int(args[1]),args[2])
