# CNF is generated automatically.
# cnfgen randkcnf 3 $i $(( $i*3 ))
# You need to install the lates tweedledum from source. Please see https://github.com/Qiskit/qiskit-terra/issues/6734

import sys
import time
from qiskit import Aer
from mqt import ddsim
from qiskit.utils import QuantumInstance
from qiskit.algorithms import Grover, AmplificationProblem
from qiskit.circuit.library.phase_oracle import PhaseOracle
from qiskit.exceptions import MissingOptionalLibraryError

def create_circuit(file_name):
    try:
        oracle = PhaseOracle.from_dimacs_file(file_name)
    except MissingOptionalLibraryError as ex:
        print(ex)
    problem = None
    if oracle is not None:
        problem = AmplificationProblem(oracle, is_good_state=oracle.evaluate_bitstring)
    return problem

def main(fname, back):
    if back == 'Aer':
        backend = Aer.get_backend('qasm_simulator')
    elif back == 'ddsim':
        backend = ddsim.DDSIMProvider().get_backend('qasm_simulator')
    else:
        raise(NotImplementedError)
    
    problem = create_circuit(fname)
    quantum_instance = QuantumInstance(backend, shots=1024)
    grover = Grover(quantum_instance=quantum_instance)
    result = None
    start = time.time()
    if problem is not None:
        result = grover.amplify(problem)
    print('Time: {:.2f} sec'.format( time.time()-start ))

if __name__ == "__main__":
    args = sys.argv
    main(args[1], args[2])
