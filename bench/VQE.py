"""
See https://github.com/cda-tum/MQTBench/blob/main/mqt/bench/benchmarks/vqe.py
"""

from qiskit import Aer
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import SLSQP
from qiskit.circuit.library import RealAmplitudes
from qiskit.utils import QuantumInstance
import networkx as nx
import sys
import time
from mqt import ddsim

def get_examplary_max_cut_qp(n_nodes: int, degree: int = 2):
    """Returns a quadratic problem formulation of a max cut problem of a random graph.
    Keyword arguments:
    n_nodes -- number of graph nodes (and also number of qubits)
    degree -- edges per node
    """
    try:
        from qiskit_optimization.applications import Maxcut
    except Exception:

        print("Please install qiskit_optimization.")
        return None
    graph = nx.random_regular_graph(d=degree, n=n_nodes, seed=111)
    maxcut = Maxcut(graph)
    return maxcut.to_quadratic_program()

def vqe(num_qubits: int, back):
    qp = get_examplary_max_cut_qp(num_qubits)
    if back == 'Aer':
        backend = Aer.get_backend('qasm_simulator')
    elif back == 'ddsim':
        backend = ddsim.DDSIMProvider().get_backend('qasm_simulator')
    else:
        raise(NotImplementedError)
    sim = QuantumInstance(backend, shots=10, seed_simulator=10)

    ansatz = RealAmplitudes(num_qubits, reps=2)
    vqe = VQE(ansatz, optimizer=SLSQP(maxiter=25), quantum_instance=sim)
    start = time.time()
    vqe_result = vqe.compute_minimum_eigenvalue(qp.to_ising()[0])
    print('Time: {:.2f} sec'.format( time.time()-start))
    
if __name__=="__main__":
    vqe(int(sys.argv[1]), sys.argv[2])
