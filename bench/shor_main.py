from qiskit.compiler import transpile
import sys, time

from qdd import QddProvider, pyQDD

import QuantumFactoring
from mpi4py import MPI
import sys

def compute_gates_for_all_ADD(N, a, nNode=1):

    #optimizer = QuantumCircuitOptimizer()

    print("##", N, len(bin(N)) - 2, a)

    print("backend creating...")
    backend = QddProvider().get_backend()
    if nNode>1:
        backend.set_options(use_mpi=True)
        print("MPI enabled")
    pyQDD.set_gc_thr(1024*1024, 1024*1024) # You can manipulate these parameters
    #backend = Aer.get_backend('qasm_simulator')

    f_gt_add = QuantumFactoring.QuantumCircuitForFactoringWithQulacs(N, "GT-ADD")
    f_gt_add.qisc.barrier(list(range(f_gt_add.qisc.num_qubits)))
    f_gt_add.find_order_algorithm(a, 100, measure=False)
    f_gt_add.qisc.measure_all()
    print(f_gt_add.qisc.num_qubits, len(f_gt_add.qisc.data))

    new_circ = transpile(f_gt_add.qisc, backend=backend, optimization_level=0)
    print(new_circ.num_qubits, len(new_circ.data))

    start = time.perf_counter()
    job = backend.run(new_circ, shots=1, show_progress = True)
    print(len(job.result().data()['counts']))
    end = time.perf_counter()
    print(end-start, " sec")
    
    del f_gt_add


if __name__ == '__main__':
    # USAGE
    # mpiexec -n 1 python -u bench/shor_main.py 51 2 1 # N a $nNodes 
    args = sys.argv
    compute_gates_for_all_ADD(int(args[1]), int(args[2]), int(args[3]))

    MPI.Finalize()
