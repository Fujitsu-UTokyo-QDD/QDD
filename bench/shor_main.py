from qiskit.compiler import transpile
import sys, time

from qdd import QddProvider, pyQDD

import QuantumFactoring
from mpi4py import MPI
import sys

def compute_gates_for_all_ADD(N, a, nNode=1, nThread=1, use_bcast=False, use_auto_swap=True, swap_ver='v1'):

    #optimizer = QuantumCircuitOptimizer()

    print("##", N, len(bin(N)) - 2, a)

    print("backend creating...")
    backend = QddProvider().get_backend()
    if nNode>1:
        backend.set_options(use_mpi=True, use_auto_swap=use_auto_swap, swap_ver=swap_ver, use_bcast=use_bcast)
        print("MPI enabled, use_bcast=",use_bcast)
        print("use_auto_swap", use_auto_swap, swap_ver)
    pyQDD.set_gc_thr(1024*1024, 1024*1024) # You can manipulate these parameters
    #backend = Aer.get_backend('qasm_simulator')
    if nThread>1:
        backend.set_options(n_threads=nThread)

    f_gt_add = QuantumFactoring.QuantumCircuitForFactoringWithQulacs(N, "GT-ADD")
    f_gt_add.qisc.barrier(list(range(f_gt_add.qisc.num_qubits)))
    f_gt_add.find_order_algorithm(a, 100, measure=False)
    f_gt_add.qisc.measure_all()
    print(f_gt_add.qisc.num_qubits, len(f_gt_add.qisc.data))

    new_circ = transpile(f_gt_add.qisc, backend=backend, optimization_level=0)
    print(new_circ.num_qubits, len(new_circ.data))

    start = time.perf_counter()
    job = backend.run(new_circ, shots=1)
    print(len(job.result().data()['counts']))
    end = time.perf_counter()
    print(end-start, " sec")
    
    del f_gt_add


if __name__ == '__main__':
    # USAGE
    # mpiexec -n 1 python -u bench/shor_main.py 51 2 1 1 0 1 v1 # N a $nNodes $nThreads(=1) $use_bcast $use_auto_swap $swap_ver
    args = sys.argv
    compute_gates_for_all_ADD(int(args[1]), int(args[2]), int(args[3]), int(args[4]), bool(int(args[5])), bool(int(args[6])), args[7])

    MPI.Finalize()
