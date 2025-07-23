from qiskit.compiler import transpile
from qiskit import qasm2
import sys, time

from qdd import QddProvider, pyQDD

import QuantumFactoring
#from mpi4py import MPI
import sys

def run_qasm(filename, nNode=1, nThread=1, use_bcast=False, use_auto_swap=True, swap_ver='v1', use_ordering=False):

    #optimizer = QuantumCircuitOptimizer()

    print("backend creating...")
    backend = QddProvider().get_backend()
    if nNode>1:
        backend.set_options(use_mpi=True, use_auto_swap=use_auto_swap, swap_ver=swap_ver, use_bcast=use_bcast)
        print("MPI enabled, use_bcast=",use_bcast)
        print("use_auto_swap", use_auto_swap, swap_ver)
    pyQDD.set_gc_thr(1024*1024, 1024*1024) # You can manipulate these parameters
    #backend = Aer.get_backend('qasm_simulator')
    if use_ordering:
        backend.set_options(ordering=True)
        print("Ordering enabled")
    if nThread>1:
        backend.set_options(n_threads=nThread)

    qisc = qasm2.load(filename, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    print(qisc.num_qubits, len(qisc.data))

    new_circ = transpile(qisc, backend=backend, optimization_level=0)
    print(new_circ.num_qubits, len(new_circ.data))

    start = time.perf_counter()
    job = backend.run(new_circ, shots=1, show_progress = True)
    print(job.result().data()['nNodes'])
    end = time.perf_counter()
    print(end-start, " sec")


if __name__ == '__main__':
    # USAGE
    # mpiexec -n 1 python -u bench/shor_main.py 51 2 1 1 0 1 v1 # qasm_path $nNodes $nThreads(=1) $use_bcast $use_auto_swap $swap_ver
    args = sys.argv
    print(args)
    run_qasm(args[1], int(args[2]), int(args[3]), bool(int(args[4])), bool(int(args[5])), args[6], bool(int(args[7])))

    #MPI.Finalize()
