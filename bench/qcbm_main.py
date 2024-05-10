from qiskit import QuantumCircuit
from qdd import QddProvider
import numpy as np
import time
import sys

def exec(backend, circuit):
    start = time.perf_counter()
    print("Start", circuit.num_qubits)
    job = backend.run(circuit, shots=1, optimization_level=0)
    res = job.result()
    end = time.perf_counter()
    print(end-start, "sec")

def native_execute(circuit, use_bcast=False, use_auto_swap=True, swap_ver='v1'):
    backend=QddProvider().get_backend()
    backend.set_options(use_mpi=True, use_auto_swap=use_auto_swap, swap_ver=swap_ver, use_bcast=use_bcast)
    print("use_auto_swap",use_auto_swap)
    #experiment = transpile(circuit, backend, optimization_level=0)
    exec(backend, circuit)

def run_bench(benchmark, nqubits, gate, args=(3, )):
    qc = QuantumCircuit(nqubits)
    getattr(qc, gate)(*args)
    native_execute(benchmark, qc)

def first_rotation(circuit, nqubits):
    circuit.rx(np.random.rand(), range(nqubits))
    circuit.rz(np.random.rand(), range(nqubits))
    return circuit


def mid_rotation(circuit, nqubits):
    circuit.rz(np.random.rand(), range(nqubits))
    circuit.rx(np.random.rand(), range(nqubits))
    circuit.rz(np.random.rand(), range(nqubits))
    return circuit


def last_rotation(circuit, nqubits):
    circuit.rz(np.random.rand(), range(nqubits))
    circuit.rx(np.random.rand(), range(nqubits))
    return circuit


def entangler(circuit, pairs):
    for a, b in pairs:
        circuit.cx(a, b)
    return circuit

def generate_qcbm_circuit(nqubits, depth, pairs):
    circuit = QuantumCircuit(nqubits)
    first_rotation(circuit, nqubits)
    entangler(circuit, pairs)
    for k in range(depth - 1):
        mid_rotation(circuit, nqubits)
        entangler(circuit, pairs)
    last_rotation(circuit, nqubits)
    return circuit

args = sys.argv
nqubit = 15
pairs = [(i, (i + 1) % nqubit) for i in range(nqubit)]
circuit = generate_qcbm_circuit(nqubit, 9, pairs)
circuit.measure_all()
native_execute(circuit, bool(int(args[1])), bool(int(args[2])), args[3])

# USAGE
# mpiexec -n 1 ./qcbm_main.py 0 1 v1 # $use_bcast $use_auto_swap $swap_ver 
