from qdd.qdd_backend import (
    _qiskit_gates_1q_0param,
    _qiskit_gates_1q_1param,
    _qiskit_gates_1q_2param,
    _qiskit_gates_1q_3param,
    _qiskit_gates_1q_4param,
    _qiskit_gates_2q_0param,
    _qiskit_gates_2q_1param,
)
from qdd.qdd_sampler import Sampler as qdd_sampler
from qiskit import QuantumCircuit
import qiskit.circuit.library as qiskit_gates
from qiskit.quantum_info import Statevector
import numpy as np
import math
import random
import scipy.stats

sampler_qdd = qdd_sampler()
# MPI backend broadcasts rank 0's circuit; every rank must build the same oracle circuit.
rng = random.Random(1234)
unitary_rng = np.random.RandomState(1234)

qc_size = 6


def assert_qdd_matches_statevector(qc):
    qdd_dist = (
        sampler_qdd.run(pubs=[(qc, [])], is_exact=True).result()[0].data.quasi_dist
    )
    statevector_circuit = qc.remove_final_measurements(inplace=False)
    statevector_probs = Statevector.from_instruction(statevector_circuit).probabilities()
    qdd_probs = np.array(
        [qdd_dist.get(i, 0.0) for i in range(2**statevector_circuit.num_qubits)]
    )

    assert np.allclose(qdd_probs, statevector_probs, rtol=0.0, atol=1e-12)


def test_1q_0param_gates():
    for qis, _ in _qiskit_gates_1q_0param.items():
        if qis == qiskit_gates.MCXGate:
            continue
        qis_gate = qis()
        num_qubits = qis_gate.num_qubits
        for _ in range(3):
            targets = rng.sample(range(qc_size), num_qubits)
            qc = QuantumCircuit(qc_size)
            qc.h(range(qc_size))
            qc.append(qis_gate, targets)
            qc.measure_all()
            assert_qdd_matches_statevector(qc)


def test_1q_1param_gates():
    for qis, _ in _qiskit_gates_1q_1param.items():
        if qis == qiskit_gates.MCPhaseGate:
            continue
        for _ in range(10):
            para = rng.uniform(-math.pi, math.pi)
            qis_gate = qis(para)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = rng.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                assert_qdd_matches_statevector(qc)


def test_1q_2param_gates():
    for qis, _ in _qiskit_gates_1q_2param.items():
        for _ in range(20):
            para1 = rng.uniform(-math.pi, math.pi)
            para2 = rng.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = rng.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                assert_qdd_matches_statevector(qc)


def test_1q_3param_gates():
    for qis, _ in _qiskit_gates_1q_3param.items():
        for _ in range(30):
            para1 = rng.uniform(-math.pi, math.pi)
            para2 = rng.uniform(-math.pi, math.pi)
            para3 = rng.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2, para3)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = rng.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                assert_qdd_matches_statevector(qc)


def test_1q_4param_gates():
    for qis, _ in _qiskit_gates_1q_4param.items():
        for _ in range(40):
            para1 = rng.uniform(-math.pi, math.pi)
            para2 = rng.uniform(-math.pi, math.pi)
            para3 = rng.uniform(-math.pi, math.pi)
            para4 = rng.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2, para3, para4)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = rng.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                assert_qdd_matches_statevector(qc)


def test_2q_0param_gates():
    for qis, _ in _qiskit_gates_2q_0param.items():
        qis_gate = qis()
        num_qubits = qis_gate.num_qubits
        for _ in range(3):
            targets = rng.sample(range(qc_size), num_qubits)
            qc = QuantumCircuit(qc_size)
            qc.h(range(qc_size))
            qc.append(qis_gate, targets)
            qc.measure_all()
            assert_qdd_matches_statevector(qc)


def test_2q_1param_gates():
    for qis, _ in _qiskit_gates_2q_1param.items():
        for _ in range(10):
            para = rng.uniform(-math.pi, math.pi)
            qis_gate = qis(para)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = rng.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                assert_qdd_matches_statevector(qc)


def test_unitary():
    bit = 3
    for _ in range(10):
        random_matrix = scipy.stats.unitary_group.rvs(2**bit, random_state=unitary_rng)
        for _ in range(3):
            targets = rng.sample(range(qc_size), bit)
            qc = QuantumCircuit(qc_size)
            qc.h(range(qc_size))
            qc.append(qiskit_gates.UnitaryGate(random_matrix), targets)
            qc.measure_all()
            assert_qdd_matches_statevector(qc)
