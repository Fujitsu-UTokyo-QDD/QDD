import pytest
from qdd import pyQDD
from qdd.qdd_backend import _qiskit_gates_1q_0param, _qiskit_gates_1q_1param, _qiskit_gates_1q_2param, _qiskit_gates_1q_3param, _qiskit_gates_1q_4param, _qiskit_gates_2q_0param, _qiskit_gates_2q_1param
import qiskit.circuit.library as qiskit_gates
import numpy as np
import math
import random
import scipy.stats


def test_1q_0param_gates():
    for qis,qdd in _qiskit_gates_1q_0param.items():
        if qis == qiskit_gates.MCXGate:
            continue
        qis_gate = qis()
        num_qubits = qis_gate.num_qubits
        qis_mat = qis_gate.to_matrix()
        controls = range(0,num_qubits-1)
        qdd_mat = pyQDD.makeGate(num_qubits, qdd, num_qubits-1, controls).getEigenMatrix()
        norm = np.linalg.norm(qdd_mat-qis_mat)
        assert(norm < 0.000001)

def test_1q_1param_gates():
    for qis,qdd in _qiskit_gates_1q_1param.items():
        if qis == qiskit_gates.MCPhaseGate:
            continue
        for _ in range(10):
            para = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para)
            num_qubits = qis_gate.num_qubits
            qis_mat = qis_gate.to_matrix()
            controls = range(0,num_qubits-1)
            qdd_mat = pyQDD.makeGate(num_qubits, qdd(para), num_qubits-1, controls).getEigenMatrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)

def test_1q_2param_gates():
    for qis,qdd in _qiskit_gates_1q_2param.items():
        for _ in range(20):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2)
            num_qubits = qis_gate.num_qubits
            qis_mat = qis_gate.to_matrix()
            controls = range(0,num_qubits-1)
            qdd_mat = pyQDD.makeGate(num_qubits, qdd(para1, para2), num_qubits-1, controls).getEigenMatrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)

def test_1q_3param_gates():
    for qis,qdd in _qiskit_gates_1q_3param.items():
        for _ in range(30):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            para3 = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2, para3)
            num_qubits = qis_gate.num_qubits
            qis_mat = qis_gate.to_matrix()
            controls = range(0,num_qubits-1)
            qdd_mat = pyQDD.makeGate(num_qubits, qdd(para1, para2, para3), num_qubits-1, controls).getEigenMatrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)

def test_1q_4param_gates():
    for qis,qdd in _qiskit_gates_1q_4param.items():
        for _ in range(40):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            para3 = random.uniform(-math.pi, math.pi)
            para4 = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2, para3, para4)
            num_qubits = qis_gate.num_qubits
            qis_mat = qis_gate.to_matrix()
            controls = range(0,num_qubits-1)
            qdd_mat = pyQDD.makeGate(num_qubits, qdd(para1, para2, para3, para4), num_qubits-1, controls).getEigenMatrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)

def test_2q_0param_gates():
    for qis,qdd in _qiskit_gates_2q_0param.items():
        qis_gate = qis()
        num_qubits = qis_gate.num_qubits
        qis_mat = qis_gate.to_matrix()
        controls = range(0,num_qubits-2)
        qdd_mat = pyQDD.makeTwoQubitGate(num_qubits, qdd(), num_qubits-1, num_qubits-2, controls).getEigenMatrix()
        norm = np.linalg.norm(qdd_mat-qis_mat)
        assert(norm < 0.000001)

def test_2q_1param_gates():
    for qis,qdd in _qiskit_gates_2q_1param.items():
        for _ in range(10):
            para = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para)
            num_qubits = qis_gate.num_qubits
            qis_mat = qis_gate.to_matrix()
            controls = range(0,num_qubits-2)
            qdd_mat = pyQDD.makeTwoQubitGate(num_qubits, qdd(para), num_qubits-1, num_qubits-2, controls).getEigenMatrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)

def test_unitary():
    random_matrix = scipy.stats.unitary_group.rvs(8)
    qis_mat = qiskit_gates.UnitaryGate(random_matrix).to_matrix()
    qdd_mat = pyQDD.unitary(random_matrix).getEigenMatrix()
    norm = np.linalg.norm(qdd_mat-qis_mat)
    assert(norm < 0.000001)
