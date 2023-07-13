import pytest
import pyQDD
import qiskit.circuit.library.standard_gates as qiskit_gates
from qiskit import QuantumCircuit as QiskitCircuit, QuantumRegister
from qiskit.quantum_info import Operator
import numpy as np
import math
import random

def test_1q_gates():
    _qiskit_gates_1q = {
        qiskit_gates.HGate: "H",
        qiskit_gates.IGate: "I",
        qiskit_gates.SdgGate: "Sdag",
        qiskit_gates.SGate: "S",
        qiskit_gates.SXdgGate: "SXdag",
        qiskit_gates.SXGate: "SX",
        qiskit_gates.TdgGate: "Tdag",
        qiskit_gates.TGate: "T",
        qiskit_gates.XGate: "X",
        qiskit_gates.YGate: "Y",
        qiskit_gates.ZGate: "Z"
    }

    for qis,str in _qiskit_gates_1q.items():
        qdd_mat = pyQDD.makeGate(1, str, 0).getEigenMatrix()
        qis_mat = qis().to_matrix()
        norm = np.linalg.norm(qdd_mat-qis_mat)
        assert(norm < 0.000001)

def test_2q_gates():
    _qiskit_gates_2q = {
        qiskit_gates.CXGate: pyQDD.CX,
        qiskit_gates.SwapGate: pyQDD.SWAP
    }

    for qis,qdd in _qiskit_gates_2q.items():
        qdd_mat = qdd(2,1,0).getEigenMatrix()
        qis_mat = qis().to_matrix()
        norm = np.linalg.norm(qdd_mat-qis_mat)
        assert(norm < 0.000001)

def test_1q_rotation_1():
    _qiskit_rotations_1q = {
        qiskit_gates.RXGate: pyQDD.rxmat,
        qiskit_gates.RYGate: pyQDD.rymat,
        qiskit_gates.RZGate: pyQDD.rzmat,
        qiskit_gates.U1Gate: pyQDD.u1,
        qiskit_gates.PhaseGate: pyQDD.p,
    }

    for qis,qdd in _qiskit_rotations_1q.items():
        for _ in range(10):
            para = random.uniform(-math.pi, math.pi)
            qdd_mat = pyQDD.makeGate(1, qdd(para), 0).getEigenMatrix()
            qis_mat = qis(para).to_matrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)

    _qiskit_rotations_1q_2 = {
        qiskit_gates.U2Gate: pyQDD.u2,
        qiskit_gates.RGate: pyQDD.r,
    }

    for qis,qdd in _qiskit_rotations_1q_2.items():
        for _ in range(20):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            qdd_mat = pyQDD.makeGate(1, qdd(para1, para2), 0).getEigenMatrix()
            qis_mat = qis(para1, para2).to_matrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)


    _qiskit_rotations_1q_3 = {
        qiskit_gates.U3Gate: pyQDD.u3,
        qiskit_gates.UGate: pyQDD.u,
    }

    for qis,qdd in _qiskit_rotations_1q_3.items():
        for _ in range(30):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            para3 = random.uniform(-math.pi, math.pi)
            qdd_mat = pyQDD.makeGate(1, qdd(para1, para2, para3), 0).getEigenMatrix()
            qis_mat = qis(para1, para2, para3).to_matrix()
            norm = np.linalg.norm(qdd_mat-qis_mat)
            assert(norm < 0.000001)

def test_1q_control():
    _qiskit_1q_control = {
        qiskit_gates.CYGate: "Y",
        qiskit_gates.CZGate: "Z",
        qiskit_gates.CSXGate: "SX",
    }

    for qis,qdd in _qiskit_1q_control.items():
        reg = QuantumRegister(2)
        circ = QiskitCircuit(reg)
        circ.append(qis(),reg)
        qis_mat = Operator(circ).data

        gate = pyQDD.makeControlGate(2, qdd, 1, [0])
        qdd_mat = gate.getEigenMatrix()
        norm = np.linalg.norm(qdd_mat-qis_mat)
        assert(norm < 0.000001)

