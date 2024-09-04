import pytest
from qiskit.circuit import (
    QuantumCircuit,
    ClassicalRegister,
    QuantumRegister,
    Clbit,
)
from qiskit.circuit.library import XGate

from qdd import QddProvider

backend = QddProvider().get_backend()


def test_classical_bit_control():
    cl_condition = Clbit()
    cl_result = Clbit()
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr, [cl_condition, cl_result])
    qc.x(qr[0])
    qc.measure(qr[0], cl_condition)
    x_gate = XGate().c_if(cl_condition, 1)
    qc.append(x_gate, [qr[1]])
    qc.measure(qr[1], cl_result)
    job = backend.run(qc, shots=1024)
    result = job.result()
    counts = result.get_counts()
    assert counts == {"11": 1024}

    cl_condition = Clbit()
    cl_result = Clbit()
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr, [cl_condition, cl_result])
    qc.x(qr[0])
    qc.measure(qr[0], cl_condition)  # cl_condition = 1
    x_gate = XGate().c_if(cl_condition, 0)  # this instruction is expected to be ignored
    circ = qc
    circ.append(x_gate, [qr[1]])
    circ.measure(qr[1], cl_result)
    job = backend.run(circ, shots=1024)
    result = job.result()
    counts = result.get_counts()
    assert counts == {"01": 1024}


def test_classical_register_control():
    cr = ClassicalRegister(4)
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr, cr)
    qc.x(qr[0])
    qc.measure(qr[0], cr[1])
    qc.measure(qr[0], cr[3])
    x_gate = XGate().c_if(cr, 10)
    qc.append(x_gate, [qr[1]])
    qc.measure(qr[1], cr[0])
    job = backend.run(qc, shots=1024)
    result = job.result()
    counts = result.get_counts()
    assert counts == {"1011": 1024}

    cr = ClassicalRegister(4)
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr, cr)
    qc.x(qr[0])
    qc.measure(qr[0], cr[1])
    qc.measure(qr[0], cr[3])  # cr = 10
    x_gate = XGate().c_if(cr, 13)  # this instruction is expected to be ignored
    qc.append(x_gate, [qr[1]])
    qc.measure(qr[1], cr[0])
    job = backend.run(qc, shots=1024)
    result = job.result()
    counts = result.get_counts()
    assert counts == {"1010": 1024}
