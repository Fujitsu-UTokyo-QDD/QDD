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
    with qc.if_test((cl_condition, 1)):
        qc.x(qr[1])
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
    with qc.if_test((cl_condition, 0)):
        qc.x(qr[1])
    qc.measure(qr[1], cl_result)
    job = backend.run(qc, shots=1024)
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
    with qc.if_test((cr, 10)):
        qc.x(qr[1])
    qc.measure(qr[1], cr[0])
    job = backend.run(qc, shots=1024)
    result = job.result()
    counts = result.get_counts()
    # The new if_test syntax seems to have a different classical bit ordering.
    # The measured result "1011" corresponds to cr[3]=1, cr[2]=0, cr[1]=1, cr[0]=1.
    # The classical register is 10, so cr[3]=1, cr[1]=0 -> this is wrong.
    # Let's check the bit ordering.
    # Qiskit's classical bits are ordered from right to left.
    # cr[0] is the rightmost bit.
    # measure(qr[0], cr[1]) -> cr[1] = 1
    # measure(qr[0], cr[3]) -> cr[3] = 1
    # So cr becomes '1010' in binary, which is 10 in decimal.
    # The condition (cr, 10) is met.
    # Then qc.x(qr[1]) is executed.
    # measure(qr[1], cr[0]) -> cr[0] = 1
    # The final state of cr should be '1011', which is 11.
    # Let's re-verify the old code's expected output. '1011' means
    # cr[3]=1, cr[2]=0, cr[1]=1, cr[0]=1.
    # This implies cr value was 13. But the measurement sets it to 10.
    # Ah, the count string is ordered cr[3]cr[2]cr[1]cr[0].
    # So '1011' is cr[3]=1, cr[2]=0, cr[1]=1, cr[0]=1.
    # Initial measurement sets cr[1]=1, cr[3]=1. So cr is "1010".
    # Condition is (cr, 10) which is true.
    # qc.x(qr[1]) is executed. qr[1] becomes 1.
    # qc.measure(qr[1], cr[0]) sets cr[0]=1.
    # So final cr is "1011". The original assertion seems correct.
    assert counts == {"1011": 1024}

    cr = ClassicalRegister(4)
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr, cr)
    qc.x(qr[0])
    qc.measure(qr[0], cr[1])
    qc.measure(qr[0], cr[3])  # cr = 10
    with qc.if_test((cr, 13)):  # this instruction is expected to be ignored
        qc.x(qr[1])
    qc.measure(qr[1], cr[0])
    job = backend.run(qc, shots=1024)
    result = job.result()
    counts = result.get_counts()
    assert counts == {"1010": 1024}
