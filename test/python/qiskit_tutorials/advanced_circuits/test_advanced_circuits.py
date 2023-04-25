# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits_advanced/01_advanced_circuits.ipynb  # noqa: E501

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import time
from itertools import combinations

import numpy as np
from qiskit import Aer, QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit import Gate, Parameter

from qdd import QddProvider
from test.python.helpers.circuit_helper import assert_probabilities_are_close, get_counts


def test_opaque_gates():
    my_gate = Gate(name='my_gate', num_qubits=2, params=[])
    qr = QuantumRegister(3, 'q')
    circ = QuantumCircuit(qr)
    circ.append(my_gate, [qr[0], qr[1]])
    circ.append(my_gate, [qr[1], qr[2]])
    print(circ)

    # The circ above cannot be executed via Both Aer simulator and Qdd simulator
    # Nothing to assert
    assert True


def test_composite_gates():
    # Build a sub-circuit
    sub_q = QuantumRegister(2)
    sub_circ = QuantumCircuit(sub_q, name='sub_circ')
    sub_circ.h(sub_q[0])
    sub_circ.crz(1, sub_q[0], sub_q[1])
    sub_circ.barrier()
    sub_circ.id(sub_q[1])
    sub_circ.u(1, 2, -2, sub_q[0])

    # Convert to a gate and stick it into an arbitrary place in the bigger circuit
    sub_inst = sub_circ.to_instruction()

    qr = QuantumRegister(3, 'q')
    circ = QuantumCircuit(qr)
    circ.h(qr[0])
    circ.cx(qr[0], qr[1])
    circ.cx(qr[1], qr[2])
    circ.append(sub_inst, [qr[1], qr[2]])
    circ.measure_all()

    decomposed_circ = circ.decompose()  # Does not modify original circuit
    print(decomposed_circ)

    aer_backend = Aer.get_backend('aer_simulator')
    qdd_backend = QddProvider().get_backend()

    _, aer_result = get_counts(circ, aer_backend, 5000, optimization_level=0)
    _, qdd_result = get_counts(circ, qdd_backend, 5000, optimization_level=0)

    assert_probabilities_are_close(aer_result, qdd_result)


def test_parameterized_circuits():
    theta = Parameter('θ')

    n = 5

    qc = QuantumCircuit(5, 1)

    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)

    qc.barrier()
    qc.rz(theta, range(5))
    qc.barrier()

    for i in reversed(range(n - 1)):
        qc.cx(i, i + 1)
    qc.h(0)
    qc.measure(0, 0)

    theta_range = np.linspace(0, 2 * np.pi, 128)

    circuits = [qc.bind_parameters({theta: theta_val})
                for theta_val in theta_range]

    aer_backend = Aer.get_backend('aer_simulator')
    qdd_backend = QddProvider().get_backend()

    aer_job = aer_backend.run(transpile(circuits, aer_backend, optimization_level=0, seed_transpiler=50),
                              shots=5000, seed_simulator=80)
    qdd_job = qdd_backend.run(transpile(circuits, qdd_backend, optimization_level=0, seed_transpiler=50),
                                    shots=5000, seed_simulator=80)

    aer_counts_list = aer_job.result().get_counts()
    qdd_counts_list = qdd_job.result().get_counts()

    for i, (aer_counts, qdd_counts) in enumerate(zip(aer_counts_list, qdd_counts_list)):
        assert_probabilities_are_close(aer_counts, qdd_counts)


def test_reducing_compilation_cost():
    start = time.time()
    qcs = []

    theta_range = np.linspace(0, 2 * np.pi, 32)

    for n in theta_range:
        qc = QuantumCircuit(5)

        for k in range(8):
            for i, j in combinations(range(5), 2):
                qc.cx(i, j)
            qc.rz(n, range(5))
            for i, j in combinations(range(5), 2):
                qc.cx(i, j)

        qcs.append(qc)

    backend = QddProvider().get_backend()
    _ = transpile(qcs, backend=backend, seed_transpiler=50)
    # QddBackend does not require assembling unlike Aer simulator
    # qobj = assemble(compiled_circuits, backend=backend)

    end = time.time()
    print('Time compiling over set of bound circuits: ', end - start)

    start = time.time()
    qc = QuantumCircuit(5)
    theta = Parameter('theta')

    for k in range(8):
        for i, j in combinations(range(5), 2):
            qc.cx(i, j)
        qc.rz(theta, range(5))
        for i, j in combinations(range(5), 2):
            qc.cx(i, j)

    transpiled_qc = transpile(qc, backend=backend, seed_transpiler=50)
    _ = [transpiled_qc.bind_parameters({theta: n}) for n in theta_range]
    # QddBackend does not require assembling unlike Aer simulator
    # qobj = assemble([transpiled_qc.bind_parameters({theta: n}) for n in theta_range], backend=backend)

    end = time.time()
    print('Time compiling over parameterized circuit, then binding: ', end - start)

    # In this test, we would like to see the above code finishes with no errors raised.
    # So, here, no specific properties to assert.
    assert True


def test_composition():
    phi = Parameter('phi')

    sub_circ1 = QuantumCircuit(2, name='sc_1')
    sub_circ1.rz(phi, 0)
    sub_circ1.rx(phi, 1)

    sub_circ2 = QuantumCircuit(2, name='sc_2')
    sub_circ2.rx(phi, 0)
    sub_circ2.rz(phi, 1)

    qc = QuantumCircuit(4)
    qr = qc.qregs[0]

    qc.append(sub_circ1.to_instruction(), [qr[0], qr[1]])
    qc.append(sub_circ2.to_instruction(), [qr[0], qr[1]])
    qc.append(sub_circ2.to_instruction(), [qr[2], qr[3]])

    qc = qc.bind_parameters({phi: 1})
    qc.measure_all()

    aer_backend = Aer.get_backend('aer_simulator')
    qdd_backend = QddProvider().get_backend()
    _, aer_result = get_counts(qc, aer_backend, 5000, optimization_level=0)
    _, qdd_result = get_counts(qc, qdd_backend, 5000, optimization_level=0)

    assert_probabilities_are_close(aer_result, qdd_result)


def test_composition_with_parameter_map():
    p = Parameter('p')
    qc = QuantumCircuit(3, name='oracle')
    qc.rz(p, 0)
    qc.cx(0, 1)
    qc.rz(p, 1)
    qc.cx(1, 2)
    qc.rz(p, 2)

    theta = Parameter('theta')
    phi = Parameter('phi')
    gamma = Parameter('gamma')

    qr = QuantumRegister(9)
    larger_qc = QuantumCircuit(qr)
    larger_qc.append(qc.to_instruction({p: theta}), qr[0:3])
    larger_qc.append(qc.to_instruction({p: phi}), qr[3:6])
    larger_qc.append(qc.to_instruction({p: gamma}), qr[6:9])

    larger_qc = larger_qc.bind_parameters({theta: 2, phi: 4, gamma: 6})
    larger_qc.measure_all()

    aer_backend = Aer.get_backend('aer_simulator')
    qdd_backend = QddProvider().get_backend()
    _, aer_result = get_counts(larger_qc, aer_backend, 5000, optimization_level=0)
    _, qdd_result = get_counts(larger_qc, qdd_backend, 5000, optimization_level=0)

    assert_probabilities_are_close(aer_result, qdd_result)
