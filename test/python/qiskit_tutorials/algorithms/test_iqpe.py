# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/algorithms/09_IQPE.ipynb  # noqa: E501

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

import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile

from qdd import QddProvider
from qdd.qdd_sampler import Sampler


def x_measurement(qc, qubit, cbit):
    """Measure 'qubit' in the X-basis, and store the result in 'cbit'"""
    qc.h(qubit)
    qc.measure(qubit, cbit)


def test_iqpe_1qubit():
    nq = 2
    m = 2
    q = QuantumRegister(nq, 'q')
    c = ClassicalRegister(m, 'c')

    qc_s = QuantumCircuit(q, c)
    qc_s.h(0)
    qc_s.x(1)
    for _ in range(2 ** (m - 1)):
        qc_s.cp(np.pi / 2, 0, 1)

    x_measurement(qc_s, q[0], c[0])

    qc_s.reset(0)
    qc_s.h(0)
    qc_s.p(-np.pi / 2, 0).c_if(c, 1)
    for _ in range(2 ** (m - 2)):
        qc_s.cp(np.pi / 2, 0, 1)

    x_measurement(qc_s, q[0], c[1])

    backend = QddProvider().get_backend()
    n_shots = 1024
    count0 = backend.run(transpile(qc_s, backend, seed_transpiler=50), shots=n_shots, seed_simulator=80)\
        .result().get_counts()

    key_new = [str(int(key, 2) / 2 ** m) for key in list(count0.keys())]
    count1 = dict(zip(key_new, count0.values()))
    print(count1)
    assert count1["0.25"] == n_shots


def test_iqpe_2qubits():
    nq = 3  # number of qubits
    m = 3  # number of classical bits
    q = QuantumRegister(nq, 'q')
    c = ClassicalRegister(m, 'c')

    qc = QuantumCircuit(q, c)
    qc.h(0)
    qc.x([1, 2])
    for _ in range(2 ** (m - 1)):
        qc.mcp(np.pi / 4, [0, 1], 2)

    x_measurement(qc, q[0], c[0])
    qc.reset(0)
    qc.h(0)
    qc.p(-np.pi / 2, 0).c_if(c, 1)
    for _ in range(2 ** (m - 2)):
        qc.mcp(np.pi / 4, [0, 1], 2)

    x_measurement(qc, q[0], c[1])

    # initialization of qubit q0
    qc.reset(0)
    qc.h(0)

    # phase correction
    qc.p(-np.pi / 4, 0).c_if(c, 1)

    qc.p(-np.pi / 2, 0).c_if(c, 2)
    qc.p(-3 * np.pi / 2, 0).c_if(c, 3)

    # c-U operations
    for _ in range(2 ** (m - 3)):
        qc.mcp(np.pi / 4, [0, 1], 2)

    # X measurement
    qc.h(0)
    qc.measure(0, 2)

    backend = QddProvider().get_backend()
    n_shots = 1024
    count0 = backend.run(transpile(qc, backend, seed_transpiler=50), shots=n_shots, seed_simulator=80)\
        .result().get_counts()

    key_new = [str(int(key, 2) / 2 ** m) for key in list(count0.keys())]
    count1 = dict(zip(key_new, count0.values()))
    print(count1)

    assert count1["0.125"] == n_shots

def test_iqpe_1qubit_with_sampler():
    nq = 2
    m = 2
    q = QuantumRegister(nq, 'q')
    c = ClassicalRegister(m, 'c')

    qc_s = QuantumCircuit(q, c)
    qc_s.h(0)
    qc_s.x(1)
    for _ in range(2 ** (m - 1)):
        qc_s.cp(np.pi / 2, 0, 1)

    x_measurement(qc_s, q[0], c[0])

    qc_s.reset(0)
    qc_s.h(0)
    qc_s.p(-np.pi / 2, 0).c_if(c, 1)
    for _ in range(2 ** (m - 2)):
        qc_s.cp(np.pi / 2, 0, 1)

    x_measurement(qc_s, q[0], c[1])

    sampler = Sampler()
    job = sampler.run(qc_s)
    result = job.result()
    dist0 = result.quasi_dists[0]

    key_new = [str(key/2 ** m) for key in list(dist0.keys())]
    dist1 = dict(zip(key_new, dist0.values()))
    print(dist1)
    assert dist1["0.25"] == 1.0

    # test below is for the case where the quantum circuit has reset gate and the sampler caluculates the exact distribution
    sampler = Sampler(run_options={"shots": None})
    job = sampler.run(qc_s)
    result = job.result()
    dist0 = result.quasi_dists[0]

    key_new = [str(key/2 ** m) for key in list(dist0.keys())]
    dist1 = dict(zip(key_new, dist0.values()))
    print(dist1)
    assert dist1["0.25"] == 1.0


def test_iqpe_2qubits_with_sampler():
    nq = 3  # number of qubits
    m = 3  # number of classical bits
    q = QuantumRegister(nq, 'q')
    c = ClassicalRegister(m, 'c')

    qc = QuantumCircuit(q, c)
    qc.h(0)
    qc.x([1, 2])
    for _ in range(2 ** (m - 1)):
        qc.mcp(np.pi / 4, [0, 1], 2)

    x_measurement(qc, q[0], c[0])
    qc.reset(0)
    qc.h(0)
    qc.p(-np.pi / 2, 0).c_if(c, 1)
    for _ in range(2 ** (m - 2)):
        qc.mcp(np.pi / 4, [0, 1], 2)

    x_measurement(qc, q[0], c[1])

    # initialization of qubit q0
    qc.reset(0)
    qc.h(0)

    # phase correction
    qc.p(-np.pi / 4, 0).c_if(c, 1)

    qc.p(-np.pi / 2, 0).c_if(c, 2)
    qc.p(-3 * np.pi / 2, 0).c_if(c, 3)

    # c-U operations
    for _ in range(2 ** (m - 3)):
        qc.mcp(np.pi / 4, [0, 1], 2)

    # X measurement
    qc.h(0)
    qc.measure(0, 2)

    sampler = Sampler()
    job = sampler.run(qc)
    result = job.result()
    dist0 = result.quasi_dists[0]

    key_new = [str(key/2 ** m) for key in list(dist0.keys())]
    dist1 = dict(zip(key_new, dist0.values()))
    print(dist1)
    assert dist1["0.125"] == 1.0

    # test below is for the case where the quantum circuit has reset gate and the sampler caluculates the exact distribution
    sampler = Sampler(run_options={"shots": None})
    job = sampler.run(qc)
    result = job.result()
    dist0 = result.quasi_dists[0]

    key_new = [str(key/2 ** m) for key in list(dist0.keys())]
    dist1 = dict(zip(key_new, dist0.values()))
    print(dist1)
    assert dist1["0.125"] == 1.0
