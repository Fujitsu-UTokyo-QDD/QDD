# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/simulators/1_aer_provider.ipynb  # noqa: E501

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
import pytest
from qiskit import Aer, QuantumCircuit, transpile
from qiskit.quantum_info import random_clifford, random_density_matrix, random_statevector, random_unitary

from qdd import QddProvider


def test_simulate_quantum_circuits():
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.cx(0, 1)
    circ.measure_all()

    simulator = QddProvider().get_backend()
    circ = transpile(circ, simulator, seed_transpiler=50)

    # get counts
    result = simulator.run(circ, seed_simulator=80).result()
    counts = result.get_counts(circ)
    assert '00' in counts
    assert '11' in counts

    # get memory
    result = simulator.run(circ, shots=100, memory=True, seed_simulator=80).result()
    memory = result.get_memory(circ)
    assert '00' in memory
    assert '11' in memory
    assert len(memory) == 100


def test_simulator_precision_options():
    simulator = QddProvider().get_backend()
    with pytest.raises(Exception):
        simulator.set_options(precision='single')


def test_custom_simulator_instructions():
    # load Aer class; this statement adds several attributes/methods to the QuantumCircuit class
    Aer.get_backend('aer_simulator')

    simulator = QddProvider().get_backend()

    # save_statevector
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.save_statevector()
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # save_unitary
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.save_unitary()
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # save multiple states
    steps = 5
    circ = QuantumCircuit(1)
    for i in range(steps):
        circ.save_statevector(label=f'psi_{i}')
        circ.rx(i * np.pi / steps, 0)
    circ.save_statevector(label=f'psi_{steps}')
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # save_state
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.save_state()
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # save_stabilizer
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.save_stabilizer()
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # save_density_matrix
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.save_density_matrix()
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # save_matrix_product_state
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.save_matrix_product_state()
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # save_superop
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.save_superop()
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # set_statevector
    num_qubits = 2
    psi = random_statevector(2 ** num_qubits, seed=100)
    circ = QuantumCircuit(num_qubits)
    circ.set_statevector(psi)
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # initialize
    circ = QuantumCircuit(num_qubits)
    circ.initialize([0, 1 / np.sqrt(2), 0, 1 / np.sqrt(2)], range(num_qubits))
    circ.measure_all()
    circ = transpile(circ, simulator, seed_transpiler=50)
    counts = simulator.run(circ, seed_simulator=80).result().get_counts()
    assert '01' in counts
    assert '11' in counts
    assert len(counts) == 2

    # set_density_matrix
    rho = random_density_matrix(2 ** num_qubits, seed=100)
    circ = QuantumCircuit(num_qubits)
    circ.set_density_matrix(rho)
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # set_stabilizer
    stab = random_clifford(num_qubits, seed=100)
    circ = QuantumCircuit(num_qubits)
    circ.set_stabilizer(stab)
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)

    # set_unitary
    unitary = random_unitary(2 ** num_qubits, seed=100)
    circ = QuantumCircuit(num_qubits)
    circ.set_unitary(unitary)
    with pytest.raises(Exception):
        transpile(circ, simulator, seed_transpiler=50)
