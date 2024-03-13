# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits_advanced/02_operators_overview.ipynb  # noqa: E501

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
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import XGate
from qiskit.quantum_info import Operator, Pauli, process_fidelity

from qdd import QddProvider
from test.python.helpers.circuit_helper import get_counts


def test_using_operators_in_circuits():
    # Create an operator
    xx = Operator(Pauli('XX'))

    # Add to a circuit
    circ = QuantumCircuit(2, 2)
    circ.append(xx, [0, 1])
    circ.measure([0, 1], [0, 1])

    backend = QddProvider().get_backend()
    circ = transpile(circ, backend, basis_gates=['rx', 'ry', 'rz', 'cx'], seed_transpiler=50)
    job = backend.run(circ, seed_simulator=80)
    assert job.result().get_counts(0) == {'00':0, '01':0, '10':0,'11': 1024}

    # Add Pauli to a circuit directly
    circ2 = QuantumCircuit(2, 2)
    circ2.append(Pauli('XX'), [0, 1])
    circ2.measure([0, 1], [0, 1])
    backend = QddProvider().get_backend()
    _, counts = get_counts(circ2, backend, 1024)
    assert counts == {'00':0, '01':0, '10':0,'11': 1024}


def test_operator_methods():
    op_x = Operator(Pauli('X'))
    op_z = Operator(Pauli('Z'))
    assert op_x.tensor(op_z) == Operator(Pauli('XZ')) == op_z.expand(op_x)
    assert op_x.compose(op_z) == Operator([[0, 1], [-1, 0]])
    assert Operator(Pauli('X')) == Operator(XGate())
    assert Operator(XGate()) != np.exp(1j * 0.5) * Operator(XGate())

    # Two operators which differ only by phase
    op_a = Operator(XGate())
    op_b = np.exp(1j * 0.5) * Operator(XGate())

    # Compute process fidelity
    fidelity = process_fidelity(op_a, op_b)
    assert fidelity == pytest.approx(1.0)
