# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/simulators/3_building_noise_models.ipynb  # noqa: E501
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/simulators/4_custom_gate_noise.ipynb  # noqa: E501

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

import pytest
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Kraus, Operator, SuperOp
from qiskit_aer import Aer
from qiskit_aer.noise import NoiseModel, QuantumError, pauli_error
from qiskit_ibm_runtime.fake_provider import FakeVigo

from qdd import QddBackend, QddProvider


def test_noise_simulation():
    device_backend = FakeVigo()
    with pytest.raises(Exception):
        QddBackend.from_backend(device_backend)

    p_gate1 = 0.05
    error_gate1 = pauli_error([('X', p_gate1), ('I', 1 - p_gate1)])
    noise_bit_flip = NoiseModel()
    noise_bit_flip.add_all_qubit_quantum_error(error_gate1, ["u1", "u2", "u3"])
    with pytest.raises(Exception):
        QddBackend(QddProvider(), noise_model=noise_bit_flip)


def test_gates_with_labels():
    """Tests the behavior of evaluating a circuit containing labeled gates."""

    # iSWAP matrix operator
    iswap_op = Operator([[1, 0, 0, 0],
                         [0, 0, 1j, 0],
                         [0, 1j, 0, 0],
                         [0, 0, 0, 1]])

    # CNOT in terms of iSWAP and single-qubit gates
    cx_circ = QuantumCircuit(2, name='cx<iSWAP>')

    # Add gates
    cx_circ.sdg(1)
    cx_circ.h(1)
    cx_circ.sdg(0)
    cx_circ.unitary(iswap_op, [0, 1], label='iswap')
    cx_circ.sdg(0)
    cx_circ.h(0)
    cx_circ.sdg(0)
    cx_circ.unitary(iswap_op, [0, 1], label='iswap')
    cx_circ.s(1)
    cx_circ.measure_all()

    qdd_backend = QddProvider().get_backend()
    circ = transpile(cx_circ, backend=qdd_backend, seed_transpiler=50)
    qdd_counts = qdd_backend.run(circ, seed_simulator=80).result().get_counts()

    aer_backend = Aer.get_backend('aer_simulator')
    circ = transpile(cx_circ, backend=aer_backend, seed_transpiler=50)
    aer_counts = aer_backend.run(circ, seed_simulator=80).result().get_counts()

    for key,value in qdd_counts:
        if key in aer_counts:
            assert value == aer_counts[key]


def test_quantum_error_and_channel_simulation():
    p_error = 0.05
    bit_flip = pauli_error([('X', p_error), ('I', 1 - p_error)])
    phase_flip = pauli_error([('Z', p_error), ('I', 1 - p_error)])
    bit_flip_kraus = Kraus(bit_flip)
    phase_flip_sop = SuperOp(phase_flip)
    superop_quantum_channel_inst = phase_flip_sop.to_instruction()

    backend = QddProvider().get_backend()

    # circuit with quantum channel-based instruction
    circ = QuantumCircuit(2)
    circ.append(superop_quantum_channel_inst, [0])
    circ.measure_all()
    with pytest.raises(Exception):
        transpile(circ, backend=backend, seed_transpiler=50)

    # circuit with quantum error-based instruction
    kraus_quantum_error_inst = QuantumError(bit_flip_kraus).to_instruction()
    circ = QuantumCircuit(2)
    circ.append(kraus_quantum_error_inst, [0])
    circ.measure_all()
    with pytest.raises(Exception):
        transpile(circ, backend=backend, seed_transpiler=50)
