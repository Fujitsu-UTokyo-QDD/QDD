# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits/1_getting_started_with_qiskit.ipynb  # noqa: E501

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

from qdd import QddProvider
from helpers.circuit_helper import assert_job_failed


def test_circuit_basics():
    circ = QuantumCircuit(3)
    circ.h(0)
    circ.cx(0, 1)
    circ.cx(0, 2)

    backend = QddProvider().get_backend()
    job = backend.run(circ, seed_simulator=80)
    # QddBackend does not support the evaluation of circuits with no measurements
    assert_job_failed(job)

    meas = QuantumCircuit(3, 3)
    meas.barrier(range(3))
    meas.measure(range(3), range(3))
    circ.add_register(meas.cregs[0])
    qc = circ.compose(meas)

    job = backend.run(transpile(qc, backend, seed_transpiler=50), shots=1024, seed_simulator=80)
    result = job.result()

    with pytest.raises(Exception):
        result.get_statevector(circ, decimals=3)
    with pytest.raises(Exception):
        result.get_unitary(circ, decimals=3)

    counts = result.get_counts(qc)
    assert '000' in counts
    assert '111' in counts
    assert len(counts) == 2
