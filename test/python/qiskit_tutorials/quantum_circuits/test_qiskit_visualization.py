# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits/2_plotting_data_in_qiskit.ipynb  # noqa: E501

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

from qiskit import QuantumCircuit, transpile

from qdd import QddProvider


def test_circuit_evaluation():
    bell = QuantumCircuit(2, 2)
    bell.h(0)
    bell.cx(0, 1)

    meas = QuantumCircuit(2, 2)
    meas.measure([0, 1], [0, 1])

    backend = QddProvider().get_backend()
    circ = bell.compose(meas)
    result = backend.run(transpile(circ, backend, seed_transpiler=50), shots=1000, seed_simulator=80).result()
    counts = result.get_counts(circ)

    assert '00' in counts
    assert '11' in counts
    assert len(counts) == 2
