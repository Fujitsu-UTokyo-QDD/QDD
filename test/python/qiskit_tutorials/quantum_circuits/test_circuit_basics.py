# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits/01_circuit_basics.ipynb  # noqa: E501

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


def test_circuit_basics():
    circ = QuantumCircuit(3)
    circ.h(0)
    circ.cx(0, 1)
    circ.cx(0, 2)

    meas = QuantumCircuit(3, 3)
    meas.barrier(range(3))
    meas.measure(range(3), range(3))
    qc = meas.compose(circ, range(3), front=True)

    backend = QddProvider().get_backend()
    qc_compiled = transpile(qc, backend, seed_transpiler=50)
    job_sim = backend.run(qc_compiled, shots=1024, seed_simulator=80)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(qc_compiled)
    print(counts)
    assert '000' in counts
    assert '111' in counts
    assert counts["000"] + counts["111"] == 1024
