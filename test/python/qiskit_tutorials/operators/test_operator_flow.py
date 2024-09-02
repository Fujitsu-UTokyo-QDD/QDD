# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/operators/01_operator_flow.ipynb  # noqa: E501

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
from qiskit.circuit import Parameter, QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import Estimator as QiskitEstimator
from qiskit.synthesis import SuzukiTrotter

from qdd.qdd_estimator import Estimator


def test_circuit_sampler():
    # Computes an expectation value through trotterlization

    two_qubit_h2 = SparsePauliOp.from_list(
        [
            ("II", -1.0523732),
            ("IZ", 0.39793742),
            ("ZI", -0.3979374),
            ("ZZ", -0.0112801),
            ("XX", 0.18093119),
        ]
    )

    evo_time = Parameter("θ")
    evo_gate = PauliEvolutionGate(
        operator=two_qubit_h2, time=evo_time, synthesis=SuzukiTrotter()
    )
    h2_measurement = two_qubit_h2.adjoint()
    bell = QuantumCircuit(2)
    bell.h(0)
    bell.cx(0, 1)
    evo_circ = bell.compose(evo_gate)

    estimator = Estimator()
    estimator_qiskit = QiskitEstimator()

    for time in range(8):
        qdd_value = (
            estimator.run(
                circuits=evo_circ, observables=h2_measurement, parameter_values=time
            )
            .result()
            .values
        )
        print(f"with QDD Estimator : {qdd_value}")
        qiskit_value = (
            estimator_qiskit.run(
                circuits=evo_circ, observables=h2_measurement, parameter_values=time
            )
            .result()
            .values
        )
        print(f"with Qiskit Estimator : {qiskit_value}")
        assert qiskit_value == pytest.approx(qdd_value, abs=0.02)
