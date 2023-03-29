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
from qiskit.circuit import Parameter
from qiskit.opflow import CX, CircuitSampler, H, I, PauliExpectation, PauliTrotterEvolution, StateFn, Suzuki, X, Z, Zero

from qdd import QddProvider


def test_circuit_sampler():
    # Computes an expectation value through trotterlization

    two_qubit_h2 = (-1.0523732 * I ^ I) + \
                   (0.39793742 * I ^ Z) + \
                   (-0.3979374 * Z ^ I) + \
                   (-0.0112801 * Z ^ Z) + \
                   (0.18093119 * X ^ X)

    evo_time = Parameter('θ')
    evolution_op = (evo_time * two_qubit_h2).exp_i()
    h2_measurement = StateFn(two_qubit_h2).adjoint()
    bell = CX @ (I ^ H) @ Zero
    evo_and_meas = h2_measurement @ evolution_op @ bell
    trotterized_op = PauliTrotterEvolution(trotter_mode=Suzuki(order=2, reps=1)).convert(evo_and_meas)
    diagonalized_meas_op = PauliExpectation().convert(trotterized_op)
    evo_time_points = list(range(8))
    h2_trotter_expectations = diagonalized_meas_op.bind_parameters({evo_time: evo_time_points})
    exact_eval = h2_trotter_expectations.eval()

    backend = QddProvider().get_backend()
    sampler = CircuitSampler(backend=backend)
    sampled_trotter_exp_op = sampler.convert(h2_trotter_expectations)
    sampled_eval = sampled_trotter_exp_op.eval()

    print(f'{exact_eval=}')
    print(f'{sampled_eval=}')
    assert exact_eval == pytest.approx(sampled_eval, abs=0.02)
