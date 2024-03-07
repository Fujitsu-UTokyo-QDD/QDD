# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/algorithms/04_vqe_advanced.ipynb  # noqa: E501

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
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA, SLSQP
from qiskit.circuit.library import TwoLocal
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.utils import algorithm_globals

from qdd.qdd_estimator_like_aer import Estimator


class TestVQEAdvanced:
    def test_vqe_initial_points(self):
        H2_op = SparsePauliOp.from_list(
            [
                ("II", -1.052373245772859),
                ("IZ", 0.39793742484318045),
                ("ZI", -0.39793742484318045),
                ("ZZ", -0.01128010425623538),
                ("XX", 0.18093119978423156),
            ]
        )
        seed = 50
        algorithm_globals.random_seed = seed

        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
        optimizer = COBYLA(maxiter=100)
        vqe = VQE(Estimator(), ansatz, optimizer)
        result = vqe.compute_minimum_eigenvalue(operator=H2_op)
        print(result)
        optimizer_evals = result.cost_function_evals
        optimal_value = result.optimal_value

        # use the optimal points obtained in the above run as the initial points of a new VQE run
        initial_pt = result.optimal_point
        algorithm_globals.random_seed = seed

        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
        optimizer = COBYLA(maxiter=100)
        vqe = VQE(Estimator(),ansatz, optimizer, initial_point=initial_pt)
        result_with_initial_points = vqe.compute_minimum_eigenvalue(operator=H2_op)
        print(result_with_initial_points)
        optimizer_evals_with_initial_points = result_with_initial_points.cost_function_evals
        optimal_value_with_initial_points = result.optimal_value
        print(f'optimizer_evals is {optimizer_evals_with_initial_points} with initial point'
              f' versus {optimizer_evals} without it.')

        reference_value = -1.85728
        assert optimal_value == pytest.approx(reference_value, abs=0.1)
        assert optimal_value_with_initial_points == pytest.approx(reference_value, abs=0.1)

