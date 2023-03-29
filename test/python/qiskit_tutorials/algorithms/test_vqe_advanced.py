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
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import COBYLA, SLSQP
from qiskit.circuit.library import TwoLocal
from qiskit.opflow import AerPauliExpectation, I, PauliExpectation, X, Z
from qiskit.utils import QuantumInstance, algorithm_globals

from qdd import QddProvider


class TestVQEAdvanced:
    def test_vqe_initial_points(self):
        h2_op = (-1.052373245772859 * I ^ I) + \
                (0.39793742484318045 * I ^ Z) + \
                (-0.39793742484318045 * Z ^ I) + \
                (-0.01128010425623538 * Z ^ Z) + \
                (0.18093119978423156 * X ^ X)

        seed = 50
        algorithm_globals.random_seed = seed
        backend = QddProvider().get_backend()
        qi = QuantumInstance(backend=backend, seed_transpiler=seed, seed_simulator=seed)

        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
        optimizer = COBYLA(maxiter=100)
        vqe = VQE(ansatz, optimizer=optimizer, quantum_instance=qi)
        result = vqe.compute_minimum_eigenvalue(operator=h2_op)
        print(result)
        optimizer_evals = result.cost_function_evals
        optimal_value = result.optimal_value

        # use the optimal points obtained in the above run as the initial points of a new VQE run
        initial_pt = result.optimal_point
        algorithm_globals.random_seed = seed
        backend = QddProvider().get_backend()
        qi = QuantumInstance(backend=backend, seed_transpiler=seed, seed_simulator=seed)

        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
        optimizer = COBYLA(maxiter=100)
        vqe = VQE(ansatz, optimizer=optimizer, initial_point=initial_pt, quantum_instance=qi)
        result_with_initial_points = vqe.compute_minimum_eigenvalue(operator=h2_op)
        print(result_with_initial_points)
        optimizer_evals_with_initial_points = result_with_initial_points.cost_function_evals
        optimal_value_with_initial_points = result.optimal_value
        print(f'optimizer_evals is {optimizer_evals_with_initial_points} with initial point'
              f' versus {optimizer_evals} without it.')

        reference_value = -1.85728
        assert optimal_value == pytest.approx(reference_value, abs=0.1)
        assert optimal_value_with_initial_points == pytest.approx(reference_value, abs=0.1)

    def test_vqe_include_custom_flag(self):
        h2_op = (-1.052373245772859 * I ^ I) + \
                (0.39793742484318045 * I ^ Z) + \
                (-0.39793742484318045 * Z ^ I) + \
                (-0.01128010425623538 * Z ^ Z) + \
                (0.18093119978423156 * X ^ X)

        seed = 50
        algorithm_globals.random_seed = seed
        backend = QddProvider().get_backend()
        qi = QuantumInstance(backend=backend, seed_transpiler=seed, seed_simulator=seed)

        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
        optimizer = COBYLA(maxiter=100)
        # Note: include_custom=True tries to use AerPauliExpectation, but PauliExpectation will be used
        # because the type of the Qdd backend does not match any Aer simulator backends.
        vqe = VQE(ansatz, optimizer=optimizer, quantum_instance=qi, include_custom=True)
        result = vqe.compute_minimum_eigenvalue(operator=h2_op)
        optimal_value1 = result.optimal_value

        reference_value = -1.85728
        assert optimal_value1 == pytest.approx(reference_value, abs=0.1)

        # Although SLSQP and SPSA examples are shown in the Qiskit tutorial, we do not test them
        # because they has been tested in other test cases: test_vqe_**.py.

    def test_vqe_with_specified_expectation(self):
        h2_op = (-1.052373245772859 * I ^ I) + \
                (0.39793742484318045 * I ^ Z) + \
                (-0.39793742484318045 * Z ^ I) + \
                (-0.01128010425623538 * Z ^ Z) + \
                (0.18093119978423156 * X ^ X)

        seed = 50
        backend = QddProvider().get_backend()
        qi = QuantumInstance(backend=backend, seed_transpiler=seed, seed_simulator=seed)

        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
        optimizer = SLSQP(maxiter=1000)
        vqe = VQE(ansatz, optimizer=optimizer, quantum_instance=qi,
                  expectation=AerPauliExpectation())

        # cannot use AerPauliExpectation along with QddBackend
        with pytest.raises(Exception):
            vqe.compute_minimum_eigenvalue(operator=h2_op)

        optimizer = COBYLA(maxiter=1000)
        vqe = VQE(ansatz, optimizer=optimizer, quantum_instance=qi, expectation=PauliExpectation(group_paulis=False))
        result = vqe.compute_minimum_eigenvalue(operator=h2_op)
        print(result)

        reference_value = -1.85728
        assert result.optimal_value == pytest.approx(reference_value, abs=0.1)
