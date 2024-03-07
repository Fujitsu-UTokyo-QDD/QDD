# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/algorithms/02_vqe_convergence.ipynb  # noqa: E501

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

import logging

import numpy as np
import pytest
from qiskit.circuit.library import TwoLocal
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import VQE, NumPyMinimumEigensolver
from qiskit_algorithms.gradients import FiniteDiffEstimatorGradient
from qiskit_algorithms.optimizers import COBYLA, L_BFGS_B, SLSQP
from qiskit_algorithms.utils import algorithm_globals

from qdd.qdd_estimator_like_aer import Estimator


def test_vqe_convergence():
    H2_op = SparsePauliOp.from_list(
        [
            ("II", -1.052373245772859),
            ("IZ", 0.39793742484318045),
            ("ZI", -0.39793742484318045),
            ("ZZ", -0.01128010425623538),
            ("XX", 0.18093119978423156),
        ]
    )

    optimizers = [COBYLA(maxiter=80), L_BFGS_B(maxiter=60), SLSQP(maxiter=60)]
    converge_cnts = np.empty([len(optimizers)], dtype=object)
    converge_vals = np.empty([len(optimizers)], dtype=object)
    results = []
    for i, optimizer in enumerate(optimizers):
        print('\rOptimizer: {}        '.format(type(optimizer).__name__), end='')
        algorithm_globals.random_seed = 50
        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')

        counts = []
        values = []

        def store_intermediate_result(eval_count, parameters, mean, std):
            counts.append(eval_count)
            values.append(mean)

        vqe = VQE(Estimator(),ansatz, optimizer, callback=store_intermediate_result)
        result = vqe.compute_minimum_eigenvalue(operator=H2_op)
        converge_cnts[i] = np.asarray(counts)
        converge_vals[i] = np.asarray(values)
        results.append(result)
        print(result)
    print('\rOptimization complete')

    npme = NumPyMinimumEigensolver()
    result = npme.compute_minimum_eigenvalue(operator=H2_op)
    ref_value = result.eigenvalue.real
    print(f'Reference value: {ref_value:.5f}')

    print(f'Qdd result: COBYLA={results[0].eigenvalue.real}, LBFGS_B={results[1].eigenvalue.real},'
          f' SLSQP={results[2].eigenvalue.real}')

    assert results[0].eigenvalue.real == pytest.approx(ref_value, abs=0.1)
    # The optimizers of L_BFGS_B and SLSQP with sampling-based simulation produce results of poor precision;
    # so, we do not assert result[1] and result[2].


def test_vqe_with_gradient_framework_and_logging():
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('qiskit.algorithms.minimum_eigen_solvers.vqe').setLevel(logging.INFO)

    H2_op = SparsePauliOp.from_list(
        [
            ("II", -1.052373245772859),
            ("IZ", 0.39793742484318045),
            ("ZI", -0.39793742484318045),
            ("ZZ", -0.01128010425623538),
            ("XX", 0.18093119978423156),
        ]
    )

    algorithm_globals.random_seed = 50
    ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')

    # although the Qiskit tutorial uses SLSQP, we use COBYLA instead because SLSQP produces imprecise results
    optimizer = COBYLA(maxiter=80)

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    estimator = Estimator()
    gradient = FiniteDiffEstimatorGradient(estimator,epsilon=0.01)
    vqe = VQE(estimator,ansatz, optimizer, callback=store_intermediate_result,
              gradient=gradient)
    result = vqe.compute_minimum_eigenvalue(operator=H2_op)
    print(f'Value using Gradient: {result.eigenvalue.real:.5f}')

    print(result)

    reference_value = -1.85728
    assert result.eigenvalue.real == pytest.approx(reference_value, abs=0.1)

def test_vqe_with_gradient_framework_and_logging_statevector():
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('qiskit.algorithms.minimum_eigen_solvers.vqe').setLevel(logging.INFO)

    H2_op = SparsePauliOp.from_list(
        [
            ("II", -1.052373245772859),
            ("IZ", 0.39793742484318045),
            ("ZI", -0.39793742484318045),
            ("ZZ", -0.01128010425623538),
            ("XX", 0.18093119978423156),
        ]
    )

    algorithm_globals.random_seed = 50
    ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')

    # although the Qiskit tutorial uses SLSQP, we use COBYLA instead because SLSQP produces imprecise results
    optimizer = COBYLA(maxiter=80)

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    estimator = Estimator()
    gradient = FiniteDiffEstimatorGradient(estimator,epsilon=0.01)
    vqe = VQE(estimator,ansatz, optimizer, callback=store_intermediate_result,
              gradient=gradient)
    result = vqe.compute_minimum_eigenvalue(operator=H2_op)
    print(f'Value using Gradient: {result.eigenvalue.real:.5f}')

    print(result)

    reference_value = -1.85728
    assert result.eigenvalue.real == pytest.approx(reference_value, abs=0.1)
