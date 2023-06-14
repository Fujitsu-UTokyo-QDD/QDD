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
from qiskit.algorithms import VQE, NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import COBYLA, L_BFGS_B, SLSQP
from qiskit.circuit.library import TwoLocal
from qiskit.opflow import Gradient, I, X, Z
from qiskit.utils import QuantumInstance, algorithm_globals

from qdd import QddProvider


@pytest.mark.slow
def test_vqe_convergence():
    h2_op = (-1.052373245772859 * I ^ I) + \
            (0.39793742484318045 * I ^ Z) + \
            (-0.39793742484318045 * Z ^ I) + \
            (-0.01128010425623538 * Z ^ Z) + \
            (0.18093119978423156 * X ^ X)

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

        backend = QddProvider().get_backend()
        vqe = VQE(ansatz, optimizer, callback=store_intermediate_result,
                  quantum_instance=QuantumInstance(backend=backend, seed_transpiler=50, seed_simulator=80))
        result = vqe.compute_minimum_eigenvalue(operator=h2_op)
        converge_cnts[i] = np.asarray(counts)
        converge_vals[i] = np.asarray(values)
        results.append(result)
        print(result)
    print('\rOptimization complete')

    npme = NumPyMinimumEigensolver()
    result = npme.compute_minimum_eigenvalue(operator=h2_op)
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

    h2_op = (-1.052373245772859 * I ^ I) + \
            (0.39793742484318045 * I ^ Z) + \
            (-0.39793742484318045 * Z ^ I) + \
            (-0.01128010425623538 * Z ^ Z) + \
            (0.18093119978423156 * X ^ X)

    algorithm_globals.random_seed = 50
    ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')

    # although the Qiskit tutorial uses SLSQP, we use COBYLA instead because SLSQP produces imprecise results
    optimizer = COBYLA(maxiter=80)

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    backend = QddProvider().get_backend()
    vqe = VQE(ansatz, optimizer, callback=store_intermediate_result,
              gradient=Gradient(grad_method='fin_diff'),
              quantum_instance=QuantumInstance(backend=backend, seed_transpiler=50, seed_simulator=80))
    result = vqe.compute_minimum_eigenvalue(operator=h2_op)
    print(f'Value using Gradient: {result.eigenvalue.real:.5f}')

    print(result)

    reference_value = -1.85728
    assert result.eigenvalue.real == pytest.approx(reference_value, abs=0.1)

def test_vqe_with_gradient_framework_and_logging_statevector():
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('qiskit.algorithms.minimum_eigen_solvers.vqe').setLevel(logging.INFO)

    h2_op = (-1.052373245772859 * I ^ I) + \
            (0.39793742484318045 * I ^ Z) + \
            (-0.39793742484318045 * Z ^ I) + \
            (-0.01128010425623538 * Z ^ Z) + \
            (0.18093119978423156 * X ^ X)

    algorithm_globals.random_seed = 50
    ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')

    # although the Qiskit tutorial uses SLSQP, we use COBYLA instead because SLSQP produces imprecise results
    optimizer = COBYLA(maxiter=80)

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    backend = QddProvider().get_backend("statevector_simulator")
    vqe = VQE(ansatz, optimizer, callback=store_intermediate_result,
              gradient=Gradient(grad_method='fin_diff'),
              quantum_instance=QuantumInstance(backend=backend, seed_transpiler=50, seed_simulator=80))
    result = vqe.compute_minimum_eigenvalue(operator=h2_op)
    print(f'Value using Gradient: {result.eigenvalue.real:.5f}')

    print(result)

    reference_value = -1.85728
    assert result.eigenvalue.real == pytest.approx(reference_value, abs=0.1)