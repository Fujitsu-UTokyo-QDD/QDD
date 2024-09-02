# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/qiskit-community/qiskit-algorithms/blob/main/docs/tutorials/03_vqe_simulation_with_noise.ipynb

# This code is a part of a Qiskit project
# (C) Copyright IBM 2017, 2024.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import time

import pytest
from qiskit_algorithms import VQE, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import SPSA
from qiskit.circuit.library import TwoLocal
from qiskit.primitives import BaseEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.utils import algorithm_globals
from qiskit_aer.primitives import Estimator as AerEstimator

from qdd.qdd_estimator import Estimator


def test_vqe():
    """Tests behavior of evaluating circuits with Qiskit's library code.

    In this test, circuits will be evaluated in VQE#compute_minimum_eigenvalue(...).
    Note: VQE#compute_minimum_eigenvalue uses CircuitSampler to evaluate the specified circuit.
    """
    H2_op = SparsePauliOp.from_list(
        [
            ("II", -1.052373245772859),
            ("IZ", 0.39793742484318045),
            ("ZI", -0.39793742484318045),
            ("ZZ", -0.01128010425623538),
            ("XX", 0.18093119978423156),
        ]
    )

    # Computes the reference value of the minimum eigenvalue via a classical algorithm.
    npme = NumPyMinimumEigensolver()
    result = npme.compute_minimum_eigenvalue(operator=H2_op)
    ref_value = result.eigenvalue.real
    print(f"Reference value: {ref_value:.5f}")  # will be around -1.84237

    # Computes the minimum eigenvalue via VQE

    # Comment out Aer simulation because it is too slow
    # aer_backend = Aer.get_backend('aer_simulator')
    aer_result = _run_vqe(H2_op, AerEstimator(), ref_value)
    # aer_result = -1.853322535656245  # an actual result of Aer simulation

    qdd_result = _run_vqe(H2_op, Estimator(), ref_value)

    # Check the equality of the results from the two backends.
    # Almost certainly, the difference will be <0.01.
    # We relax the strictness of the closeness check to avoid test flakiness stemming from the random nature of VQE.
    assert aer_result == pytest.approx(qdd_result, abs=0.1)
    print(f"VQE result: Aer={aer_result}, Qdd={qdd_result}")


def _run_vqe(
    hamiltonian: SparsePauliOp, estimator: BaseEstimator, ref_value: float
) -> float:
    """Executes VQE to compute the minimum eigenvalue of the given hamiltonian, and returns the computed result.
    This method also prints the difference b/w the computed value and the given reference value.
    """

    seed = 170
    iterations = 125
    algorithm_globals.random_seed = seed

    def print_intermediate_result(eval_count, parameters, mean, std):
        if eval_count % 20 == 0:
            print(f"Progress: (count, energy) = ({eval_count, mean})")

    ansatz = TwoLocal(rotation_blocks="ry", entanglement_blocks="cz")
    spsa = SPSA(maxiter=iterations)
    vqe = VQE(estimator, ansatz, optimizer=spsa, callback=print_intermediate_result)
    start = time.time()
    result = vqe.compute_minimum_eigenvalue(operator=hamiltonian)
    print(f"Elapsed time: {time.time() - start}")
    print(f"VQE result is: {result.eigenvalue.real:.5f}")
    print(
        f"Delta from reference energy value is {(result.eigenvalue.real - ref_value):.5f}"
    )

    return result.eigenvalue.real
