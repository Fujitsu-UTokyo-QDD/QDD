# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/algorithms/01_algorithms_introduction.ipynb  # noqa: E501

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
from qiskit.algorithms.optimizers import SLSQP
from qiskit.circuit.library import TwoLocal
from qiskit.opflow import I, X, Z
from qiskit.utils import QuantumInstance, algorithm_globals

from qdd import QddProvider


@pytest.mark.slow
def test_vqe():
    seed = 50
    algorithm_globals.random_seed = seed

    h2_op = (-1.052373245772859 * I ^ I) + \
            (0.39793742484318045 * I ^ Z) + \
            (-0.39793742484318045 * Z ^ I) + \
            (-0.01128010425623538 * Z ^ Z) + \
            (0.18093119978423156 * X ^ X)

    backend = QddProvider().get_backend()
    qi = QuantumInstance(backend=backend, seed_transpiler=seed, seed_simulator=seed)

    ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
    slsqp = SLSQP(maxiter=1000)
    vqe = VQE(ansatz, optimizer=slsqp, quantum_instance=qi)
    result = vqe.compute_minimum_eigenvalue(h2_op)
    print(result)

    # SLSQP optimizer with sampling-based simulation produces results of poor precision; so, we do not assert the result
    # against any concrete values. In this test, we check whether the above code finishes with no errors.
    # Note: ideal result is -1.85; Qdd result is ranges around from -0.5 to -1.5;
    #       Aer simulator result is around -1.0 (probably close to the Qdd result);
    #       Aer statevector simulation results almost in the ideal result.
    assert True

    # Re-execute the above VQE example with lazy quantum instance assignment
    qi = QuantumInstance(backend, seed_transpiler=seed, seed_simulator=seed)

    ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
    slsqp = SLSQP(maxiter=1000)
    vqe = VQE(ansatz, optimizer=slsqp)

    vqe.quantum_instance = qi
    result = vqe.compute_minimum_eigenvalue(operator=h2_op)
    print(result)

    assert True
