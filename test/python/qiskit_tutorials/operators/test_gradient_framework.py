# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/operators/02_gradients_framework.ipynb  # noqa: E501

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
from qiskit import QuantumCircuit
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import CG
from qiskit.circuit import ParameterVector
from qiskit.opflow import Gradient, I, X, Z
from qiskit.utils import QuantumInstance

from qdd import QddProvider


@pytest.mark.slow
def test_vqe_with_gradient_based_optimization():
    # Instantiate the system Hamiltonian
    h2_hamiltonian = -1.05 * (I ^ I) + 0.39 * (I ^ Z) - 0.39 * (Z ^ I) - 0.01 * (Z ^ Z) + 0.18 * (X ^ X)

    # This is the target energy
    h2_energy = -1.85727503

    # Define the Ansatz
    wavefunction = QuantumCircuit(2)
    params = ParameterVector('theta', length=8)
    it = iter(params)
    wavefunction.ry(next(it), 0)
    wavefunction.ry(next(it), 1)
    wavefunction.rz(next(it), 0)
    wavefunction.rz(next(it), 1)
    wavefunction.cx(0, 1)
    wavefunction.ry(next(it), 0)
    wavefunction.ry(next(it), 1)
    wavefunction.rz(next(it), 0)
    wavefunction.rz(next(it), 1)

    grad = Gradient(grad_method='lin_comb')

    q_backend = QddProvider().get_backend()
    qi_sv = QuantumInstance(backend=q_backend, seed_transpiler=2, seed_simulator=80)

    # Conjugate Gradient algorithm
    optimizer = CG(maxiter=50)

    # Gradient callable
    vqe = VQE(wavefunction, optimizer=optimizer, gradient=grad, quantum_instance=qi_sv)

    result = vqe.compute_minimum_eigenvalue(h2_hamiltonian)
    print('Result:', result.optimal_value, 'Reference:', h2_energy)
    assert result.optimal_value == pytest.approx(h2_energy, abs=0.1)
