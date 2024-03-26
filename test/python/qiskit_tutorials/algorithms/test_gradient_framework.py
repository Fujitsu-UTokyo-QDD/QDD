# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/qiskit-community/qiskit-algorithms/blob/main/docs/tutorials/12_gradients_framework.ipynb

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

import pytest
from qiskit import QuantumCircuit
from qiskit_algorithms import VQE
from qiskit_algorithms.gradients import LinCombEstimatorGradient
from qiskit_algorithms.optimizers import CG
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp

from qdd.qdd_estimator import Estimator


def test_vqe_with_gradient_based_optimization():
    # Instantiate the system Hamiltonian
    h2_hamiltonian = SparsePauliOp.from_list(
        [
            ("II", -1.05),
            ("IZ", 0.39),
            ("ZI", -0.39),
            ("ZZ", -0.01),
            ("XX", 0.18),
        ]
    )


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

    estimator = Estimator()
    grad = LinCombEstimatorGradient(estimator=estimator)

    # Conjugate Gradient algorithm
    optimizer = CG()
    #optimizer = COBYLA()

    # Gradient callable
    vqe = VQE(estimator=estimator, ansatz=wavefunction, optimizer=optimizer, gradient=grad)

    result = vqe.compute_minimum_eigenvalue(h2_hamiltonian)
    print('Result:', result.optimal_value, 'Reference:', h2_energy)
    assert result.optimal_value == pytest.approx(h2_energy, abs=0.1)
