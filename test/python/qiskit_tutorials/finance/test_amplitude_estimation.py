# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://qiskit-community.github.io/qiskit-finance/tutorials/00_amplitude_estimation.html

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

import numpy as np
import pytest
from qiskit_algorithms import (
    AmplitudeEstimation,
    EstimationProblem,
    FasterAmplitudeEstimation,
    IterativeAmplitudeEstimation,
    MaximumLikelihoodAmplitudeEstimation,
)
from qiskit.circuit import QuantumCircuit

from qdd.qdd_sampler import Sampler


class BernoulliA(QuantumCircuit):
    """A circuit representing the Bernoulli A operator."""

    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(theta_p, 0)


class BernoulliQ(QuantumCircuit):
    """A circuit representing the Bernoulli Q operator."""

    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        self._theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(2 * self._theta_p, 0)

    def power(self, k, matrix_power=False):
        # implement the efficient power of Q
        q_k = QuantumCircuit(1)
        q_k.ry(2 * k * self._theta_p, 0)
        return q_k


@pytest.mark.no_mpi_support
def test_qae():
    sampler = Sampler()
    # Canonical AE
    p = 0.2
    circ_a = BernoulliA(p)
    circ_q = BernoulliQ(p)
    problem = EstimationProblem(
        state_preparation=circ_a,  # A operator
        grover_operator=circ_q,  # Q operator
        objective_qubits=[
            0
        ],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
    )

    ae = AmplitudeEstimation(
        num_eval_qubits=3,  # the number of evaluation qubits specifies circuit width and accuracy
        sampler=sampler,
    )
    ae_result = ae.estimate(problem)
    print("Raw estimate:", ae_result.estimation)
    print("Interpolated MLE estimator:", ae_result.mle)
    assert ae_result.mle == pytest.approx(0.2, abs=0.01)

    # Iterative Amplitude Estimation
    iae = IterativeAmplitudeEstimation(
        epsilon_target=0.01,  # target accuracy
        alpha=0.05,  # width of the confidence interval
        sampler=sampler,
    )
    iae_result = iae.estimate(problem)
    print("Estimate:", iae_result.estimation)
    assert iae_result.estimation == pytest.approx(0.2, abs=0.01)

    # Maximum Likelihood Amplitude Estimation
    mlae = MaximumLikelihoodAmplitudeEstimation(
        evaluation_schedule=3,  # log2 of the maximal Grover power
        sampler=sampler,
    )
    mlae_result = mlae.estimate(problem)
    print("Estimate:", mlae_result.estimation)
    assert mlae_result.estimation == pytest.approx(0.2, abs=0.01)

    # Faster Amplitude Estimation
    # redefine the problem without the grover operator because FasterAmplitudeEstimation with rescaling discards it
    problem = EstimationProblem(
        state_preparation=circ_a,  # A operator
        objective_qubits=[
            0
        ],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
    )
    fae = FasterAmplitudeEstimation(
        delta=0.01,  # target accuracy
        maxiter=3,  # determines the maximal power of the Grover operator
        sampler=sampler,
    )
    fae_result = fae.estimate(problem)
    print("Estimate:", fae_result.estimation)
    assert fae_result.estimation == pytest.approx(0.2, abs=0.01)
