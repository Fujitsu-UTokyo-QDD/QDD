# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/algorithms/05_qaoa.ipynb  # noqa: E501

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


import numpy as np
from qiskit_algorithms import QAOA, SamplingVQE, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.result import QuasiDistribution
from qiskit_algorithms.utils import algorithm_globals
from qiskit.primitives import Sampler as QiskitSampler

from qdd.qdd_sampler_like_aer import Sampler


def test_qaoa():
    def get_operator(weight_matrix):
        r"""Generate Hamiltonian for the graph partitioning
        Notes:
            Goals:
                1 Separate the vertices into two set of the same size.
                2 Make sure the number of edges between the two set is minimized.
            Hamiltonian:
                H = H_A + H_B
                H_A = sum\_{(i,j)\in E}{(1-ZiZj)/2}
                H_B = (sum_{i}{Zi})^2 = sum_{i}{Zi^2}+sum_{i!=j}{ZiZj}
                H_A is for achieving goal 2 and H_B is for achieving goal 1.
        Args:
            weight_matrix: Adjacency matrix.
        Returns:
            Operator for the Hamiltonian
        A constant shift for the obj function.
        """
        num_nodes = len(weight_matrix)
        pauli_list = []
        coeffs = []
        shift = 0

        for i in range(num_nodes):
            for j in range(i):
                if weight_matrix[i, j] != 0:
                    x_p = np.zeros(num_nodes, dtype=bool)
                    z_p = np.zeros(num_nodes, dtype=bool)
                    z_p[i] = True
                    z_p[j] = True
                    pauli_list.append(Pauli((z_p, x_p)))
                    coeffs.append(-0.5)
                    shift += 0.5

        for i in range(num_nodes):
            for j in range(num_nodes):
                if i != j:
                    x_p = np.zeros(num_nodes, dtype=bool)
                    z_p = np.zeros(num_nodes, dtype=bool)
                    z_p[i] = True
                    z_p[j] = True
                    pauli_list.append(Pauli((z_p, x_p)))
                    coeffs.append(1.0)
                else:
                    shift += 1

        return SparsePauliOp(pauli_list, coeffs=coeffs), shift

    def sample_most_likely(state_vector):
        """Compute the most likely binary string from state vector.
        Args:
            state_vector (numpy.ndarray or dict): state vector or counts.
        Returns:
            numpy.ndarray: binary string as numpy.ndarray of ints.
        """
        if isinstance(state_vector, QuasiDistribution):
            values = list(state_vector.values())
        else:
            values = state_vector
        n = 4
        k = np.argmax(np.abs(values))
        x = bitfield(k,n)
        x.reverse()
        return np.asarray(x)

    def objective_value(x, w):
        """Compute the value of a cut.
        Args:
            x (numpy.ndarray): binary string as numpy array.
            w (numpy.ndarray): adjacency matrix.
        Returns:
            float: value of the cut.
        """
        outer_prod = np.outer(x, (1 - x))
        w_01 = np.where(w != 0, 1, 0)
        return np.sum(w_01 * outer_prod)

    def bitfield(n, L):
        result = np.binary_repr(n,L)
        return [int(digit) for digit in result] # [2:] to chop off the "0b" part

    algorithm_globals.random_seed = 10598
    adjacency_matrix = np.array([[0., 1., 1., 0.],
                                 [1., 0., 1., 1.],
                                 [1., 1., 0., 1.],
                                 [0., 1., 1., 0.]])

    # compute via QAOA
    optimizer = COBYLA()
    qaoa = QAOA(sampler=Sampler(run_options={"seed_simulator":80},transpile_options={"seed_transpiler":50}),optimizer=optimizer,reps=8)
    qubit_op, offset = get_operator(adjacency_matrix)
    result = qaoa.compute_minimum_eigenvalue(qubit_op)
    x_most_likely = sample_most_likely(result.eigenstate)
    qaoa_result = objective_value(x_most_likely, adjacency_matrix)
    print(f'Objective value computed by QAOA is {qaoa_result}')

    # compute via QAOA (Qiskit)
    optimizer = COBYLA()
    qaoa_qiskit = QAOA(sampler=QiskitSampler(),optimizer=optimizer,reps=2)
    qubit_op_qiskit, offset = get_operator(adjacency_matrix)
    result_qiskit = qaoa_qiskit.compute_minimum_eigenvalue(qubit_op_qiskit)
    x_most_likely_aer = sample_most_likely(result_qiskit.eigenstate)
    qaoa_result_qiskit = objective_value(x_most_likely_aer, adjacency_matrix)
    print(f'Objective value computed by QAOA_qiskit is {qaoa_result_qiskit}')

    assert qaoa_result == qaoa_result_qiskit

    # compute via a classical algorithm
    npme = NumPyMinimumEigensolver()
    result = npme.compute_minimum_eigenvalue(qubit_op)
    x_most_likely = sample_most_likely(result.eigenstate)
    reference_value = objective_value(x_most_likely, adjacency_matrix)
    print(f'Objective value computed by the NumPyMinimumEigensolver is {reference_value}')

    # compute via VQE
    optimizer = COBYLA()
    ansatz = TwoLocal(qubit_op.num_qubits, 'ry', 'cz', reps=2, entanglement='linear')
    vqe = SamplingVQE(sampler=Sampler(run_options={"seed_simulator":80},transpile_options={"seed_transpiler":50}),ansatz=ansatz,optimizer=optimizer)
    result = vqe.compute_minimum_eigenvalue(qubit_op)
    x_most_likely = sample_most_likely(result.eigenstate)
    vqe_result = objective_value(x_most_likely, adjacency_matrix)
    print(f'Objective value computed by VQE is {vqe_result}')

    assert vqe_result == reference_value
