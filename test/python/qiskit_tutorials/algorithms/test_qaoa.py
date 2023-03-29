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

from collections import OrderedDict

import numpy as np
from qiskit.algorithms import QAOA, VQE, NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit.opflow import PauliSumOp, StateFn
from qiskit.quantum_info import Pauli
from qiskit.utils import QuantumInstance, algorithm_globals

from qdd import QddProvider


def test_qaoa():
    def get_operator(weight_matrix):
        r"""Generate Hamiltonian for the graph partitioning
        Notes:
            Goals:
                1 separate the vertices into two set of the same size
                2 make sure the number of edges between the two set is minimized.
            Hamiltonian:
                H = H_A + H_B
                H_A = sum\_{(i,j)\in E}{(1-ZiZj)/2}
                H_B = (sum_{i}{Zi})^2 = sum_{i}{Zi^2}+sum_{i!=j}{ZiZj}
                H_A is for achieving goal 2 and H_B is for achieving goal 1.
        Args:
            weight_matrix (numpy.ndarray) : adjacency matrix.
        Returns:
            PauliSumOp: operator for the Hamiltonian
            float: a constant shift for the obj function.
        """
        num_nodes = len(weight_matrix)
        pauli_list = []
        shift = 0

        for i in range(num_nodes):
            for j in range(i):
                if weight_matrix[i, j] != 0:
                    x_p = np.zeros(num_nodes, dtype=bool)
                    z_p = np.zeros(num_nodes, dtype=bool)
                    z_p[i] = True
                    z_p[j] = True
                    pauli_list.append([-0.5, Pauli((z_p, x_p))])
                    shift += 0.5

        for i in range(num_nodes):
            for j in range(num_nodes):
                if i != j:
                    x_p = np.zeros(num_nodes, dtype=bool)
                    z_p = np.zeros(num_nodes, dtype=bool)
                    z_p[i] = True
                    z_p[j] = True
                    pauli_list.append([1, Pauli((z_p, x_p))])
                else:
                    shift += 1

        pauli_list = [(pauli[1].to_label(), pauli[0]) for pauli in pauli_list]
        return PauliSumOp.from_list(pauli_list), shift

    def sample_most_likely(state_vector):
        """Compute the most likely binary string from state vector.
        Args:
            state_vector (numpy.ndarray or dict): state vector or counts.
        Returns:
            numpy.ndarray: binary string as numpy.ndarray of ints.
        """
        if isinstance(state_vector, (OrderedDict, dict)):
            # get the binary string with the largest count
            binary_string = sorted(state_vector.items(), key=lambda kv: kv[1])[-1][0]
            x = np.asarray([int(y) for y in reversed(list(binary_string))])
            return x
        elif isinstance(state_vector, StateFn):
            binary_string = list(state_vector.sample().keys())[0]
            x = np.asarray([int(y) for y in reversed(list(binary_string))])
            return x
        else:
            n = int(np.log2(state_vector.shape[0]))
            k = np.argmax(np.abs(state_vector))
            x = np.zeros(n)
            for i in range(n):
                x[i] = k % 2
                k >>= 1
            return x

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

    algorithm_globals.random_seed = 10598
    adjacency_matrix = np.array([[0., 1., 1., 0.],
                                 [1., 0., 1., 1.],
                                 [1., 1., 0., 1.],
                                 [0., 1., 1., 0.]])

    # compute via QAOA
    optimizer = COBYLA()
    backend = QddProvider().get_backend()
    qi = QuantumInstance(backend, seed_transpiler=50, seed_simulator=80)
    qaoa = QAOA(optimizer, quantum_instance=qi)
    qubit_op, offset = get_operator(adjacency_matrix)
    result = qaoa.compute_minimum_eigenvalue(qubit_op)
    x_most_likely = sample_most_likely(result.eigenstate)
    qaoa_result = objective_value(x_most_likely, adjacency_matrix)
    print(f'Objective value computed by QAOA is {qaoa_result}')

    # compute via a classical algorithm
    npme = NumPyMinimumEigensolver()
    result = npme.compute_minimum_eigenvalue(qubit_op)
    x_most_likely = sample_most_likely(result.eigenstate)
    reference_value = objective_value(x_most_likely, adjacency_matrix)
    print(f'Objective value computed by the NumPyMinimumEigensolver is {reference_value}')

    # compute via VQE
    optimizer = COBYLA()
    ansatz = TwoLocal(qubit_op.num_qubits, 'ry', 'cz', reps=5, entanglement='linear')
    backend = QddProvider().get_backend()
    qi = QuantumInstance(backend, seed_transpiler=50, seed_simulator=80)
    vqe = VQE(ansatz, optimizer, quantum_instance=qi)
    result = vqe.compute_minimum_eigenvalue(qubit_op)
    x_most_likely = sample_most_likely(result.eigenstate)
    vqe_result = objective_value(x_most_likely, adjacency_matrix)
    print(f'Objective value computed by VQE is {vqe_result}')

    assert qaoa_result == vqe_result == reference_value
