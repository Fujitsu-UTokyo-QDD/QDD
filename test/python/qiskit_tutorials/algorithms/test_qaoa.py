# The code in this file has been written using part of the code in the Qiskit tutorial:
# https://learning.quantum.ibm.com/tutorial/quantum-approximate-optimization-algorithm
import numpy as np
from qiskit.circuit.library import QAOAAnsatz
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from scipy.optimize import minimize
from typing import Sequence
import rustworkx as rx

from qdd import QddProvider
from qdd.qdd_estimator import Estimator
from qdd.qdd_sampler import Sampler


def test_qaoa():
    n = 5

    graph = rx.PyGraph()
    graph.add_nodes_from(np.arange(0, n, 1))
    edge_list = [
        (0, 1, 1.0),
        (0, 2, 1.0),
        (0, 4, 1.0),
        (1, 2, 1.0),
        (2, 3, 1.0),
        (3, 4, 1.0),
    ]
    graph.add_edges_from(edge_list)

    def build_max_cut_paulis(graph: rx.PyGraph):
        """Convert the graph to Pauli list.

        This function does the inverse of `build_max_cut_graph`
        """
        pauli_list = []
        for edge in list(graph.edge_list()):
            paulis = ["I"] * len(graph)
            paulis[edge[0]], paulis[edge[1]] = "Z", "Z"

            weight = graph.get_edge_data(edge[0], edge[1])

            pauli_list.append(("".join(paulis)[::-1], weight))

        return pauli_list

    max_cut_paulis = build_max_cut_paulis(graph)

    cost_hamiltonian = SparsePauliOp.from_list(max_cut_paulis)

    circuit = QAOAAnsatz(cost_operator=cost_hamiltonian, reps=2)
    circuit.measure_all()

    backend = QddProvider().get_backend()

    # Create pass manager for transpilation
    pm = generate_preset_pass_manager(optimization_level=3, backend=backend)

    candidate_circuit = pm.run(circuit)

    initial_gamma = np.pi
    initial_beta = np.pi / 2
    init_params = [initial_gamma, initial_beta, initial_gamma, initial_beta]

    def cost_func_estimator(params, ansatz, hamiltonian, estimator):

        # transform the observable defined on virtual qubits to
        # an observable defined on all physical qubits
        isa_hamiltonian = hamiltonian.apply_layout(ansatz.layout)

        job = estimator.run([ansatz], [isa_hamiltonian], [params])

        results = job.result()
        cost = results.values[0]

        return cost

    estimator = Estimator(run_options={"shots": None}, approximation=True)

    result = minimize(
        cost_func_estimator,
        init_params,
        args=(candidate_circuit, cost_hamiltonian, estimator),
        method="COBYLA",
        tol=1e-2,
    )

    optimized_circuit = candidate_circuit.assign_parameters(result.x)

    sampler = Sampler(run_options={"shots": None})

    job = sampler.run(optimized_circuit)
    result = job.result()
    final_distribution_int = result.quasi_dists[0]

    # auxiliary functions to sample most likely bitstring
    def to_bitstring(integer, num_bits):
        result = np.binary_repr(integer, width=num_bits)
        return [int(digit) for digit in result]

    keys = list(final_distribution_int.keys())
    values = list(final_distribution_int.values())
    most_likely = keys[np.argmax(np.abs(values))]
    most_likely_bitstring = to_bitstring(most_likely, len(graph))
    most_likely_bitstring.reverse()

    def evaluate_sample(x: Sequence[int], graph: rx.PyGraph) -> float:
        assert len(x) == len(
            list(graph.nodes())
        ), "The length of x must coincide with the number of nodes in the graph."
        return sum(
            x[u] * (1 - x[v]) + x[v] * (1 - x[u]) for u, v in list(graph.edge_list())
        )

    cut_value = evaluate_sample(most_likely_bitstring, graph)

    assert cut_value == 5
