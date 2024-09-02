from itertools import chain
from typing import Dict

import pytest
from qiskit import QuantumCircuit, transpile
from qiskit.providers import Backend, JobStatus, JobV1

from qdd import QddProvider


def get_simple_circuit(num_qubits: int) -> QuantumCircuit:
    """Returns a simple circuit that contains no stochastic operations."""

    circ = QuantumCircuit(num_qubits)
    for i in range(num_qubits):
        circ.x(i)
    for i in range(num_qubits):
        if i % 2 != 0:
            circ.cx(i - 1, i)
    circ.measure_all()
    return circ


def run_simple_circuit(num_qubits: int, shots: int) -> JobV1:
    """Executes a simple circuit that contains no stochastic operations and returns a 'counts' result."""

    circ = get_simple_circuit(num_qubits)
    backend = QddProvider().get_backend()
    job = backend.run(
        transpile(circuits=circ, backend=backend, seed_transpiler=50),
        shots=shots,
        seed_simulator=80,
    )
    return job


def get_oracle_counts_of_simple_circuit_run(num_qubit: int, shots: int) -> dict:
    """Returns the counts (oracle value) obtained via measurement of the circuit created by run_simple_circuit(...)."""

    if shots <= 0:
        assert 0, "Illegal arguments. shots must be positive."

    reversed_measured_values = []
    for i in range(num_qubit):
        if i % 2 == 0:
            reversed_measured_values.append("1")
        else:
            reversed_measured_values.append("0")

    measured_value_binary = "".join(reversed(reversed_measured_values))
    return {measured_value_binary: shots}


def get_counts(
    circuit: QuantumCircuit, backend: Backend, n_shots: int, optimization_level: int = 1
):
    """Simulates the given circuit via the given backend {n_shots} times.
    Returns:
        (dict, dict): the first dict is raw data (i.e., result.data()['counts']).
                      the second one is formatted data (i.e., result.get_counts()).
    """

    job = backend.run(
        transpile(
            circuits=circuit,
            backend=backend,
            optimization_level=optimization_level,
            seed_transpiler=50,
        ),
        shots=n_shots,
        seed_simulator=80,
    )
    result = job.result()
    counts = result.data()["counts"]
    formatted_counts = result.get_counts()

    return counts, formatted_counts


def assert_probabilities_are_close(
    counts_1st: Dict[str, int], counts_2nd: Dict[str, int], atol=0.05
):
    measured_count = sum(counts_1st.values())
    counts_1st_sorted = sorted(counts_1st.items())
    counts_2nd_sorted = sorted(counts_2nd.items())

    orig_prob_1st = {kv[0]: kv[1] / float(measured_count) for kv in counts_1st_sorted}
    prob_1st = orig_prob_1st.copy()
    orig_prob_2nd = {kv[0]: kv[1] / float(measured_count) for kv in counts_2nd_sorted}
    prob_2nd = orig_prob_2nd.copy()

    # If the probability for a measured value is small, the value may not be contained in the counts dicts.
    # To avoid assertion error due to missing keys, we ignore values measured with small probabilities (< 0.01).
    kv_1st_of_small_prob = filter(lambda kv: kv[1] < 0.01, prob_1st.items())
    kv_2nd_of_small_prob = filter(lambda kv: kv[1] < 0.01, prob_2nd.items())
    for kv_of_small_prob in set(chain(kv_1st_of_small_prob, kv_2nd_of_small_prob)):
        prob_1st.pop(kv_of_small_prob[0], None)
        prob_2nd.pop(kv_of_small_prob[0], None)

    assert prob_1st == pytest.approx(prob_2nd, abs=atol)


def assert_job_failed(job: JobV1):
    assert job.status() == JobStatus.ERROR
    result = job.result()
    assert result.success is False
    assert len(result.results) == 0
    assert len(result.status) > 0  # checks the error message is not empty


def assert_job_succeeded(job: JobV1):
    assert job.status() == JobStatus.DONE
    result = job.result()
    assert result.success is True
    assert len(result.results) >= 1
