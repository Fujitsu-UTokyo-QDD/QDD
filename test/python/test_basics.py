"""This file contains the basic tests.
Especially, the following test viewpoints are tested.
- #qubits in a circuit
  - Note: Test value of #qubits=#max-qubit is tested in production tests.
- #shots
- measurement targets
- memory
"""

import cmath, math
import pytest
from qiskit import QiskitError, QuantumCircuit, transpile
from qiskit_aer import Aer
from qiskit.circuit import Parameter
from qiskit.circuit.random import random_circuit
import numpy as np

from qdd import QddBackend, QddProvider
from test.python.helpers.circuit_helper import (
    assert_job_failed,
    get_oracle_counts_of_simple_circuit_run,
    run_simple_circuit,
)
from qdd import pyQDD


def test_num_qubits_2():
    """
    Note: This test also checks the followings test values.
    - #shots = 20
    - #measurements = #qubits (all qubits are measured)
    - memory = False
    """

    num_qubits = 2
    shots = 20
    counts = (
        run_simple_circuit(num_qubits=num_qubits, shots=shots).result().get_counts()
    )
    counts = {k: v for k, v in counts.items() if v != 0}
    expected = get_oracle_counts_of_simple_circuit_run(
        num_qubit=num_qubits, shots=shots
    )
    assert counts == expected


def test_num_qubits_max_plus_1():
    max_qubits = QddBackend._DEFAULT_CONFIG["n_qubits"]
    job = run_simple_circuit(num_qubits=max_qubits + 1, shots=20, skip_transpile=True)
    assert_job_failed(job)


def test_shots_minus_1():
    job = run_simple_circuit(num_qubits=2, shots=-1)
    assert_job_failed(job)


def test_shots_0():
    job = run_simple_circuit(num_qubits=2, shots=0)
    assert_job_failed(job)


def test_shots_max_shots():
    max_shots = QddBackend._DEFAULT_CONFIG["max_shots"]
    counts = run_simple_circuit(num_qubits=2, shots=max_shots).result().get_counts()
    expected = get_oracle_counts_of_simple_circuit_run(num_qubit=2, shots=max_shots)
    counts = {k: v for k, v in counts.items() if v != 0}
    assert counts == expected


def test_shots_max_shots_plus_1():
    max_shots = QddBackend._DEFAULT_CONFIG["max_shots"]
    if max_shots < 1024:
        assert (
            0
        ), "#max-shots is too small; The backend configuration value seems to be wrong."

    job = run_simple_circuit(num_qubits=2, shots=max_shots + 1)
    assert_job_failed(job)


def test_memory_true():
    """Tests the behavior of memory=True"""

    # Check whether AerSimulator and qddBackend return the same result of memory
    circ_x = QuantumCircuit(3)
    circ_x.x(0)
    circ_x.measure_all()

    aer_backend = Aer.get_backend("aer_simulator")
    aer_job = aer_backend.run(
        transpile(circuits=circ_x, backend=aer_backend, seed_transpiler=50),
        shots=20,
        memory=True,
        seed_simulator=80,
    )
    aer_memory = aer_job.result().get_memory()

    qdd_backend = QddProvider().get_backend()
    qdd_job = qdd_backend.run(
        transpile(circuits=circ_x, backend=qdd_backend, seed_transpiler=50),
        shots=20,
        memory=True,
        seed_simulator=80,
    )
    qdd_memory = qdd_job.result().get_memory()

    assert qdd_memory == aer_memory

    # Test the case where the elements in the memory list are not unique (i.e., at least one qubit is in superposition)
    circ_h = QuantumCircuit(3)
    circ_h.h(0)
    circ_h.measure_all()
    qdd_job_h = qdd_backend.run(
        transpile(circuits=circ_h, backend=qdd_backend, seed_transpiler=50),
        shots=10000,
        memory=True,
        seed_simulator=80,
    )
    qdd_memory_h = qdd_job_h.result().get_memory()

    assert "000" in qdd_memory_h
    assert "001" in qdd_memory_h
    assert len(qdd_memory_h) == 10000


def test_sv():
    """Tests the behavior of statevector simulator"""

    qdd_backend = QddProvider().get_backend("statevector_simulator")

    # Test the case where the elements in the memory list are not unique (i.e., at least one qubit is in superposition)
    circ_h = QuantumCircuit(3)
    circ_h.h(0)
    qdd_job_h = qdd_backend.run(
        transpile(circuits=circ_h, backend=qdd_backend, seed_transpiler=50),
        seed_simulator=80,
    )
    global_phase = qdd_job_h.result().results[0].header.global_phase
    sv = qdd_job_h.result().get_statevector()
    assert cmath.isclose(sv[0] * cmath.exp(1j * global_phase), 1 / math.sqrt(2))
    assert cmath.isclose(sv[1] * cmath.exp(1j * global_phase), 1 / math.sqrt(2))


def test_get_counts():
    """Tests behavior of Result.get_counts(...)."""

    # execute a single experiment
    qc1 = QuantumCircuit(2)
    qc1.x(0)
    qc1.measure_all()

    backend = QddProvider().get_backend()
    single_exp_job = backend.run(
        transpile(circuits=qc1, backend=backend, seed_transpiler=50),
        shots=20,
        seed_simulator=80,
    )
    single_exp_counts_with_circ = single_exp_job.result().get_counts(qc1)
    single_exp_counts_with_circ = {
        k: v for k, v in single_exp_counts_with_circ.items() if v != 0
    }
    single_exp_counts_without_circ = single_exp_job.result().get_counts()
    single_exp_counts_without_circ = {
        k: v for k, v in single_exp_counts_without_circ.items() if v != 0
    }
    print(
        f"Result: {single_exp_counts_with_circ=} and {single_exp_counts_without_circ}"
    )

    assert single_exp_counts_with_circ == single_exp_counts_without_circ == {"01": 20}

    # execute multiple experiments
    qc2 = QuantumCircuit(2)
    qc2.x(1)
    qc2.measure_all()

    qc_list = [qc1, qc2]
    multi_exp_job = backend.run(
        transpile(circuits=qc_list, backend=backend, seed_transpiler=50),
        shots=20,
        seed_simulator=80,
    )
    multi_exp_job_counts_all = multi_exp_job.result().get_counts()
    multi_exp_job_counts_all = [
        {k: v for k, v in conunts.items() if v != 0}
        for conunts in multi_exp_job_counts_all
    ]
    multi_exp_job_counts_qc1 = multi_exp_job.result().get_counts(qc1)
    multi_exp_job_counts_qc1 = {
        k: v for k, v in multi_exp_job_counts_qc1.items() if v != 0
    }
    multi_exp_job_counts_qc2 = multi_exp_job.result().get_counts(qc2)
    multi_exp_job_counts_qc2 = {
        k: v for k, v in multi_exp_job_counts_qc2.items() if v != 0
    }
    print(
        f"Result: {multi_exp_job_counts_all=}, {multi_exp_job_counts_qc1=} and {multi_exp_job_counts_qc2=}"
    )

    assert multi_exp_job_counts_qc1 == {"01": 20}
    assert multi_exp_job_counts_qc2 == {"10": 20}
    assert multi_exp_job_counts_all == [{"01": 20}, {"10": 20}]


def test_warn_unsupported_options():
    circ = QuantumCircuit(2)
    circ.x(0)
    circ.measure_all()
    backend = QddProvider().get_backend()

    with pytest.warns(UserWarning, match="unsupported_opt=foo"):
        job = backend.run(
            transpile(circuits=circ, backend=backend, seed_transpiler=50),
            shots=20,
            unsupported_opt="foo",
            seed_simulator=80,
        )
        print(job.result().get_counts())


def test_get_backend_by_name():
    backend = QddProvider().get_backend("qasm_simulator")
    assert type(backend) == QddBackend


def test_default_options():
    circ = QuantumCircuit(2)
    circ.x(0)
    circ.measure_all()

    backend = QddProvider().get_backend()
    job = backend.run(
        transpile(circuits=circ, backend=backend, seed_transpiler=50), seed_simulator=80
    )
    result = job.result()

    # default shots is 1024
    assert result.results[0].shots == 1024  # header
    counts = result.get_counts()
    assert sum(counts.values()) == 1024  # the number of sampled values

    # memory option is disabled by default
    with pytest.raises(QiskitError):
        result.get_memory()


def test_parametric_circuit_with_unbound_parameter():
    theta = Parameter("θ")
    n = 5
    qc = QuantumCircuit(n, 1)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)

    qc.barrier()
    qc.rz(theta, range(n))
    qc.barrier()

    for i in reversed(range(n - 1)):
        qc.cx(i, i + 1)
    qc.h(0)
    qc.measure(0, 0)

    backend = QddProvider().get_backend()
    job = backend.run(
        transpile(circuits=qc, backend=backend, seed_transpiler=50), seed_simulator=80
    )
    assert_job_failed(job)

def create_random(nQubits):
    circ = QuantumCircuit(nQubits);
    # TODO: random circuit creation
    return circ

def test_mv_multiply():
    for i in range(10):
        nQubits = 5
        circ = random_circuit(num_qubits=nQubits, depth=1,)
        circ.measure_all()

        backend = QddProvider().get_backend()
        circ2 = transpile(circ, backend=backend)
        circ2.measure_all()
        result = backend.run(circ2).result().to_dict()["results"][0]["edge"]    
        qdd_vector = result.getEigenVector()

        circ.save_statevector()
        aer_backend = Aer.get_backend("aer_simulator")
        circ2 = transpile(circ, backend=aer_backend)
        aer_result = aer_backend.run(circ2).result()
        aer_vector = aer_result.get_statevector(circ2)
        aer_vector_np = np.asarray(aer_vector)
        if np.allclose(qdd_vector, aer_vector_np) == False:
            print("###",i,"###")
            print(circ)
            print(qdd_vector)
            print(aer_vector_np)
            assert(False)

def test_mm_multiply_x():
    for _ in range(1):
        nQubits = 2
        circ = QuantumCircuit(nQubits)
        circ.x(1)
        #circ.cx(1,3)
        print(circ)

        backend = QddProvider().get_backend()
        result_medge = backend.merge_circuit(transpile(circ, backend=backend), 100000)
        qdd_unitary = result_medge.getEigenMatrix(nQubits)

        aer_backend = Aer.get_backend("unitary_simulator")
        circ2 = transpile(circ, backend=aer_backend)
        aer_result = aer_backend.run(circ2).result()
        aer_unitary = aer_result.get_unitary(circ2)
        aer_unitary_np = np.asarray(aer_unitary)
        if np.allclose(qdd_unitary, aer_unitary_np) == False:
            print(qdd_unitary)
            print(aer_unitary_np)
            assert(0)


def test_mm_multiply():
    backend = QddProvider().get_backend()
    aer_backend = Aer.get_backend("unitary_simulator")
    for i in range(10):
        nQubits = 10
        circ = random_circuit(num_qubits=nQubits, depth=2, max_operands=2)
        print(circ)

        result_medge = backend.merge_circuit(transpile(circ, backend=backend), 100000)
        qdd_unitary = result_medge.getEigenMatrix(nQubits)

        circ2 = transpile(circ, backend=aer_backend)
        aer_result = aer_backend.run(circ2).result()
        aer_unitary = aer_result.get_unitary(circ2)
        aer_unitary_np = np.asarray(aer_unitary)
        if np.allclose(qdd_unitary, aer_unitary_np, atol=0.1) == False:
            print("###", i, "###")
            print(qdd_unitary)
            print(aer_unitary_np)
            assert(0)
