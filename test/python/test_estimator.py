# The code in this file has been written using part of the code in the Qiskit API documentation.
# https://docs.quantum.ibm.com/api/qiskit/qiskit.primitives.BaseEstimatorV1
from math import sqrt
from threading import Event
from unittest.mock import Mock
import warnings

import pytest
from qiskit import QuantumCircuit
from qiskit.primitives import StatevectorEstimator as QiskitEstimator
from qiskit.primitives.containers.estimator_pub import EstimatorPub
from qiskit.circuit.library import real_amplitudes
from qiskit.quantum_info import SparsePauliOp
from qdd.qdd_estimator import Estimator
from qdd import QddBackend


@pytest.mark.parametrize(
    ("pauli", "expected"), [("ZI", -1), ("IZ", 1), ("XI", 0), ("YI", 0), ("ZZ", -1)]
)
def test_exact_estimator_partial_observable(pauli, expected):
    circuit = QuantumCircuit(2)
    circuit.x(1)
    pubs = [(circuit, [SparsePauliOp(pauli)])]

    result = Estimator().run(pubs, precision=0).result()[0]
    reference = QiskitEstimator().run(pubs).result()[0]

    assert result.data.evs == pytest.approx(expected, abs=1e-12)
    assert result.data.evs == pytest.approx(reference.data.evs.item(), abs=1e-12)


@pytest.mark.parametrize("pub_format", ["tuple", "object"])
@pytest.mark.parametrize("run_precision", [None, 0.5])
def test_estimator_pub_precision_precedence(monkeypatch, pub_format, run_precision):
    run_options = {"shots": 99, "memory": True}
    estimator = Estimator(default_precision=0.125, run_options=run_options)
    circuit = QuantumCircuit(2)
    circuit.x(1)
    precisions = [0.25, None, 0.0, 0.5, 0.25]
    paulis = ["ZI", "IZ", "ZZ", "ZI", "ZI"]
    pubs = [
        (circuit, [SparsePauliOp(pauli)], None, precision)
        for pauli, precision in zip(paulis, precisions)
    ]
    if pub_format == "object":
        pubs = [EstimatorPub.coerce(pub) for pub in pubs]
    backend_run = Mock(wraps=estimator._backend.run)
    monkeypatch.setattr(estimator._backend, "run", backend_run)

    results = estimator.run(pubs, precision=run_precision).result()

    fallback = run_precision if run_precision is not None else 0.125
    expected_precisions = [0.25, fallback, 0.0, 0.5, 0.25]
    expected_shots = [16, 4 if run_precision is not None else 64, None, 4, 16]
    assert len(results) == len(pubs)
    for result, precision, shots, ev in zip(
        results, expected_precisions, expected_shots, [-1, 1, -1, -1, -1]
    ):
        assert result.data.metadata["shots"] == shots
        assert result.metadata["target_precision"] == precision
        assert result.metadata["shots"] == shots
        assert result.data.evs == pytest.approx(ev, abs=1e-12)
        assert result.data.stds == pytest.approx(0, abs=1e-12)
    assert [call.kwargs["shots"] for call in backend_run.call_args_list] == expected_shots
    assert all(call.kwargs["memory"] for call in backend_run.call_args_list)
    assert dict(estimator._run_options) == run_options

    following = estimator.run([(circuit, [SparsePauliOp("ZI")])]).result()[0]
    assert following.data.metadata["shots"] == 64


@pytest.mark.parametrize("precision", [0.25, 1e-300])
@pytest.mark.parametrize("source", ["run", "pub"])
def test_estimator_precision_shot_limit(monkeypatch, precision, source):
    monkeypatch.setattr(QddBackend, "_MAX_SHOTS", 4)
    estimator = Estimator()
    circuit = QuantumCircuit(2)
    circuit.h(1)
    backend_run = Mock(wraps=estimator._backend.run)
    monkeypatch.setattr(estimator._backend, "run", backend_run)
    pubs = [(circuit, [SparsePauliOp("ZI")], None, precision if source == "pub" else None)]

    with pytest.warns(UserWarning, match="maximum"):
        result = estimator.run(pubs, precision=precision if source == "run" else None).result()[0]

    assert backend_run.call_args.kwargs["shots"] == 4
    assert result.data.metadata["shots"] == 4
    assert result.metadata["shots"] == 4
    assert result.metadata["target_precision"] == precision
    assert result.data.stds == pytest.approx(sqrt(result.data.metadata["variance"] / 4))


@pytest.mark.parametrize(("precision", "shots"), [(0.5, 4), (1e300, 1)])
@pytest.mark.parametrize("source", ["run", "pub"])
def test_estimator_precision_shot_limit_boundaries(monkeypatch, precision, shots, source):
    monkeypatch.setattr(QddBackend, "_MAX_SHOTS", 4)
    estimator = Estimator()
    circuit = QuantumCircuit(2)
    circuit.x(1)
    pubs = [(circuit, [SparsePauliOp("ZI")], None, precision if source == "pub" else None)]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = estimator.run(pubs, precision=precision if source == "run" else None).result()[0]

    assert not caught
    assert result.data.metadata["shots"] == shots
    assert result.data.evs == pytest.approx(-1)


@pytest.mark.parametrize(("precision", "shots"), [(0.3, 12), (1 / sqrt(5), 5)])
def test_estimator_precision_shot_rounding(precision, shots):
    circuit = QuantumCircuit(2)

    result = Estimator().run(
        [(circuit, [SparsePauliOp("ZI")])], precision=precision
    ).result()[0]

    assert result.data.metadata["shots"] == shots


@pytest.mark.parametrize("precision", [-0.1, float("nan"), float("inf")])
@pytest.mark.parametrize("source", ["run", "pub", "default"])
def test_estimator_invalid_precision(precision, source):
    circuit = QuantumCircuit(2)
    pubs = [(circuit, [SparsePauliOp("ZI")])]
    estimator = Estimator(default_precision=precision if source == "default" else 0)
    if source == "pub":
        pubs = [(circuit, [SparsePauliOp("ZI")], None, precision)]

    with pytest.raises(ValueError):
        estimator.run(pubs, precision=precision if source == "run" else None).result()


def test_estimator_jobs_keep_their_precision(monkeypatch):
    estimator = Estimator()
    circuit = QuantumCircuit(2)
    circuit.x(1)
    pubs = [(circuit, [SparsePauliOp("ZI")])]
    first_started = Event()
    resume_first = Event()
    run = estimator._run

    def defer_first(*args, **kwargs):
        if args[0][0].precision == 0.25:
            first_started.set()
            if not resume_first.wait(timeout=10):
                raise TimeoutError("First estimator job was not resumed")
        return run(*args, **kwargs)

    monkeypatch.setattr(estimator, "_run", defer_first)
    first_job = estimator.run(pubs, precision=0.25)
    try:
        assert first_started.wait(timeout=5)
        second = estimator.run(pubs, precision=0.5).result()[0]
    finally:
        # Native simulations run sequentially; only job submission overlaps.
        resume_first.set()
    first = first_job.result()[0]

    assert first.data.metadata["shots"] == 16
    assert second.data.metadata["shots"] == 4


def test_estimator():
    psi1 = real_amplitudes(num_qubits=2, reps=2)

    H1 = SparsePauliOp.from_list([("II", 1), ("IZ", 2), ("XI", 3)])

    theta1 = [0, 1, 1, 2, 3, 5]

    pubs = [(psi1, [H1], theta1)]

    estimator = Estimator()
    estimator_qiskit = QiskitEstimator()

    # calculate [ <psi1(theta1)|H1|psi1(theta1)> ]
    job = estimator.run(pubs, precision=0.0)
    job_result = job.result()[0]  # It will block until the job finishes.
    print(f"The qdd-job with sampling finished with result {job_result}")
    job = estimator.run(pubs, precision=0.01)
    job_result_low_precision = job.result()[0]
    print(
        f"The qdd-job with approximation finished with result {job_result_low_precision}"
    )

    job = estimator_qiskit.run(pubs)
    job_result_qiskit = job.result()[0]
    print(f"The qiskit-job finished with result {job_result_qiskit}")

    assert job_result.data.evs == pytest.approx(job_result_qiskit.data.evs, rel=1e-6)

    std = job_result_low_precision.data.stds
    assert job_result_low_precision.data.evs == pytest.approx(
        job_result_qiskit.data.evs, rel=0.2, abs=std * 6
    )


def test_estimator2():
    # calculate [ <psi1(theta1)|H1|psi1(theta1)>,
    #             <psi2(theta2)|H2|psi2(theta2)>,
    #             <psi1(theta3)|H3|psi1(theta3)> ]

    psi1 = real_amplitudes(num_qubits=2, reps=2)
    psi2 = real_amplitudes(num_qubits=2, reps=3)

    H1 = SparsePauliOp.from_list([("II", 1), ("IZ", 2), ("XI", 3)])
    H2 = SparsePauliOp.from_list([("IZ", 1)])
    H3 = SparsePauliOp.from_list([("ZI", 1), ("ZZ", 1)])

    theta1 = [0, 1, 1, 2, 3, 5]
    theta2 = [0, 1, 1, 2, 3, 5, 8, 13]
    theta3 = [1, 2, 3, 4, 5, 6]

    pubs = [(psi1, [H1], theta1), (psi2, [H2], theta2), (psi1, [H3], theta3)]

    estimator = Estimator()
    estimator_qiskit = QiskitEstimator()

    job2 = estimator.run(pubs, precision=0.0)
    job_result = job2.result()
    print(f"The qdd-job with sampling finished with result {job_result}")
    job2 = estimator.run(pubs, precision=0.01)
    job_result_low_precision = job2.result()
    print(
        f"The qdd-job with approximation finished with result {job_result_low_precision}"
    )

    job2 = estimator_qiskit.run(pubs)
    job_result_qiskit = job2.result()
    print(f"The qiskit-job finished with result {job_result_qiskit}")

    for i in range(3):
        assert job_result[i].data.evs == pytest.approx(
            job_result_qiskit[i].data.evs, rel=1e-6
        )

    stds = [job_result_low_precision[i].data.stds for i in range(3)]
    for i in range(3):
        assert job_result_low_precision[i].data.evs == pytest.approx(
            job_result_qiskit[i].data.evs, rel=0.2, abs=stds[i] * 6
        )
