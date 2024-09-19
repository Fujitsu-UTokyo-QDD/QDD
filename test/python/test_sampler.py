# The code in this file has been written using part of the code in the Qiskit API documentation.
# https://docs.quantum.ibm.com/api/qiskit/qiskit.primitives.BaseSamplerV1
import pytest
from qiskit.primitives import StatevectorSampler as QiskitSampler
from qiskit import QuantumCircuit
from qiskit.circuit.library import RealAmplitudes
from qiskit.result import QuasiDistribution

from qdd.qdd_sampler import Sampler


def test_sampler():
    # a Bell circuit
    bell = QuantumCircuit(2)
    bell.h(0)
    bell.cx(0, 1)
    bell.measure_all()

    # initialization of the sampler
    sampler_qiskit = QiskitSampler()
    sampler = Sampler(
        run_options={"shots": 4096}
    )  # if shots is default value, tests sometimes fail because precision is not enough
    sampler_exact = Sampler(run_options={"shots": None})

    # Sampler runs a job on the Bell circuit
    job = sampler.run(circuits=[bell], parameter_values=[[]], parameters=[[]])
    job_result = job.result()
    print("QDD Sampler with shots:")
    print([q.binary_probabilities() for q in job_result.quasi_dists])
    job_exact = sampler_exact.run(
        circuits=[bell], parameter_values=[[]], parameters=[[]]
    )
    job_result_exact = job_exact.result()
    print("QDD Sampler without shots:")
    print([q.binary_probabilities() for q in job_result_exact.quasi_dists])
    pubs = [(bell, [])]
    job_qiskit = sampler_qiskit.run(pubs=pubs, shots=4096)
    qiskit_dists = []
    for result in job_qiskit.result():
        shots = result.metadata["shots"]
        qiskit_counts = result.data.meas.get_counts()
        quasidist_dict = {k: v / shots for k, v in qiskit_counts.items()}
        qiskit_dists.append(QuasiDistribution(quasidist_dict))
    print("Qiskit Sampler:")
    print([q.binary_probabilities() for q in qiskit_dists])
    for dist, dist_exact, dist_qiskit in zip(
        job_result.quasi_dists,
        job_result_exact.quasi_dists,
        qiskit_dists,
    ):
        # delete 0s from the distribution
        dist = {k: v for k, v in dist.items() if v != 0}
        dist_exact = {k: v for k, v in dist_exact.items() if v != 0}
        assert dist == pytest.approx(dist_qiskit, rel=0.2, abs=0.01)
        assert dist_exact == pytest.approx(dist_qiskit, rel=0.2, abs=0.01)


def test_sampler_param():

    # two parameterized circuits
    pqc = RealAmplitudes(num_qubits=2, reps=2)
    pqc.measure_all()
    pqc2 = RealAmplitudes(num_qubits=2, reps=3)
    pqc2.measure_all()

    theta1 = [0, 1, 1, 2, 3, 5]
    theta2 = [0, 1, 2, 3, 4, 5, 6, 7]

    # initialization of the sampler
    sampler_qiskit = QiskitSampler()
    sampler = Sampler(
        run_options={"shots": 4096}
    )  # if shots is default value, tests sometimes fail because precision is not enough
    sampler_exact = Sampler(run_options={"shots": None})

    # Sampler runs a job on the parameterized circuits
    job2 = sampler.run(
        circuits=[pqc, pqc2],
        parameter_values=[theta1, theta2],
        parameters=[pqc.parameters, pqc2.parameters],
    )
    job_result = job2.result()
    print("QDD Sampler with shots:")
    print([q.binary_probabilities() for q in job_result.quasi_dists])
    job2_exact = sampler_exact.run(
        circuits=[pqc, pqc2],
        parameter_values=[theta1, theta2],
        parameters=[pqc.parameters, pqc2.parameters],
    )
    job_result_exact = job2_exact.result()
    print("QDD Sampler without shots:")
    print([q.binary_probabilities() for q in job_result_exact.quasi_dists])
    pubs = [(pqc, theta1), (pqc2, theta2)]
    job2_qiskit = sampler_qiskit.run(pubs=pubs, shots=4096)
    job_result_qiskit = job2_qiskit.result()
    qiskit_dists = []
    for result in job_result_qiskit:
        shots = result.metadata["shots"]
        qiskit_counts = result.data.meas.get_counts()
        quasidist_dict = {k: v / shots for k, v in qiskit_counts.items()}
        qiskit_dists.append(QuasiDistribution(quasidist_dict))
    print("Qiskit Sampler:")
    print([q.binary_probabilities() for q in qiskit_dists])
    for dist, dist_exact, dist_qiskit in zip(
        job_result.quasi_dists,
        job_result_exact.quasi_dists,
        qiskit_dists,
    ):
        assert dist == pytest.approx(dist_qiskit, rel=0.2, abs=0.01)
        assert dist_exact == pytest.approx(dist_qiskit, rel=0.2, abs=0.01)
