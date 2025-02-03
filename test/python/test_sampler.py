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
    pubs = [(bell, [])]

    # initialization of the sampler
    sampler_qiskit = QiskitSampler()
    sampler = (
        Sampler()
    )  # if shots is default value, tests sometimes fail because precision is not enough

    # Sampler runs a job on the Bell circuit
    job = sampler.run(pubs=pubs, shots=4096)
    job_result = job.result()
    print("QDD Sampler with shots:")
    print([result.data.quasi_dist.binary_probabilities() for result in job_result])
    dists = [result.data.quasi_dist for result in job_result]
    job_exact = sampler.run(pubs=pubs, is_exact=True)
    job_result_exact = job_exact.result()
    print("QDD Sampler without shots:")
    print(
        [result.data.quasi_dist.binary_probabilities() for result in job_result_exact]
    )
    dists_exact = [result.data.quasi_dist for result in job_result_exact]
    job_qiskit = sampler_qiskit.run(pubs=pubs, shots=4096)
    qiskit_dists = []
    for result in job_qiskit.result():
        shots = result.metadata["shots"]
        qiskit_counts = result.data.meas.get_counts()
        quasidist_dict = {k: v / shots for k, v in qiskit_counts.items()}
        qiskit_dists.append(QuasiDistribution(quasidist_dict))
    print("Qiskit Sampler:")
    print([q.binary_probabilities() for q in qiskit_dists])
    for dist, dist_exact, dist_qiskit in zip(dists, dists_exact, qiskit_dists):
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

    pubs = [(pqc, theta1), (pqc2, theta2)]

    # initialization of the sampler
    sampler_qiskit = QiskitSampler()
    sampler = Sampler()

    # Sampler runs a job on the parameterized circuits
    job2 = sampler.run(
        pubs=pubs,
        shots=4096,
    )
    job_result = job2.result()
    print("QDD Sampler with shots:")
    print([result.data.quasi_dist.binary_probabilities() for result in job_result])
    dists = [result.data.quasi_dist for result in job_result]
    job2_exact = sampler.run(
        pubs=pubs,
        is_exact=True,
    )
    job_result_exact = job2_exact.result()
    print("QDD Sampler without shots:")
    print(
        [result.data.quasi_dist.binary_probabilities() for result in job_result_exact]
    )
    dists_exact = [result.data.quasi_dist for result in job_result_exact]
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
        dists,
        dists_exact,
        qiskit_dists,
    ):
        assert dist == pytest.approx(dist_qiskit, rel=0.2, abs=0.01)
        assert dist_exact == pytest.approx(dist_qiskit, rel=0.2, abs=0.01)
