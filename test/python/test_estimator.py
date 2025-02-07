# The code in this file has been written using part of the code in the Qiskit API documentation.
# https://docs.quantum.ibm.com/api/qiskit/qiskit.primitives.BaseEstimatorV1
from math import sqrt
import pytest
from qiskit.primitives import StatevectorEstimator as QiskitEstimator
from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import SparsePauliOp
from qdd.qdd_estimator import Estimator


def test_estimator():
    psi1 = RealAmplitudes(num_qubits=2, reps=2)

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

    psi1 = RealAmplitudes(num_qubits=2, reps=2)
    psi2 = RealAmplitudes(num_qubits=2, reps=3)

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
