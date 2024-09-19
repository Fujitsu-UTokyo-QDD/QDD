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

    estimator = Estimator(
        run_options={"shots": 4096}, approximation=False
    )  # if shots is default value, test sometimes fail because precision is not enough
    estimator_approx = Estimator(run_options={"shots": 4096}, approximation=True)
    estimator_exact = Estimator(run_options={"shots": None}, approximation=True)
    estimator_qiskit = QiskitEstimator()

    # calculate [ <psi1(theta1)|H1|psi1(theta1)> ]
    job = estimator.run([psi1], [H1], [theta1])
    job_result = job.result()  # It will block until the job finishes.
    print(f"The qdd-job with sampling finished with result {job_result}")
    job = estimator_approx.run([psi1], [H1], [theta1])
    job_result_approx = job.result()
    print(f"The qdd-job with approximation finished with result {job_result_approx}")
    job = estimator_exact.run([psi1], [H1], [theta1])
    job_result_exact = job.result()
    print(f"The qdd-job without sampling finished with result {job_result_exact}")
    job = estimator_qiskit.run(pubs=[(psi1, H1, theta1)])
    job_result_qiskit = job.result()[0]
    print(f"The qiskit-job finished with result {job_result_qiskit}")
    var = job_result.metadata[0]["variance"]
    std = sqrt(var / 4096)
    assert job_result.values[0] == pytest.approx(
        job_result_qiskit.data.evs, rel=0.2, abs=std * 6
    )
    var = job_result_approx.metadata[0]["variance"]
    std = sqrt(var / 4096)
    assert job_result_approx.values[0] == pytest.approx(
        job_result_qiskit.data.evs, rel=0.2, abs=std * 6
    )
    assert job_result_exact.values[0] == pytest.approx(
        job_result_qiskit.data.evs, rel=1e-6
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

    estimator = Estimator(
        run_options={"shots": 4096}, approximation=False
    )  # if shots is default value, test sometimes fail because precision is not enough
    estimator_approx = Estimator(run_options={"shots": 4096}, approximation=True)
    estimator_exact = Estimator(run_options={"shots": None}, approximation=True)
    estimator_qiskit = QiskitEstimator()

    job2 = estimator.run([psi1, psi2, psi1], [H1, H2, H3], [theta1, theta2, theta3])
    job_result = job2.result()
    print(f"The qdd-job with sampling finished with result {job_result}")
    job2 = estimator_approx.run(
        [psi1, psi2, psi1], [H1, H2, H3], [theta1, theta2, theta3]
    )
    job_result_approx = job2.result()
    print(f"The qdd-job with approximation finished with result {job_result_approx}")
    job2 = estimator_exact.run(
        [psi1, psi2, psi1], [H1, H2, H3], [theta1, theta2, theta3]
    )
    job_result_exact = job2.result()
    print(f"The qdd-job without sampling finished with result {job_result_exact}")
    pubs = [(psi1, H1, theta1), (psi2, H2, theta2), (psi1, H3, theta3)]
    job2 = estimator_qiskit.run(pubs=pubs)
    job_result_qiskit = job2.result()
    print(f"The qiskit-job finished with result {job_result_qiskit}")

    vars = [job_result.metadata[i]["variance"] for i in range(3)]
    stds = [sqrt(var / 4096) for var in vars]
    for i in range(3):
        assert job_result.values[i] == pytest.approx(
            job_result_qiskit[i].data.evs, rel=0.2, abs=stds[i] * 6
        )
    vars = [job_result_approx.metadata[i]["variance"] for i in range(3)]
    stds = [sqrt(var / 4096) for var in vars]
    for i in range(3):
        assert job_result_approx.values[i] == pytest.approx(
            job_result_qiskit[i].data.evs, rel=0.2, abs=stds[i] * 6
        )
    for i in range(3):
        assert job_result_exact.values[i] == pytest.approx(
            job_result_qiskit[i].data.evs, rel=1e-6
        )
