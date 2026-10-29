# The code in this file has been written using part of the code in the Qiskit API documentation.
# https://docs.quantum.ibm.com/api/qiskit/qiskit.primitives.BaseSamplerV1
from collections import Counter

import pytest
from qiskit.primitives import StatevectorSampler as QiskitSampler
from qiskit import ClassicalRegister, QuantumCircuit
from qiskit.circuit.library import real_amplitudes
from qiskit.result import QuasiDistribution

from qdd.qdd_sampler import Sampler


def _basis_state_measurement_circuit(bits, register_name="meas"):
    # Keep at least two qubits so one-bit measurements also work with two MPI ranks.
    circuit = QuantumCircuit(max(2, len(bits)))
    circuit.add_register(ClassicalRegister(len(bits), register_name))
    for qubit, bit in enumerate(reversed(bits)):
        if bit == "1":
            circuit.x(qubit)
    circuit.measure(range(len(bits)), range(len(bits)))
    return circuit


@pytest.mark.parametrize(
    "bits",
    [
        "0",
        "1",
        "00",
        "100",
        "10000000",
        "000000000",
        "000000001",
        "100000000",
        "100000001",
        "10000000000000001",
    ],
)
def test_sampler_bit_array_width(bits):
    circuit = _basis_state_measurement_circuit(bits)
    shots = 8

    result = Sampler().run([circuit], shots=shots).result()[0]

    bit_array = result.data.meas
    assert bit_array.num_bits == len(bits)
    assert bit_array.num_shots == shots
    assert bit_array.get_bitstrings() == [bits] * shots
    assert bit_array.get_counts() == {bits: shots}
    assert result.data.quasi_dist.binary_probabilities(num_bits=len(bits)) == {bits: 1.0}


def test_sampler_bit_array_width_per_circuit():
    bitstrings = ["1", "100", "100000000"]
    circuits = [_basis_state_measurement_circuit(bits) for bits in bitstrings]
    shots = 8

    results = Sampler().run(circuits, shots=shots).result()

    assert len(results) == len(bitstrings)
    for result, bits in zip(results, bitstrings):
        assert result.data.meas.num_bits == len(bits)
        assert result.data.meas.get_bitstrings() == [bits] * shots


@pytest.mark.parametrize("register_name", ["readout", "quasi_dist"])
def test_sampler_register_name(register_name):
    circuit = _basis_state_measurement_circuit("100", register_name)
    shots = 8

    result = Sampler().run([circuit], shots=shots).result()[0]

    assert "meas" not in result.data
    assert result.data[register_name].num_bits == 3
    assert result.data[register_name].get_bitstrings() == ["100"] * shots
    if register_name == "quasi_dist":
        assert result.metadata["quasi_dist"].binary_probabilities() == {"100": 1.0}
    else:
        assert result.data.quasi_dist.binary_probabilities() == {"100": 1.0}


@pytest.mark.parametrize("with_reset", [False, True])
def test_sampler_multiple_registers(with_reset):
    circuit = QuantumCircuit(3)
    left = ClassicalRegister(2, "left")
    right = ClassicalRegister(9, "right")
    circuit.add_register(left, right)
    if with_reset:
        # Select the per-shot execution path without changing the initial state.
        circuit.reset(0)
    circuit.x([0, 2])
    circuit.measure(0, right[8])
    circuit.measure(0, right[0])
    circuit.measure(1, right[4])
    circuit.measure(2, left[1])
    circuit.measure(1, left[0])
    shots = 8

    result = Sampler().run([circuit], shots=shots).result()[0]
    reference = QiskitSampler().run([circuit], shots=shots).result()[0]

    assert set(result.data) == {"left", "right", "quasi_dist"}
    for name, width in [("left", 2), ("right", 9)]:
        assert result.data[name].num_bits == width
        assert result.data[name].num_shots == shots
        assert result.data[name].get_bitstrings() == reference.data[name].get_bitstrings()
    assert result.data.left.get_counts() == {"10": shots}
    assert result.data.right.get_counts() == {"100000001": shots}
    assert result.data.quasi_dist.binary_probabilities() == {"10000000110": 1.0}


def test_sampler_register_shot_alignment():
    circuit = QuantumCircuit(2)
    left = ClassicalRegister(1, "left")
    right = ClassicalRegister(1, "right")
    circuit.add_register(left, right)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.measure(0, left[0])
    circuit.measure(1, right[0])
    shots = 32

    result = Sampler().run([circuit], shots=shots).result()[0]

    left_bits = result.data.left.get_bitstrings()
    right_bits = result.data.right.get_bitstrings()
    assert left_bits == right_bits
    joint_counts = Counter(right + left for left, right in zip(left_bits, right_bits))
    assert result.data.quasi_dist.binary_probabilities(num_bits=2) == {
        bits: count / shots for bits, count in joint_counts.items()
    }


def test_sampler_register_names_per_circuit():
    circuits = [
        _basis_state_measurement_circuit("1", "first"),
        _basis_state_measurement_circuit("100000000", "second"),
    ]
    shots = 8

    results = Sampler().run(circuits, shots=shots).result()

    assert set(results[0].data) == {"first", "quasi_dist"}
    assert set(results[1].data) == {"second", "quasi_dist"}
    assert results[0].data.first.get_bitstrings() == ["1"] * shots
    assert results[1].data.second.get_bitstrings() == ["100000000"] * shots


def test_sampler_register_bit_mapping():
    circuit = QuantumCircuit(2)
    base = ClassicalRegister(3, "base")
    # Register-local bit order can differ from the circuit's classical bit order.
    alias = ClassicalRegister(bits=[base[2], base[0]], name="alias")
    circuit.add_register(base, alias)
    circuit.x(0)
    circuit.measure(0, base[2])
    circuit.measure(1, base[0])
    shots = 8

    result = Sampler().run([circuit], shots=shots).result()[0]

    assert result.data.base.num_bits == 3
    assert result.data.alias.num_bits == 2
    assert result.data.base.get_bitstrings() == ["100"] * shots
    assert result.data.alias.get_bitstrings() == ["01"] * shots


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
        assert dist == pytest.approx(dist_qiskit, rel=0.2, abs=0.015)
        assert dist_exact == pytest.approx(dist_qiskit, rel=0.2, abs=0.015)


def test_sampler_param():

    # two parameterized circuits
    pqc = real_amplitudes(num_qubits=2, reps=2)
    pqc.measure_all()
    pqc2 = real_amplitudes(num_qubits=2, reps=3)
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
        assert dist == pytest.approx(dist_qiskit, rel=0.2, abs=0.02)
        assert dist_exact == pytest.approx(dist_qiskit, rel=0.2, abs=0.02)
