import pytest
from qiskit import ClassicalRegister, QiskitError, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp

from qdd import QddProvider
from qdd.qdd_estimator import Estimator
from qdd.qdd_sampler import Sampler
from test.python.helpers.circuit_helper import assert_job_failed


def _assert_exact_measurement_probabilities(circuit, expected):
    backend = QddProvider().get_backend()
    backend_result = backend.run(circuit, shots=None).result()
    assert backend_result.success
    backend_probabilities = backend_result.data()["probabilities"]
    sampler_probabilities = (
        Sampler().run([circuit], is_exact=True).result()[0].data.quasi_dist
    )

    for probabilities in (backend_probabilities, sampler_probabilities):
        assert all(0 <= outcome < 2**circuit.num_clbits for outcome in probabilities)
        assert sum(probabilities.values()) == pytest.approx(1.0, abs=1e-12)
        nonzero = {outcome: p for outcome, p in probabilities.items() if p != 0}
        assert nonzero == pytest.approx(expected, abs=1e-12)


def test_exact_probabilities_partial_measurement():
    circuit = QuantumCircuit(2, 1)
    circuit.x(1)
    circuit.measure(1, 0)

    _assert_exact_measurement_probabilities(circuit, {1: 1.0})


def test_exact_probabilities_permuted_measurements():
    circuit = QuantumCircuit(2, 2)
    circuit.x(0)
    circuit.measure(0, 1)
    circuit.measure(1, 0)

    _assert_exact_measurement_probabilities(circuit, {2: 1.0})


def test_exact_probabilities_marginalize_unmeasured_qubits():
    circuit = QuantumCircuit(3, 1)
    circuit.h([0, 2])
    circuit.measure(2, 0)

    _assert_exact_measurement_probabilities(circuit, {0: 0.5, 1: 0.5})


def test_exact_probabilities_repeated_qubit_measurement():
    circuit = QuantumCircuit(2, 3)
    circuit.h([0, 1])
    circuit.measure(0, 0)
    circuit.measure(0, 2)

    # Both classical bits refer to the same quantum outcome and are correlated.
    _assert_exact_measurement_probabilities(circuit, {0: 0.5, 5: 0.5})


def test_exact_probabilities_overwritten_classical_bit():
    circuit = QuantumCircuit(2, 1)
    circuit.x(0)
    circuit.measure(0, 0)
    circuit.measure(1, 0)

    _assert_exact_measurement_probabilities(circuit, {0: 1.0})


def test_exact_probabilities_unwritten_classical_bits():
    circuit = QuantumCircuit(2, 4)
    circuit.x(1)
    circuit.measure(1, 2)

    _assert_exact_measurement_probabilities(circuit, {4: 1.0})


def test_exact_probabilities_identity_measurement_mapping():
    circuit = QuantumCircuit(2, 2)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.measure([0, 1], [0, 1])

    _assert_exact_measurement_probabilities(circuit, {0: 0.5, 3: 0.5})


def test_exact_probabilities_multiple_registers():
    circuit = QuantumCircuit(2)
    left = ClassicalRegister(1, "left")
    right = ClassicalRegister(2, "right")
    circuit.add_register(left, right)
    circuit.x(0)
    circuit.h(1)
    circuit.measure(0, left[0])
    circuit.measure(1, right[1])

    _assert_exact_measurement_probabilities(circuit, {1: 0.5, 5: 0.5})


def test_statevector_exact_probabilities_without_measurements():
    circuit = QuantumCircuit(2)
    circuit.x(1)
    backend = QddProvider().get_backend("statevector_simulator")

    result = backend.run(circuit, shots=None).result()

    assert result.success
    assert result.data()["probabilities"] == pytest.approx({0: 0, 1: 0, 2: 1, 3: 0})


def _dynamic_measurement_circuit(operation):
    circuit = QuantumCircuit(2, 2)
    circuit.h(0)
    circuit.cx(0, 1)
    if operation == "reset":
        circuit.reset(0)
    elif operation == "mid_measurement":
        circuit.measure(0, 0)
        circuit.x(0)
    elif operation == "classical_control":
        circuit.measure(0, 0)
        with circuit.if_test((circuit.clbits[0], 1)):
            circuit.x(1)
    circuit.measure([0, 1], [0, 1])
    return circuit


@pytest.mark.parametrize("backend_name", ["qasm_simulator", "statevector_simulator"])
@pytest.mark.parametrize("operation", ["reset", "mid_measurement", "classical_control"])
def test_exact_probabilities_reject_dynamic_circuits(backend_name, operation):
    circuit = _dynamic_measurement_circuit(operation)
    backend = QddProvider().get_backend(backend_name)

    job = backend.run(circuit, shots=None)

    assert_job_failed(job)
    assert "shots=None" in job.result().status
    assert "finite shot count" in job.result().status


@pytest.mark.parametrize("primitive", ["sampler", "estimator"])
def test_exact_primitives_report_unsupported_reset(primitive):
    circuit = _dynamic_measurement_circuit("reset")

    with pytest.raises(QiskitError, match="shots=None.*reset"):
        if primitive == "sampler":
            Sampler().run([circuit], is_exact=True).result()
        else:
            circuit.remove_final_measurements()
            Estimator().run([(circuit, [SparsePauliOp("ZI")])], precision=0).result()


def test_finite_shots_still_support_reset():
    circuit = _dynamic_measurement_circuit("reset")
    backend = QddProvider().get_backend()
    shots = 32

    result = backend.run(circuit, shots=shots).result()

    assert result.success
    counts = result.get_counts()
    assert sum(counts.values()) == shots
    assert set(counts) <= {"00", "10"}


def test_exact_probabilities_terminal_barriers():
    circuit = QuantumCircuit(2, 2)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.measure(0, 0)
    circuit.barrier()
    circuit.measure(1, 1)
    circuit.barrier()

    _assert_exact_measurement_probabilities(circuit, {0: 0.5, 3: 0.5})


def test_exact_probabilities_deferred_measurement():
    circuit = QuantumCircuit(2, 2)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.measure(0, 0)
    # This gate acts only on an unmeasured qubit, so q0's measurement can be deferred.
    circuit.x(1)
    circuit.measure(1, 1)

    _assert_exact_measurement_probabilities(circuit, {1: 0.5, 2: 0.5})


def test_exact_probabilities_reject_batch_with_dynamic_circuit():
    supported = QuantumCircuit(2, 2)
    supported.measure([0, 1], [0, 1])
    unsupported = _dynamic_measurement_circuit("reset")
    backend = QddProvider().get_backend()

    job = backend.run([supported, unsupported], shots=None)

    assert_job_failed(job)
    assert "shots=None" in job.result().status
