import pytest
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.primitives import BackendSamplerV2

from qdd import QddProvider


@pytest.mark.parametrize("bits", ["000", "001", "11", "100", "10000000", "100000000"])
@pytest.mark.parametrize("with_reset", [False, True])
def test_backend_memory_format(bits, with_reset):
    """Raw memory must work with both Result and the standard V2 wrapper."""
    circuit = QuantumCircuit(len(bits))
    if with_reset:
        # Resetting the initial |0> leaves the expected output unchanged but
        # selects the backend's per-shot execution path instead of sampling
        # a shared final state, so both paths' raw memory formats are tested.
        circuit.reset(0)
    for qubit, bit in enumerate(reversed(bits)):
        if bit == "1":
            circuit.x(qubit)
    circuit.measure_all()

    shots = 8
    backend = QddProvider().get_backend()
    result = backend.run(circuit, shots=shots, memory=True).result()

    assert result.data()["memory"] == [hex(int(bits, 2))] * shots
    assert result.get_memory() == [bits] * shots
    assert result.get_counts() == {bits: shots}

    pub_result = BackendSamplerV2(backend=backend).run([circuit], shots=shots).result()[0]
    assert pub_result.data.meas.num_bits == len(bits)
    assert pub_result.data.meas.get_bitstrings() == [bits] * shots
    assert pub_result.data.meas.get_counts() == {bits: shots}


@pytest.mark.parametrize("with_reset", [False, True])
def test_backend_memory_register_mapping(with_reset):
    """Hex memory preserves classical bit order and register boundaries."""
    qubits = QuantumRegister(2, "q")
    left = ClassicalRegister(1, "left")
    right = ClassicalRegister(2, "right")
    circuit = QuantumCircuit(qubits, left, right)
    if with_reset:
        # Exercise the per-shot execution path without changing the expected
        # output: the qubit is already |0> before this reset.
        circuit.reset(qubits[0])
    circuit.x(qubits)
    circuit.measure(qubits[1], left[0])
    circuit.measure(qubits[0], right[1])

    shots = 8
    backend = QddProvider().get_backend()
    result = backend.run(circuit, shots=shots, memory=True).result()

    assert result.data()["memory"] == ["0x5"] * shots
    assert result.get_memory() == ["10 1"] * shots
    assert result.get_counts() == {"10 1": shots}

    pub_result = BackendSamplerV2(backend=backend).run([circuit], shots=shots).result()[0]
    assert pub_result.data.left.get_bitstrings() == ["1"] * shots
    assert pub_result.data.right.get_bitstrings() == ["10"] * shots


def test_statevector_memory_without_classical_bits():
    circuit = QuantumCircuit(2)
    circuit.x(0)
    backend = QddProvider().get_backend("statevector_simulator")

    result = backend.run(circuit, shots=8, memory=True).result()

    assert result.data()["memory"] == ["0x0"] * 8
