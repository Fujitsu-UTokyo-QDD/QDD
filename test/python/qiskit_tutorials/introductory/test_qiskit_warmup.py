# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit/blob/7a7d6f19de4c68f4337370c631048ca18aa91c26/docs/intro_tutorial1.rst

# Copyright 2018 IBM and its contributors
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from qiskit import QuantumCircuit, transpile

from qdd import QddProvider


def test_warmup():
    simulator = QddProvider().get_backend()

    # Create a Quantum Circuit acting on the q register
    circuit = QuantumCircuit(2, 2)

    # Add a H gate on qubit 0
    circuit.h(0)

    # Add a CX (CNOT) gate on control qubit 0 and target qubit 1
    circuit.cx(0, 1)

    # Map the quantum measurement to the classical bits
    circuit.measure([0, 1], [0, 1])

    # compile the circuit down to low-level QASM instructions
    # supported by the backend (not needed for simple circuits)
    compiled_circuit = transpile(circuit, simulator, seed_transpiler=50)

    # Execute the circuit on the qasm simulator
    job = simulator.run(compiled_circuit, shots=1000, seed_simulator=80)

    # Grab results from the job
    result = job.result()

    # Returns counts
    counts = result.get_counts(compiled_circuit)
    print("\nTotal count for 00 and 11 are:", counts)
    assert '00' in counts
    assert '11' in counts
    assert counts["00"] + counts["11"] == 1000
