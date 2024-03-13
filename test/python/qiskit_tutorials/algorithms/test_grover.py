# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/qiskit-community/qiskit-algorithms/blob/main/docs/tutorials/06_grover.ipynb

# This code is a part of a Qiskit project
# (C) Copyright IBM 2017, 2024.
# 
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
# 
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import numpy as np
import pytest
from qiskit import MissingOptionalLibraryError, QuantumCircuit
from qiskit_algorithms import AmplificationProblem, Grover
from qiskit.circuit.library import GroverOperator, PhaseOracle
from qiskit.quantum_info import Statevector

from qdd.qdd_sampler import Sampler


def test_grover():
    # the state we desire to find is '11'
    good_state = ['11']

    # define an oracle with QuantumCircuit
    oracle = QuantumCircuit(2)
    oracle.cz(0, 1)

    sampler = Sampler(run_options={"shots":1024,"seed_simulator":80},transpile_options={"seed_transpiler":50})

    problem = AmplificationProblem(oracle, is_good_state=good_state)
    grover = Grover(sampler=sampler)
    result = grover.amplify(problem)
    print('Result type:', type(result))
    print('Success!' if result.oracle_evaluation else 'Failure!')
    print('Top measurement:', result.top_measurement)

    assert result.oracle_evaluation is True
    assert result.top_measurement == '11'

    # define an oracle with Statevector
    oracle = Statevector.from_label('11')
    problem = AmplificationProblem(oracle, is_good_state=good_state)
    grover = Grover(sampler=sampler)
    result = grover.amplify(problem)
    print('Result type:', type(result))
    print('Success!' if result.oracle_evaluation else 'Failure!')
    print('Top measurement:', result.top_measurement)

    assert result.oracle_evaluation is True
    assert result.top_measurement == '11'

    # define an oracle with PhaseOracle
    expression = '(a & b)'
    try:
        oracle = PhaseOracle(expression)
        problem = AmplificationProblem(oracle)
        grover = Grover(sampler=sampler)
        result = grover.amplify(problem)
        print('Result type:', type(result))
        print('Success!' if result.oracle_evaluation else 'Failure!')
        print('Top measurement:', result.top_measurement)

        assert result.oracle_evaluation is True
        assert result.top_measurement == '11'
    except MissingOptionalLibraryError as ex:
        print(ex)
        pytest.fail('Error')


def test_amplitude_amplification():
    # Specifying `state_preparation`
    # to prepare a superposition of |01>, |10>, and |11>
    oracle = QuantumCircuit(3)
    oracle.ccz(0, 1, 2)

    theta = 2 * np.arccos(1 / np.sqrt(3))
    state_preparation = QuantumCircuit(3)
    state_preparation.ry(theta, 0)
    state_preparation.ch(0, 1)
    state_preparation.x(1)
    state_preparation.h(2)

    sampler = Sampler(run_options={"shots":1024,"seed_simulator":80},transpile_options={"seed_transpiler":50})

    # we only care about the first two bits being in state 1, thus add both possibilities for the last qubit
    problem = AmplificationProblem(oracle, state_preparation=state_preparation, is_good_state=['110', '111'])

    grover = Grover(sampler=sampler)
    result = grover.amplify(problem)
    print('Success!' if result.oracle_evaluation else 'Failure!')
    print('Top measurement:', result.top_measurement)

    assert result.oracle_evaluation is True
    assert result.top_measurement == '110' or result.top_measurement == '111'


def test_grover_operator():
    oracle = QuantumCircuit(5)
    oracle.ccz(0,1,2)
    grover_op = GroverOperator(oracle, reflection_qubits=[0, 1, 2], insert_barriers=True)

    sampler = Sampler(run_options={"shots":1024,"seed_simulator":80},transpile_options={"seed_transpiler":50})

    problem = AmplificationProblem(oracle=oracle, grover_operator=grover_op,
                                   is_good_state=['00111', '01111', '10111', '11111'])
    grover = Grover(sampler=sampler)
    result = grover.amplify(problem)
    print('Success!' if result.oracle_evaluation else 'Failure!')
    print('Top measurement:', result.top_measurement)

    amplified_values = list(sorted(map(lambda kv: kv[0],
                                       filter(lambda kv: kv[1] >= 150/1024, result.circuit_results[0].items()))))

    assert result.oracle_evaluation is True
    assert amplified_values == ['00111', '01111', '10111', '11111']
