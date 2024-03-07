# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/d0d8e42afa0e42f6742b06bfb58e267cb9137599/tutorials/algorithms/07_grover_examples.ipynb  # noqa: E501

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import os
import tempfile

import pytest
from qiskit import MissingOptionalLibraryError
from qiskit_algorithms import AmplificationProblem, Grover
from qiskit.circuit.library import PhaseOracle

from qdd.qdd_sampler_like_aer import Sampler


def test_sat():
    input_3sat_instance = '''
    c example DIMACS-CNF 3-SAT
    p cnf 3 5
    -1 -2 -3 0
    1 -2 3 0
    1 2 -3 0
    1 -2 -3 0
    -1 2 3 0
    '''

    fp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
    fp.write(input_3sat_instance)
    file_name = fp.name
    fp.close()
    oracle = None
    is_error = False
    try:
        oracle = PhaseOracle.from_dimacs_file(file_name)
    except MissingOptionalLibraryError as ex:
        print(ex)
        is_error = True
    finally:
        os.remove(file_name)

    if is_error:
        pytest.fail('Error')

    problem = AmplificationProblem(oracle, is_good_state=oracle.evaluate_bitstring)
    grover = Grover(sampler=Sampler(run_options={"shots":1024,"seed_simulator":80},transpile_options={"seed_transpiler":50}))
    result = grover.amplify(problem)
    print(result.assignment)
    assert result.assignment == '000' or result.assignment == '011' or result.assignment == '101'

    amplified_values = list(sorted(map(lambda kv: kv[0],
                                       filter(lambda kv: kv[1] >= 150/1024, result.circuit_results[0].items()))))
    assert amplified_values == ['000', '011', '101']


def test_grover_oracle_boolean_exp():
    expression = '(w ^ x) & ~(y ^ z) & (x & y & z)'
    try:
        oracle = PhaseOracle(expression)
        problem = AmplificationProblem(oracle, is_good_state=oracle.evaluate_bitstring)
        grover = Grover(sampler=Sampler(run_options={"shots":1024,"seed_simulator":80},transpile_options={"seed_transpiler":50}))
        result = grover.amplify(problem)
        assert result.top_measurement == '1110'
    except MissingOptionalLibraryError as ex:
        print(ex)
        pytest.fail('Error')
