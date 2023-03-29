# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-nature/blob/0.2.2/docs/tutorials/02_vibrational_structure.ipynb

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

import os.path

from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.drivers.second_quantization import GaussianForcesDriver
from qiskit_nature.mappers.second_quantization import DirectMapper
from qiskit_nature.problems.second_quantization import VibrationalStructureProblem


def test_vibrational_structure():
    driver = GaussianForcesDriver(logfile=os.path.dirname(__file__) + '/aux_files/CO2_freq_B3LYP_ccpVDZ.log')

    # setup Hamiltonian
    vibrational_problem = VibrationalStructureProblem(driver, num_modals=2, truncation_order=2)
    second_q_ops = vibrational_problem.second_q_ops()
    print(second_q_ops[0])

    # transform Hamiltonian
    qubit_converter = QubitConverter(mapper=DirectMapper())
    qubit_op = qubit_converter.convert(second_q_ops[0])
    print(qubit_op)
    assert qubit_op.num_qubits == 8

    # try different number of modals
    vibrational_problem = VibrationalStructureProblem(driver, num_modals=3, truncation_order=2)
    second_q_ops = vibrational_problem.second_q_ops()
    qubit_converter = QubitConverter(mapper=DirectMapper())
    qubit_op = qubit_converter.convert(second_q_ops[0])
    print(qubit_op)
    print(type(qubit_op))
    assert qubit_op.num_qubits == 12
