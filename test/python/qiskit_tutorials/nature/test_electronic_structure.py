# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-nature/blob/0.2.2/docs/tutorials/01_electronic_structure.ipynb

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

from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.drivers import Molecule
from qiskit_nature.drivers.second_quantization import ElectronicStructureDriverType, ElectronicStructureMoleculeDriver
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper
from qiskit_nature.problems.second_quantization import ElectronicStructureProblem


def test_electronic_structure():
    molecule = Molecule(
        geometry=[["H", [0.0, 0.0, 0.0]], ["H", [0.0, 0.0, 0.735]]], charge=0, multiplicity=1
    )
    driver = ElectronicStructureMoleculeDriver(
        molecule, basis="sto3g", driver_type=ElectronicStructureDriverType.PYSCF
    )

    # setup Hamiltonian
    es_problem = ElectronicStructureProblem(driver)
    second_q_op = es_problem.second_q_ops()
    print(second_q_op[0])

    # transform Hamiltonian
    qubit_converter = QubitConverter(mapper=JordanWignerMapper())
    qubit_op = qubit_converter.convert(second_q_op[0])
    print(qubit_op)
    assert qubit_op.num_qubits == 4

    # reduce the number of qubits
    qubit_converter = QubitConverter(mapper=ParityMapper(), two_qubit_reduction=True)
    qubit_op = qubit_converter.convert(second_q_op[0], num_particles=es_problem.num_particles)
    print(qubit_op)
    assert qubit_op.num_qubits == 2
