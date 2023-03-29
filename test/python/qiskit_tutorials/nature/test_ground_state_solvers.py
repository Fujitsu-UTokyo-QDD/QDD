# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-nature/blob/0.2.2/docs/tutorials/03_ground_state_solvers.ipynb

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

import pytest
from qiskit.algorithms import VQE, NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit_nature.algorithms import GroundStateEigensolver, VQEUCCFactory
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.drivers import Molecule
from qiskit_nature.drivers.second_quantization import ElectronicStructureDriverType, ElectronicStructureMoleculeDriver
from qiskit_nature.mappers.second_quantization import JordanWignerMapper
from qiskit_nature.problems.second_quantization import ElectronicStructureProblem

from qdd import QddProvider


@pytest.mark.slow
def test_solver():
    seed = 50
    algorithm_globals.random_seed = seed

    molecule = Molecule(geometry=[['H', [0., 0., 0.]],
                                  ['H', [0., 0., 0.735]]],
                        charge=0, multiplicity=1)
    driver = ElectronicStructureMoleculeDriver(molecule, basis='sto3g', driver_type=ElectronicStructureDriverType.PYSCF)
    es_problem = ElectronicStructureProblem(driver)
    qubit_converter = QubitConverter(JordanWignerMapper())

    # numpy
    numpy_solver = NumPyMinimumEigensolver()
    calc = GroundStateEigensolver(qubit_converter, numpy_solver)
    res_numpy = calc.solve(es_problem)
    print(res_numpy)

    # VQE
    quantum_instance = QuantumInstance(QddProvider().get_backend(), seed_transpiler=seed, seed_simulator=seed)
    optimizer = COBYLA(maxiter=100)
    vqe_solver = VQEUCCFactory(quantum_instance, optimizer=optimizer)
    calc = GroundStateEigensolver(qubit_converter, vqe_solver)
    res_vqe = calc.solve(es_problem)
    print(res_vqe)
    assert res_vqe.electronic_energies == pytest.approx(res_numpy.electronic_energies, abs=0.2)

    # VQE2
    tl_circuit = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
    optimizer = COBYLA(maxiter=100)
    another_solver = VQE(ansatz=tl_circuit, quantum_instance=quantum_instance, optimizer=optimizer)
    calc = GroundStateEigensolver(qubit_converter, another_solver)
    res_another = calc.solve(es_problem)
    print(res_another)
    assert res_another.electronic_energies == pytest.approx(res_numpy.electronic_energies, abs=0.2)

    # The section "Using a filter function" is skipped, because it does not include Quantum calculation.
