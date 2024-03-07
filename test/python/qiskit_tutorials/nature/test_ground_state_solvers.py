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
from qiskit.circuit.library import TwoLocal
from qiskit_algorithms import VQE, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms.utils import algorithm_globals
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.units import DistanceUnit

from qdd.qdd_estimator_like_aer import Estimator


def test_solver():
    seed = 50
    algorithm_globals.random_seed = seed

    driver = PySCFDriver(
        atom="H 0 0 0; H 0 0 0.735",
        basis="sto3g",
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM,
    )
    es_problem = driver.run()
    mapper = JordanWignerMapper()


    ansatz = UCCSD(
        es_problem.num_spatial_orbitals,
        es_problem.num_particles,
        mapper,
        initial_state=HartreeFock(
            es_problem.num_spatial_orbitals,
            es_problem.num_particles,
            mapper,
        ),
    )


    # numpy
    numpy_solver = NumPyMinimumEigensolver()
    calc = GroundStateEigensolver(mapper, numpy_solver)
    res_numpy = calc.solve(es_problem)
    print(res_numpy)

    # VQE
    optimizer = COBYLA(maxiter=100)
    vqe_solver = VQE(estimator=Estimator(), ansatz=ansatz, optimizer=optimizer)
    calc = GroundStateEigensolver(mapper, vqe_solver)
    res_vqe = calc.solve(es_problem)
    print(res_vqe)
    assert res_vqe.groundenergy == pytest.approx(res_numpy.groundenergy, abs=0.2)

    # VQE2
    tl_circuit = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
    optimizer = COBYLA(maxiter=100)
    another_solver = VQE(estimator=Estimator(), ansatz=tl_circuit, optimizer=optimizer)
    calc = GroundStateEigensolver(mapper, another_solver)
    res_another = calc.solve(es_problem)
    print(res_another)
    assert res_another.groundenergy == pytest.approx(res_numpy.groundenergy, abs=0.2)

    # The section "Using a filter function" is skipped, because it does not include Quantum calculation.
