# The code in this file has been written using part of the code in the Qiskit tutorial:
# https://learning.quantum.ibm.com/tutorial/variational-quantum-eigensolver

import numpy as np
import pytest
from qiskit.circuit.library import EfficientSU2
from qiskit.primitives import StatevectorEstimator as QiskitEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from scipy.optimize import minimize

from qdd import QddProvider
from qdd.qdd_estimator import Estimator as QddEstimator


def test_vqe():
    backend = QddProvider().get_backend()
    estimator_qdd = QddEstimator(run_options={"shots": None}, approximation=True)

    estimator_qiskit = QiskitEstimator()

    hamiltonian = SparsePauliOp.from_list(
        [("YZ", 0.3980), ("ZI", -0.3980), ("ZZ", -0.0113), ("XX", 0.1810)]
    )

    ansatz = EfficientSU2(hamiltonian.num_qubits)
    num_params = ansatz.num_parameters

    pm = generate_preset_pass_manager(backend=backend, optimization_level=3)

    ansatz_isa = pm.run(ansatz)

    def cost_func_qdd(params, ansatz, hamiltonian, estimator):
        result = estimator.run(ansatz, hamiltonian, params).result()
        energy = result.values[0]

        return energy

    x0 = 2 * np.pi * np.random.random(num_params)
    result_qdd = minimize(
        cost_func_qdd,
        x0,
        args=(ansatz_isa, hamiltonian, estimator_qdd),
        method="cobyla",
    )
    energy_qdd = result_qdd.fun

    def cost_func_qiskit(params, ansatz, hamiltonian, estimator):
        result = estimator.run(pubs=[(ansatz, hamiltonian, params)]).result()
        energy = result[0].data.evs

        return energy

    x0 = 2 * np.pi * np.random.random(num_params)
    result_qiskit = minimize(
        cost_func_qiskit,
        x0,
        args=(ansatz_isa, hamiltonian, estimator_qiskit),
        method="cobyla",
    )
    energy_qiskit = result_qiskit.fun

    assert pytest.approx(energy_qdd, abs=1e-2) == energy_qiskit
