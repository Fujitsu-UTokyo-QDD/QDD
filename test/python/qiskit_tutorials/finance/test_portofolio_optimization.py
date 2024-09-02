# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://qiskit-community.github.io/qiskit-finance/tutorials/01_portfolio_optimization.html

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

import datetime

import numpy as np
from qiskit_algorithms import QAOA, NumPyMinimumEigensolver, SamplingVQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit_algorithms.utils import algorithm_globals
from qiskit_finance.applications.optimization import PortfolioOptimization
from qiskit_finance.data_providers import RandomDataProvider
from qiskit_optimization.algorithms import MinimumEigenOptimizer

from qdd.qdd_sampler import Sampler


class TestPortofolioOptimization:

    def test_portofolio_optimization(self):
        # set number of assets (= number of qubits)
        num_assets = 4
        seed = 123

        # Generate expected return and covariance matrix from (random) time-series
        stocks = [("TICKER%s" % i) for i in range(num_assets)]
        data = RandomDataProvider(
            tickers=stocks,
            start=datetime.datetime(2016, 1, 1),
            end=datetime.datetime(2016, 1, 30),
            seed=seed,
        )
        data.run()
        mu = data.get_period_return_mean_vector()
        sigma = data.get_period_return_covariance_matrix()

        q = 0.5  # set risk factor
        budget = num_assets // 2  # set budget

        portfolio = PortfolioOptimization(
            expected_returns=mu, covariances=sigma, risk_factor=q, budget=budget
        )
        qp = portfolio.to_quadratic_program()

        # Numpy
        exact_mes = NumPyMinimumEigensolver()
        exact_eigensolver = MinimumEigenOptimizer(exact_mes)
        result_numpy = exact_eigensolver.solve(qp)

        # VQE
        algorithm_globals.random_seed = 1234
        cobyla = COBYLA()
        cobyla.set_options(maxiter=500)
        ry = TwoLocal(num_assets, "ry", "cz", reps=3, entanglement="full")
        vqe_mes = SamplingVQE(sampler=Sampler(), ansatz=ry, optimizer=cobyla)
        vqe = MinimumEigenOptimizer(vqe_mes)
        result_vqe = vqe.solve(qp)
        assert np.all(result_vqe.x == result_numpy.x)

        # QAOA
        cobyla = COBYLA()
        cobyla.set_options(maxiter=250)
        qaoa_mes = QAOA(sampler=Sampler(), optimizer=cobyla, reps=3)
        qaoa = MinimumEigenOptimizer(qaoa_mes)
        result_qaoa = qaoa.solve(qp)
        assert np.all(result_qaoa.x == result_numpy.x)
