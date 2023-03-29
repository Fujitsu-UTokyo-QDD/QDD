# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-machine-learning/blob/stable/0.3/docs/tutorials/02_neural_network_classifier_and_regressor.ipynb

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

import numpy as np
import pytest
from qiskit import QuantumCircuit
from qiskit.algorithms.optimizers import COBYLA, L_BFGS_B
from qiskit.circuit import Parameter
from qiskit.circuit.library import RealAmplitudes, ZZFeatureMap
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit_machine_learning.algorithms.classifiers import VQC, NeuralNetworkClassifier
from qiskit_machine_learning.algorithms.regressors import VQR, NeuralNetworkRegressor
from qiskit_machine_learning.neural_networks import CircuitQNN, TwoLayerQNN

from qdd import QddProvider


@pytest.mark.slow
class TestClassifierRegressor:
    algorithm_globals.random_seed = 42
    quantum_instance = QuantumInstance(QddProvider().get_backend(), seed_transpiler=50, seed_simulator=80)

    def test_classification(self):
        num_inputs = 2
        num_samples = 20
        X = 2 * algorithm_globals.random.random([num_samples, num_inputs]) - 1
        y01 = 1 * (np.sum(X, axis=1) >= 0)  # in { 0,  1}
        y = 2 * y01 - 1  # in {-1, +1}
        y_one_hot = np.zeros((num_samples, 2))
        for i in range(num_samples):
            y_one_hot[i, y01[i]] = 1

        # OpflowQNN
        # construct QNN
        opflow_qnn = TwoLayerQNN(num_inputs, quantum_instance=self.quantum_instance)
        # QNN maps inputs to [-1, +1]
        opflow_qnn.forward(X[0, :], algorithm_globals.random.random(opflow_qnn.num_weights))
        # construct neural network classifier
        opflow_classifier = NeuralNetworkClassifier(opflow_qnn, optimizer=COBYLA())
        # fit classifier to data
        opflow_classifier.fit(X, y)
        # score classifier
        opflow_classifier.score(X, y)
        # evaluate data points
        y_predict = opflow_classifier.predict(X)
        wrong = 0
        for x, y_target, y_p in zip(X, y, y_predict):
            if y_target != y_p:
                wrong += 1
        assert wrong < 10

        # CircuitQNN
        # construct feature map
        feature_map = ZZFeatureMap(num_inputs)
        # construct ansatz
        ansatz = RealAmplitudes(num_inputs, reps=1)
        # construct quantum circuit
        qc = QuantumCircuit(num_inputs)
        qc.append(feature_map, range(num_inputs))
        qc.append(ansatz, range(num_inputs))

        def parity(value):
            return "{:b}".format(value).count("1") % 2
        output_shape = 2  # corresponds to the number of classes, possible outcomes of the (parity) mapping.
        # construct QNN
        circuit_qnn = CircuitQNN(
            circuit=qc,
            input_params=feature_map.parameters,
            weight_params=ansatz.parameters,
            interpret=parity,
            output_shape=output_shape,
            quantum_instance=self.quantum_instance,
        )
        # construct classifier
        circuit_classifier = NeuralNetworkClassifier(neural_network=circuit_qnn, optimizer=COBYLA())
        # fit classifier to data
        circuit_classifier.fit(X, y01)
        # score classifier
        circuit_classifier.score(X, y01)
        # evaluate data points
        y_predict = circuit_classifier.predict(X)
        wrong = 0
        for x, y_target, y_p in zip(X, y01, y_predict):
            if y_target != y_p:
                wrong += 1
        assert wrong < 10

        # VQC
        # construct feature map, ansatz, and optimizer
        feature_map = ZZFeatureMap(num_inputs)
        ansatz = RealAmplitudes(num_inputs, reps=1)
        # construct variational quantum classifier
        vqc = VQC(
            feature_map=feature_map,
            ansatz=ansatz,
            loss="cross_entropy",
            optimizer=COBYLA(),
            quantum_instance=self.quantum_instance
        )
        # fit classifier to data
        vqc.fit(X, y_one_hot)
        # score classifier
        vqc.score(X, y_one_hot)
        # evaluate data points
        y_predict = vqc.predict(X)
        wrong = 0
        for x, y_target, y_p in zip(X, y_one_hot, y_predict):
            if not np.all(y_target == y_p):
                wrong += 1
        assert wrong < 15

    @pytest.mark.filterwarnings("ignore:The loss argument value \"l2\" is deprecated:DeprecationWarning")
    def test_regression(self):
        num_samples = 20
        eps = 0.2
        lb, ub = -np.pi, np.pi
        X_ = np.linspace(lb, ub, num=50).reshape(50, 1)
        def f(x): return np.sin(x)
        X = (ub - lb) * algorithm_globals.random.random([num_samples, 1]) + lb
        y = f(X[:, 0]) + eps * (2 * algorithm_globals.random.random(num_samples) - 1)

        # OpflowQNN
        # construct simple feature map
        param_x = Parameter("x")
        feature_map = QuantumCircuit(1, name="fm")
        feature_map.ry(param_x, 0)
        # construct simple ansatz
        param_y = Parameter("y")
        ansatz = QuantumCircuit(1, name="vf")
        ansatz.ry(param_y, 0)
        # construct QNN
        regression_opflow_qnn = TwoLayerQNN(1, feature_map, ansatz, quantum_instance=self.quantum_instance)
        # construct the regressor from the neural network
        regressor = NeuralNetworkRegressor(
            neural_network=regression_opflow_qnn, loss="l2", optimizer=L_BFGS_B()
        )
        # fit to data
        regressor.fit(X, y)
        # score the result
        regressor.score(X, y)
        # evaluation
        y_ = regressor.predict(X_)
        for x, y_p in zip(X_, y_):
            assert y_p == pytest.approx(f(x), abs=0.2)

        # VQR
        vqr = VQR(
            feature_map=feature_map,
            ansatz=ansatz,
            optimizer=L_BFGS_B(),
            quantum_instance=self.quantum_instance
        )
        # fit regressor
        vqr.fit(X, y)
        # score result
        vqr.score(X, y)
        # evaluation
        y_ = vqr.predict(X_)
        for x, y_p in zip(X_, y_):
            assert y_p == pytest.approx(f(x), abs=0.2)
