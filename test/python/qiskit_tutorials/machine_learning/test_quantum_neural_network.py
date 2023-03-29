# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-machine-learning/blob/stable/0.3/docs/tutorials/01_neural_networks.ipynb

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
from qiskit import Aer, QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.circuit.library import RealAmplitudes, ZZFeatureMap
from qiskit.opflow import AerPauliExpectation, Gradient, ListOp, PauliExpectation, PauliSumOp, StateFn
from qiskit.utils import QuantumInstance, algorithm_globals

from qdd import QddProvider


@pytest.mark.filterwarnings('ignore:`np.long` is a deprecated alias for `np.compat.long`:DeprecationWarning:numba.*')
@pytest.mark.filterwarnings('ignore:`np.bool` is a deprecated alias for the builtin `bool`:DeprecationWarning:numba.*')
@pytest.mark.filterwarnings('ignore:`np.MachAr` is deprecated:DeprecationWarning:numba.*')
class TestQNN:

    # OpflowQNN
    def test_opflow_qnn(self):
        # This import statement emits DeprecationWarning in namba.*
        # To apply filterwarnings directive, the import statement must be within a test method instead of at the header.
        from qiskit_machine_learning.neural_networks import OpflowQNN

        algorithm_globals.random_seed = 42
        # define quantum instances (statevector and sample based)
        qi_sv = QuantumInstance(Aer.get_backend('aer_simulator_statevector'), seed_transpiler=50, seed_simulator=80)
        # qdd
        qi_qdd = QuantumInstance(QddProvider().get_backend(), seed_transpiler=50, seed_simulator=80)

        # construct parametrized circuit
        params1 = [Parameter('input1'), Parameter('weight1')]
        qc1 = QuantumCircuit(1)
        qc1.h(0)
        qc1.ry(params1[0], 0)
        qc1.rx(params1[1], 0)
        qc_sfn1 = StateFn(qc1)

        # construct cost operator
        h1 = StateFn(PauliSumOp.from_list([('Z', 1.0), ('X', 1.0)]))

        # combine operator and circuit to objective function
        op1 = ~h1 @ qc_sfn1

        # construct OpflowQNN with the operator, the input parameters, the weight parameters,
        # the expected value, gradient, and quantum instance.
        qnn1_aer = OpflowQNN(op1, [params1[0]], [params1[1]], AerPauliExpectation(), Gradient(), qi_sv)
        qnn1_qdd = OpflowQNN(op1, [params1[0]], [params1[1]], PauliExpectation(), Gradient(), qi_qdd)
        # define (random) input and weights
        input1 = algorithm_globals.random.random(qnn1_aer.num_inputs)
        weights1 = algorithm_globals.random.random(qnn1_aer.num_weights)

        # QNN batched forward pass
        out_forward1_aer = qnn1_aer.forward(input1, weights1)
        out_forward1_qdd = qnn1_qdd.forward(input1, weights1)
        # QNN batched backward pass
        out_backward1_aer = qnn1_aer.backward(input1, weights1)
        out_backward1_qdd = qnn1_qdd.backward(input1, weights1)
        assert np.allclose(out_forward1_aer, out_forward1_qdd, atol=0.1)
        assert out_backward1_aer[0] is None and out_backward1_qdd[0] is None
        assert np.allclose(out_backward1_aer[1], out_backward1_qdd[1], atol=0.1)

        # another qnn
        op2 = ListOp([op1, op1])
        qnn2_aer = OpflowQNN(op2, [params1[0]], [params1[1]], AerPauliExpectation(), Gradient(), qi_sv)
        qnn2_qdd = OpflowQNN(op2, [params1[0]], [params1[1]], PauliExpectation(), Gradient(), qi_qdd)
        # QNN forward pass
        out_forward2_aer = qnn2_aer.forward(input1, weights1)
        out_forward2_qdd = qnn2_qdd.forward(input1, weights1)
        # QNN backward pass
        out_backward2_aer = qnn2_aer.backward(input1, weights1)
        out_backward2_qdd = qnn2_qdd.backward(input1, weights1)
        assert np.allclose(out_forward2_aer, out_forward2_qdd, atol=0.1)
        assert out_backward2_aer[0] is None and out_backward2_qdd[0] is None
        assert np.allclose(out_backward2_aer[1], out_backward2_qdd[1], atol=0.1)

    def test_two_layer_qnn(self):
        # This import statement emits DeprecationWarning in namba.*
        # To apply filterwarnings directive, the import statement must be within a test method instead of at the header.
        from qiskit_machine_learning.neural_networks import TwoLayerQNN

        algorithm_globals.random_seed = 42
        # specify the number of qubits
        num_qubits = 3
        # define quantum instances (statevector and sample based)
        qi_sv = QuantumInstance(Aer.get_backend('aer_simulator_statevector'), seed_transpiler=50, seed_simulator=80)
        # qdd
        qi_qdd = QuantumInstance(QddProvider().get_backend(), seed_transpiler=50, seed_simulator=80)

        # specify the feature map
        fm = ZZFeatureMap(num_qubits, reps=2)
        # specify the ansatz
        ansatz = RealAmplitudes(num_qubits, reps=1)
        # specify the observable
        observable = PauliSumOp.from_list([("Z" * num_qubits, 1)])
        # define two layer QNN
        qnn3_aer = TwoLayerQNN(num_qubits, feature_map=fm, ansatz=ansatz, observable=observable, quantum_instance=qi_sv)
        qnn3_qdd = TwoLayerQNN(num_qubits, feature_map=fm, ansatz=ansatz, observable=observable,
                                  quantum_instance=qi_qdd)
        # define (random) input and weights
        input3 = algorithm_globals.random.random(qnn3_aer.num_inputs)
        weights3 = algorithm_globals.random.random(qnn3_aer.num_weights)
        # QNN forward pass
        out_forward3_aer = qnn3_aer.forward(input3, weights3)
        out_forward3_qdd = qnn3_qdd.forward(input3, weights3)
        # QNN backward pass
        out_backward3_aer = qnn3_aer.backward(input3, weights3)
        out_backward3_qdd = qnn3_qdd.backward(input3, weights3)
        assert np.allclose(out_forward3_aer, out_forward3_qdd, atol=0.1)
        assert out_backward3_aer[0] is None and out_backward3_qdd[0] is None
        assert np.allclose(out_backward3_aer[1], out_backward3_qdd[1], atol=0.1)

    def test_circuit_qnn(self):
        # This import statement emits DeprecationWarning in namba.*
        # To apply filterwarnings directive, the import statement must be within a test method instead of at the header.
        from qiskit_machine_learning.neural_networks import CircuitQNN

        algorithm_globals.random_seed = 42
        # specify the number of qubits
        num_qubits = 3
        # qiskit
        # To compare with qdd, shots is 1024 instead of original 10.
        n_shots = 1024
        qi_qasm = QuantumInstance(Aer.get_backend('aer_simulator'), shots=n_shots, seed_transpiler=50,
                                  seed_simulator=80)
        # qdd
        qi_qdd = QuantumInstance(QddProvider().get_backend(), seed_transpiler=50, seed_simulator=80)

        # 4.1
        qc = RealAmplitudes(num_qubits, entanglement="linear", reps=1)
        # specify circuit QNN
        qnn4_aer = CircuitQNN(qc, [], qc.parameters, sparse=True, quantum_instance=qi_qasm)
        qnn4_qdd = CircuitQNN(qc, [], qc.parameters, sparse=True, quantum_instance=qi_qdd)
        # define (random) input and weights
        input4 = algorithm_globals.random.random(qnn4_aer.num_inputs)
        weights4 = algorithm_globals.random.random(qnn4_aer.num_weights)
        # QNN forward pass
        out_forward4_aer = qnn4_aer.forward(input4, weights4).todense()
        out_forward4_qdd = qnn4_qdd.forward(input4, weights4).todense()
        # QNN backward pass, returns a tuple of sparse matrices
        out_backward4_aer = qnn4_aer.backward(input4, weights4)
        out_backward4_qdd = qnn4_qdd.backward(input4, weights4)
        assert np.allclose(out_forward4_aer, out_forward4_qdd, atol=0.1)
        assert out_backward4_aer[0] is None and out_backward4_qdd[0] is None
        assert np.allclose(out_backward4_aer[1].todense(), out_backward4_qdd[1].todense(), atol=0.1)

        # 4.2
        def parity(x): return '{:b}'.format(x).count('1') % 2
        output_shape = 2  # this is required in case of a callable with dense output
        qnn6_aer = CircuitQNN(qc, [], qc.parameters, sparse=False, interpret=parity, output_shape=output_shape,
                              quantum_instance=qi_qasm)
        qnn6_qdd = CircuitQNN(qc, [], qc.parameters, sparse=False, interpret=parity, output_shape=output_shape,
                                 quantum_instance=qi_qdd)
        # define (random) input and weights
        input6 = algorithm_globals.random.random(qnn6_aer.num_inputs)
        weights6 = algorithm_globals.random.random(qnn6_aer.num_weights)
        # QNN forward pass
        out_forward6_aer = qnn6_aer.forward(input6, weights6)
        out_forward6_qdd = qnn6_qdd.forward(input6, weights6)
        # QNN backward pass
        out_backward6_aer = qnn6_aer.backward(input6, weights6)
        out_backward6_qdd = qnn6_qdd.backward(input6, weights6)
        assert np.allclose(out_forward6_aer, out_forward6_qdd, atol=0.1)
        assert out_backward6_aer[0] is None and out_backward6_qdd[0] is None
        assert np.allclose(out_backward6_aer[1], out_backward6_qdd[1], atol=0.1)

        # 4.3
        # specify circuit QNN
        qnn7_aer = CircuitQNN(qc, [], qc.parameters, sampling=True, quantum_instance=qi_qasm)
        qnn7_qdd = CircuitQNN(qc, [], qc.parameters, sampling=True, quantum_instance=qi_qdd)
        # define (random) input and weights
        input7 = algorithm_globals.random.random(qnn7_aer.num_inputs)
        weights7 = algorithm_globals.random.random(qnn7_aer.num_weights)
        # QNN forward pass, results in samples of measured bit strings mapped to integers
        out_forward7_aer = qnn7_aer.forward(input7, weights7)
        out_forward7_qdd = qnn7_qdd.forward(input7, weights7)
        # QNN backward pass
        out_backward7_aer = qnn7_aer.backward(input7, weights7)
        out_backward7_qdd = qnn7_qdd.backward(input7, weights7)

        # compare the sampled probability for each value in possible range of 0--2^(#qubits)
        for i in range(0, 2**num_qubits):
            assert np.count_nonzero(out_forward7_aer == i) / float(n_shots)\
                == pytest.approx(np.count_nonzero(out_forward7_qdd == i) / float(n_shots), abs=0.1)
        assert out_backward7_aer == out_backward7_qdd

        # 4.4
        # specify circuit QNN
        qnn8_aer = CircuitQNN(qc, [], qc.parameters, sampling=True, interpret=parity, quantum_instance=qi_qasm)
        qnn8_qdd = CircuitQNN(qc, [], qc.parameters, sampling=True, interpret=parity, quantum_instance=qi_qdd)
        # define (random) input and weights
        input8 = algorithm_globals.random.random(qnn8_aer.num_inputs)
        weights8 = algorithm_globals.random.random(qnn8_aer.num_weights)
        # QNN forward pass, results in samples of measured bit strings
        out_forward8_aer = qnn8_aer.forward(input8, weights8)
        out_forward8_qdd = qnn8_qdd.forward(input8, weights8)
        # QNN backward pass
        out_backward8_aer = qnn8_aer.backward(input8, weights8)
        out_backward8_qdd = qnn8_qdd.backward(input8, weights8)

        # compare the sampled probability for each parity value (0 or 1)
        assert np.count_nonzero(out_forward8_aer == 0) / float(n_shots)\
            == pytest.approx(np.count_nonzero(out_forward8_qdd == 0) / float(n_shots), abs=0.1)
        assert np.count_nonzero(out_forward8_aer == 1) / float(n_shots)\
            == pytest.approx(np.count_nonzero(out_forward8_qdd == 1) / float(n_shots), abs=0.1)
        assert out_backward8_aer == out_backward8_qdd
