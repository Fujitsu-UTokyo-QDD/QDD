# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits_advanced/04_transpiler_passes_and_passmanager.ipynb  # noqa: E501

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

import math
from logging import DEBUG, INFO

from qiskit import QuantumCircuit, transpile

from qdd import QddProvider


def test_preset_pass_manager():
    """Tests the behavior of preset pass managers in the transpiler."""

    qc = QuantumCircuit(10)

    random_state = [
        1 / math.sqrt(4) * complex(0, 1),
        1 / math.sqrt(8) * complex(1, 0),
        0,
        0,
        0,
        0,
        0,
        0,
        1 / math.sqrt(8) * complex(1, 0),
        1 / math.sqrt(8) * complex(0, 1),
        0,
        0,
        0,
        0,
        1 / math.sqrt(4) * complex(1, 0),
        1 / math.sqrt(8) * complex(1, 0),
    ]

    qc.initialize(random_state, range(4))

    backend = QddProvider().get_backend()
    optimized_0 = transpile(
        qc, backend=backend, seed_transpiler=11, optimization_level=0
    )
    print("gates = ", optimized_0.count_ops())
    print("depth = ", optimized_0.depth())

    optimized_1 = transpile(
        qc, backend=backend, seed_transpiler=11, optimization_level=1
    )
    print("gates = ", optimized_1.count_ops())
    print("depth = ", optimized_1.depth())

    optimized_2 = transpile(
        qc, backend=backend, seed_transpiler=11, optimization_level=2
    )
    print("gates = ", optimized_2.count_ops())
    print("depth = ", optimized_2.depth())

    optimized_3 = transpile(
        qc, backend=backend, seed_transpiler=11, optimization_level=3
    )
    print("gates = ", optimized_3.count_ops())
    print("depth = ", optimized_3.depth())

    # In this test, we would like to check whether the above transpilations finish with no errors raised.
    # So, here, there are no specific properties to assert.
    assert True


def test_transpiler_log(caplog):
    """Tests whether we can see logger messages from the transpiler."""

    caplog.set_level(DEBUG)

    log_circ = QuantumCircuit(2, 2)
    log_circ.h(0)
    log_circ.h(1)
    log_circ.h(1)
    log_circ.x(1)
    log_circ.cx(0, 1)
    log_circ.measure([0, 1], [0, 1])

    backend = QddProvider().get_backend()
    transpile(log_circ, backend, seed_transpiler=50)

    assert any(t[0].startswith("qiskit.transpiler.") for t in caplog.record_tuples)
    assert any(t[1] == DEBUG for t in caplog.record_tuples)
    assert any(t[1] == INFO for t in caplog.record_tuples)
