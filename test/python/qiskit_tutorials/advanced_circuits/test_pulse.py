# The code in this file has been written using part of the code in the Qiskit tutorial below.
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits_advanced/06_building_pulse_schedules.ipynb  # noqa: E501
# https://github.com/Qiskit/qiskit-tutorials/blob/2393ed1c47952ecf7802f478d99a0a2c105ea9b1/tutorials/circuits_advanced/07_pulse_scheduler.ipynb  # noqa: E501

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
from qiskit import QuantumCircuit, pulse, schedule
from qiskit.pulse import DriveChannel

from qdd import QddProvider
from test.python.helpers.circuit_helper import assert_job_failed


def test_pulse_backend_configuration():
    """Tests that Pulse is not supported in QddBackend."""

    backend = QddProvider().get_backend()
    assert backend.configuration().open_pulse is False


def test_schedule():
    """Tests schedule-related backend behavior."""

    backend = QddProvider().get_backend()
    channel = DriveChannel(0)
    with pulse.build(backend) as delay_5dt:
        pulse.delay(5, channel)

    job = backend.run(delay_5dt, seed_simulator=80)
    assert_job_failed(job)

    circ = QuantumCircuit(2, 2)
    circ.h(0)
    circ.cx(0, 1)
    circ.measure([0, 1], [0, 1])

    with pulse.build() as h_q0:
        pulse.play(
            pulse.library.Gaussian(duration=256, amp=0.2, sigma=50, name="custom"),
            pulse.DriveChannel(0),
        )

    circ.add_calibration("h", [0], h_q0)
    with pytest.raises(Exception):
        schedule(circ, backend)
