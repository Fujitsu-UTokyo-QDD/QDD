# This code is part of Qiskit.
#
# (C) Copyright IBM 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# Code adapted from
# https://github.com/Qiskit/qiskit-aer/blob/main/qiskit_aer/primitives/sampler.py

"""
Sampler class.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Iterable

import numpy as np
from qiskit.compiler import transpile
from qiskit.exceptions import QiskitError
from qiskit.primitives import SamplerResult
from qiskit.primitives import (
    BaseSamplerV2,
    SamplerPubResult,
    SamplerPubLike,
    DataBin,
    BitArray,
)
from qiskit.primitives.containers.sampler_pub import SamplerPub
from qiskit.result import QuasiDistribution

from qdd import QddProvider
from .qdd_utils import _circuit_key


class Sampler(BaseSamplerV2):
    """
    QDD implementation of SamplerV2 class.
    """

    def __init__(
        self,
        *,
        backend_options: dict | None = None,
        transpile_options: dict | None = None,
        skip_transpilation: bool = False,
    ):
        """
        Args:
            backend_options: Options passed to QDDSimulator.
            transpile_options: Options passed to transpile.
            skip_transpilation: if True, transpilation is skipped.
        """
        self._circuits = []
        self._parameters = []

        self._backend = QddProvider().get_backend()
        backend_options = {} if backend_options is None else backend_options
        self._backend.set_options(**backend_options)
        self._default_shots = self._backend.options.get("shots")
        self._transpile_options = {} if transpile_options is None else transpile_options
        self._skip_transpilation = skip_transpilation

        self._transpiled_circuits = {}
        self._circuit_ids = {}

    def run(self, pubs: Iterable[SamplerPubLike], *, shots=None, is_exact=False):
        if shots is None:
            shots = self._default_shots
        if is_exact:
            shots = None
        coerced_pubs = [SamplerPub.coerce(pub, shots) for pub in pubs]
        circuits = [pub.circuit for pub in coerced_pubs]
        parameter_values = [pub.parameter_values.as_array() for pub in coerced_pubs]

        from typing import List

        from qiskit.primitives.primitive_job import PrimitiveJob

        circuit_indices: List[int] = []
        for circuit in circuits:
            key = _circuit_key(circuit)
            is_parameterized = circuit.num_parameters > 0

            index = None
            if not is_parameterized:
                index = self._circuit_ids.get(key)

            if index is None:
                index = len(self._circuits)
                self._circuits.append(circuit)
                self._parameters.append(circuit.parameters)
                if not is_parameterized:
                    self._circuit_ids[key] = index

            circuit_indices.append(index)
        job = PrimitiveJob(self._call, circuit_indices, parameter_values, shots)
        job._submit()
        return job

    def _call(
        self,
        circuits: Sequence[int],
        parameter_values,
        shots,
    ) -> SamplerResult:

        is_shots_none = shots is None
        self._transpile(circuits)

        experiment_circuits = []
        parameter_binds = []
        for i, value in zip(circuits, parameter_values):
            if len(value) != len(self._parameters[i]):
                raise QiskitError(
                    f"The number of values ({len(value)}) does not match "
                    f"the number of parameters ({len(self._parameters[i])})."
                )

            experiment_circuits.append(self._transpiled_circuits[i])
            parameter_binds.append(dict(zip(self._parameters[i], value)))

        result = self._backend.run(
            experiment_circuits,
            parameter_binds=parameter_binds,
            shots=shots,
            memory=True,
        ).result()

        results = []
        # Postprocessing
        for i in range(len(experiment_circuits)):
            if is_shots_none:
                probabilities = result.data(i)["probabilities"]
                quasi_dist = QuasiDistribution(probabilities)
                data = DataBin(quasi_dist=quasi_dist)
                metadata = {
                    "shots": None,
                    "simulator_metadata": result.results[i]._metadata,
                }
                results.append(SamplerPubResult(data, metadata))
            else:
                counts = result.get_counts(i)
                quasi_dist = QuasiDistribution(
                    {k.replace(" ", ""): v / shots for k, v in counts.items()},
                    shots=shots,
                )
                memory = result.get_memory(i)
                memory_uint8 = [[np.uint8(int(mem, 2))] for mem in memory]
                bit_array = BitArray(
                    np.array(memory_uint8), result.results[i].meas_level
                )

                data = DataBin(quasi_dist=quasi_dist, meas=bit_array)
                metadata = {
                    "shots": shots,
                    "simulator_metadata": result.results[i]._metadata,
                }
                results.append(SamplerPubResult(data, metadata))

        return results

    def _transpile(self, circuit_indices: Sequence[int]):
        to_handle = [
            i for i in set(circuit_indices) if i not in self._transpiled_circuits
        ]
        if to_handle:
            circuits = (self._circuits[i] for i in to_handle)
            if not self._skip_transpilation:
                circuits = transpile(
                    list(circuits),
                    self._backend,
                    **self._transpile_options,
                )
            for i, circuit in zip(to_handle, circuits):
                self._transpiled_circuits[i] = circuit
