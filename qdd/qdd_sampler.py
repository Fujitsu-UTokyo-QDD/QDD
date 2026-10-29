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
    """QDD implementation of Qiskit's Sampler V2 primitive.

    The sampler evaluates circuits with the QDD qasm backend and returns
    quasi-distributions, and shot memory when finite shots are requested.
    Shot memory is exposed as a BitArray under each classical register's name.
    If a register is named ``quasi_dist``, the joint quasi-distribution is
    stored in the pub result's metadata instead of its data.
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
        """Submit sampler primitive inputs for execution.

        Args:
            pubs: Iterable of sampler primitive inputs accepted by Qiskit.
            shots: Default number of shots for pubs without their own shot count.
                If omitted, the backend default is used.
            is_exact: If true, compute exact probabilities by running with
                ``shots=None``, overriding both pub and run shot counts.

        Returns:
            A Qiskit ``PrimitiveJob`` whose result contains sampler results.
        """
        if shots is None:
            shots = self._default_shots
        if is_exact:
            shots = None
        coerced_pubs = [SamplerPub.coerce(pub, shots) for pub in pubs]
        circuits = [pub.circuit for pub in coerced_pubs]
        parameter_values = [pub.parameter_values.as_array() for pub in coerced_pubs]
        pub_shots = [None if is_exact else pub.shots for pub in coerced_pubs]

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
        job = PrimitiveJob(self._run, circuit_indices, parameter_values, pub_shots)
        job._submit()
        return job

    def _run(self, circuits, parameter_values, pub_shots):
        # The backend accepts one shot count per run; group matching pubs.
        shot_groups = {}
        for index, shots in enumerate(pub_shots):
            shot_groups.setdefault(shots, []).append(index)

        results = [None] * len(circuits)
        for shots, indices in shot_groups.items():
            pub_results = self._call(
                [circuits[index] for index in indices],
                [parameter_values[index] for index in indices],
                shots,
            )
            # Preserve input order even when equal-shot pubs are interleaved.
            for index, result in zip(indices, pub_results):
                results[index] = result
        return results

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
        if not result.success:
            raise QiskitError(result.status)

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
                circuit = experiment_circuits[i]
                # Raw memory has no register separators and follows circuit bit order.
                memory = result.data(i)["memory"]
                bit_array = BitArray.from_samples(
                    memory, num_bits=circuit.num_clbits
                )
                register_data = {
                    register.name: bit_array.slice_bits(
                        [circuit.find_bit(bit).index for bit in register]
                    )
                    for register in circuit.cregs
                }

                metadata = {
                    "shots": shots,
                    "simulator_metadata": result.results[i]._metadata,
                }
                if "quasi_dist" in register_data:
                    metadata["quasi_dist"] = quasi_dist
                else:
                    register_data["quasi_dist"] = quasi_dist
                data = DataBin(**register_data)
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
