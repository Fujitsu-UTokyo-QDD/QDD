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

import numpy as np
from qiskit.circuit import ParameterExpression, QuantumCircuit
from qiskit.compiler import transpile
from qiskit.exceptions import QiskitError
from qiskit.primitives import BaseSampler, SamplerResult
from qiskit.primitives.utils import final_measurement_mapping, init_circuit
from qiskit.result import QuasiDistribution

from qdd import QddProvider


class Sampler(BaseSampler):
    """
    QDD implementation of Sampler class.

    :Run Options:

        - **shots** (None or int) --
          The number of shots. If None, it calculates the probabilities exactly.
          Otherwise, it samples from multinomial distributions.

        - **seed** (int) --
          Set a fixed seed for ``seed_simulator``. If shots is None, this option is ignored.

    .. note::
        Precedence of seeding is as follows:

        1. ``seed_simulator`` in runtime (i.e. in :meth:`__call__`)
        2. ``seed`` in runtime (i.e. in :meth:`__call__`)
        3. ``seed_simulator`` of ``backend_options``.
        4. default.
    """

    def __init__(
        self,
        *,
        backend_options: dict | None = None,
        transpile_options: dict | None = None,
        run_options: dict | None = None,
        skip_transpilation: bool = False,
    ):
        """
        Args:
            backend_options: Options passed to QDDSimulator.
            transpile_options: Options passed to transpile.
            run_options: Options passed to run.
            skip_transpilation: if True, transpilation is skipped.
        """
        super().__init__(options=run_options)
        self._circuits = []
        self._parameters = []

        self._backend = QddProvider().get_backend()
        backend_options = {} if backend_options is None else backend_options
        self._backend.set_options(**backend_options)
        self._transpile_options = {} if transpile_options is None else transpile_options
        self._skip_transpilation = skip_transpilation

        self._transpiled_circuits = {}
        self._circuit_ids = {}

    def _call(
        self,
        circuits: Sequence[int],
        parameter_values: Sequence[Sequence[float]],
        **run_options,
    ) -> SamplerResult:
        seed = run_options.pop("seed", None)
        if seed is not None:
            run_options.setdefault("seed_simulator", seed)

        is_shots_none = "shots" in run_options and run_options["shots"] is None
        self._transpile(circuits)

        experiment_manager = _ExperimentManager()
        for i, value in zip(circuits, parameter_values):
            if len(value) != len(self._parameters[i]):
                raise QiskitError(
                    f"The number of values ({len(value)}) does not match "
                    f"the number of parameters ({len(self._parameters[i])})."
                )

            experiment_manager.append(
                key=i,
                parameter_bind=dict(zip(self._parameters[i], value)),
                experiment_circuit=self._transpiled_circuits[i],
            )

        result = self._backend.run(
            experiment_manager.experiment_circuits,
            parameter_binds=experiment_manager.parameter_binds,
            **run_options,
        ).result()

        # Postprocessing
        metadata = []
        quasis = []
        for i in experiment_manager.experiment_indices:
            if is_shots_none:
                probabilities = result.data(i)["probabilities"]
                num_qubits = result.results[i].header.n_qubits
                quasi_dist = QuasiDistribution(probabilities)
                quasis.append(quasi_dist)
                metadata.append(
                    {"shots": None, "simulator_metadata": result.results[i]._metadata}
                )
            else:
                counts = result.get_counts(i)
                shots = sum(counts.values())
                quasis.append(
                    QuasiDistribution(
                        {k.replace(" ", ""): v / shots for k, v in counts.items()},
                        shots=shots,
                    )
                )
                metadata.append(
                    {"shots": shots, "simulator_metadata": result.results[i]._metadata}
                )

        return SamplerResult(quasis, metadata)

    def _run(
        self,
        circuits: Sequence[QuantumCircuit],
        parameter_values: Sequence[Sequence[float]],
        **run_options,
    ):
        from typing import List

        from qiskit.primitives.primitive_job import PrimitiveJob
        from qiskit.primitives.utils import _circuit_key

        circuit_indices: List[int] = []
        for circuit in circuits:
            index = self._circuit_ids.get(_circuit_key(circuit))
            if index is not None:
                circuit_indices.append(index)
            else:
                circuit_indices.append(len(self._circuits))
                self._circuit_ids[_circuit_key(circuit)] = len(self._circuits)
                self._circuits.append(circuit)
                self._parameters.append(circuit.parameters)
        job = PrimitiveJob(self._call, circuit_indices, parameter_values, **run_options)
        job._submit()
        return job

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


class _ExperimentManager:
    def __init__(self):
        self.keys: list[int] = []
        self.experiment_circuits: list[QuantumCircuit] = []
        self.parameter_binds: list[dict[ParameterExpression, float]] = []
        self._input_indices: list[list[int]] = []
        self._num_experiment: int = 0

    def __len__(self):
        return self._num_experiment

    @property
    def experiment_indices(self):
        """indices of experiments"""
        return np.argsort(sum(self._input_indices, [])).tolist()

    def append(
        self,
        key: tuple[int, int],
        parameter_bind: dict[ParameterExpression, float],
        experiment_circuit: QuantumCircuit,
    ):
        """append experiments"""
        if parameter_bind and key in self.keys:
            key_index = self.keys.index(key)
            for k, vs in self.parameter_binds[key_index].items():
                vs.append(parameter_bind[k])
            self._input_indices[key_index].append(self._num_experiment)
        else:
            self.experiment_circuits.append(experiment_circuit)
            self.keys.append(key)
            self.parameter_binds.append({k: v for k, v in parameter_bind.items()})
            self._input_indices.append([self._num_experiment])

        self._num_experiment += 1
