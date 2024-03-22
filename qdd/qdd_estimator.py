# This code is part of Qiskit.
#
# (C) Copyright IBM 2022, 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# Code adapted from:
# https://github.com/Qiskit/qiskit-aer/blob/main/qiskit_aer/primitives/estimator.py

"""
Estimator class.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from warnings import warn

import numpy as np
from qiskit.circuit import ParameterExpression, QuantumCircuit
from qiskit.compiler import transpile
from qiskit.primitives import BaseEstimator, EstimatorResult
from qiskit.primitives.primitive_job import PrimitiveJob
from qiskit.primitives.utils import _circuit_key, _observable_key, init_observable
from qiskit.providers import Options,convert_to_target
from qiskit.quantum_info import Pauli, PauliList
from qiskit.quantum_info.operators.base_operator import BaseOperator
from qiskit.result import QuasiDistribution
from qiskit.result.models import ExperimentResult
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import (
    ApplyLayout,
    EnlargeWithAncilla,
    FullAncillaAllocation,
    Optimize1qGatesDecomposition,
    SetLayout,
)

from qdd import QddProvider

class Estimator(BaseEstimator):
    """
    QDD implmentation of Estimator.

    :Run Options:
        - **shots** (None or int) --
          The number of shots. If None and approximation is True, it calculates the exact
          expectation values. Otherwise, it calculates expectation values with sampling.

        - **seed** (int) --
          Set a fixed seed for the sampling.

    .. note::
        Precedence of seeding for ``seed_simulator`` is as follows:

        1. ``seed_simulator`` in runtime (i.e. in :meth:`run`)
        2. ``seed`` in runtime (i.e. in :meth:`run`)
        3. ``seed_simulator`` of ``backend_options``.
        4. default.

        ``seed`` is also used for sampling from a normal distribution when approximation is True.

        When combined with the approximation option, we get the expectation values as follows:

        * shots is None and approximation=False: Return an expectation value with sampling-noise w/
          warning.
        * shots is int and approximation=False: Return an expectation value with sampling-noise.
        * shots is None and approximation=True: Return an exact expectation value.
        * shots is int and approximation=True: Return expectation value with sampling-noise using a
          normal distribution approximation.
    """

    def __init__(
        self,
        *,
        backend_options: dict | None = None,
        transpile_options: dict | None = None,
        run_options: dict | None = None,
        approximation: bool = False,
        skip_transpilation: bool = False,
        abelian_grouping: bool = True,
    ):
        """
        Args:
            backend_options: Options passed to QDD backend
            transpile_options: Options passed to transpile.
            run_options: Options passed to run.
            approximation: If True, it calculates expectation values with normal distribution
                approximation. Note that this appproximation ignores readout errors.
            skip_transpilation: If True, transpilation is skipped.
            abelian_grouping: Whether the observable should be grouped into commuting.
                If approximation is True, this parameter is ignored and assumed to be False.
        """
        super().__init__(options=run_options)
        # These three private attributes used to be created by super, but were deprecated in Qiskit
        # 0.46. See https://github.com/Qiskit/qiskit/pull/11051
        self._circuits = []
        self._parameters = []
        self._observables = []

        backend_options = {} if backend_options is None else backend_options
        self._backend = QddProvider().get_backend()
        self._backend.set_options(**backend_options)
        self._transpile_options = Options()
        if transpile_options is not None:
            self._transpile_options.update_options(**transpile_options)
        self._approximation = approximation
        if self._approximation and "shots" in run_options:
            self._approximation_shots = run_options["shots"]
        else:
            self._approximation_shots = None
        self._skip_transpilation = skip_transpilation
        self._cache: dict[tuple[tuple[int], tuple[int], bool], tuple[dict, dict]] = {}
        self._transpiled_circuits: dict[int, QuantumCircuit] = {}
        self._layouts: dict[int, list[int]] = {}
        self._circuit_ids: dict[tuple, int] = {}
        self._observable_ids: dict[tuple, int] = {}
        self._abelian_grouping = abelian_grouping

    def _call(
        self,
        circuits: Sequence[int],
        observables: Sequence[int],
        parameter_values: Sequence[Sequence[float]],
        **run_options,
    ) -> EstimatorResult:
        seed = run_options.pop("seed", None)
        if seed is not None:
            run_options.setdefault("seed_simulator", seed)

        return self._compute(circuits, observables, parameter_values, run_options, seed)

    def _run(
        self,
        circuits: Sequence[QuantumCircuit],
        observables: Sequence[BaseOperator],
        parameter_values: Sequence[Sequence[float]],
        **run_options,
    ) -> PrimitiveJob:
        circuit_indices: list = []
        for circuit in circuits:
            index = self._circuit_ids.get(_circuit_key(circuit))
            if index is not None:
                circuit_indices.append(index)
            else:
                circuit_indices.append(len(self._circuits))
                self._circuit_ids[_circuit_key(circuit)] = len(self._circuits)
                self._circuits.append(circuit)
                self._parameters.append(circuit.parameters)
        observable_indices: list = []
        for observable in observables:
            observable = init_observable(observable)
            index = self._observable_ids.get(_observable_key(observable))
            if index is not None:
                observable_indices.append(index)
            else:
                observable_indices.append(len(self._observables))
                self._observable_ids[_observable_key(observable)] = len(self._observables)
                self._observables.append(observable)
        job = PrimitiveJob(
            self._call,
            circuit_indices,
            observable_indices,
            parameter_values,
            **run_options,
        )
        job._submit()
        return job

    def _compute(self, circuits, observables, parameter_values, run_options, seed):

        # Key for cache
        key = (tuple(circuits), tuple(observables), self._approximation)

        # Create expectation value experiments.
        if key in self._cache:  # Use a cache
            experiments_dict, obs_maps = self._cache[key]
            exp_map = self._pre_process_params(circuits, observables, parameter_values, obs_maps)
            experiments, parameter_binds = self._flatten(experiments_dict, exp_map)
            post_processings = self._create_post_processing(
                circuits, observables, parameter_values, obs_maps, exp_map, self._approximation, self._approximation_shots, seed
            )
        else:
            self._transpile_circuits(circuits)
            circ_obs_map = defaultdict(list)
            # Aggregate observables
            for circ_ind, obs_ind in zip(circuits, observables):
                circ_obs_map[circ_ind].append(obs_ind)
            experiments_dict = {}
            obs_maps = {}  # circ_ind => obs_ind => term_ind (Original Pauli) => basis_ind
            # Group and create measurement circuit
            for circ_ind, obs_indices in circ_obs_map.items():
                pauli_list = sum(
                    [self._observables[obs_ind].paulis for obs_ind in obs_indices]
                ).unique()
                if self._abelian_grouping:
                    pauli_lists = pauli_list.group_commuting(qubit_wise=True)
                else:
                    pauli_lists = [PauliList(pauli) for pauli in pauli_list]
                obs_map = defaultdict(list)
                for obs_ind in obs_indices:
                    for pauli in self._observables[obs_ind].paulis:
                        for basis_ind, pauli_list in enumerate(pauli_lists):
                            if pauli in pauli_list:
                                obs_map[obs_ind].append(basis_ind)
                                break
                obs_maps[circ_ind] = obs_map
                bases = [_paulis2basis(pauli_list) for pauli_list in pauli_lists]
                if len(bases) == 1 and not bases[0].x.any() and not bases[0].z.any():  # identity
                    break
                meas_circuits = [self._create_meas_circuit(basis, circ_ind) for basis in bases]
                circuit = (
                    self._circuits[circ_ind]
                    if self._skip_transpilation
                    else self._transpiled_circuits[circ_ind]
                )
                experiments_dict[circ_ind] = self._combine_circs(circuit, meas_circuits)
            self._cache[key] = experiments_dict, obs_maps

            exp_map = self._pre_process_params(circuits, observables, parameter_values, obs_maps)

            # Flatten
            experiments, parameter_binds = self._flatten(experiments_dict, exp_map)

            # Create PostProcessing
            post_processings = self._create_post_processing(
                circuits, observables, parameter_values, obs_maps, exp_map, self._approximation, self._approximation_shots, seed
            )

        if self._approximation:
            run_options["shots"] = None
        if (not self._approximation) and ("shots" not in run_options or run_options["shots"] is None):
            warn(
                "The number of shots is not specified. "
                f"It sets the number of shots to default value ({self._backend._DEFAULT_SHOTS})."
            )
            run_options["shots"] = self._backend._DEFAULT_SHOTS
        # Run experiments
        if experiments:
            results = (
                self._backend.run(
                    circuits=experiments,
                    parameter_binds=parameter_binds if any(parameter_binds) else None,
                    **run_options,
                )
                .result()
                .results
            )
        else:
            results = []

       # Post processing (calculate expectation values)
        expectation_values, metadata = zip(
            *(post_processing.run(results) for post_processing in post_processings)
        )
        return EstimatorResult(np.real_if_close(expectation_values), list(metadata))

    def _pre_process_params(self, circuits, observables, parameter_values, obs_maps):
        exp_map = defaultdict(dict)  # circ_ind => basis_ind => (parameter, parameter_values)
        for circ_ind, obs_ind, param_val in zip(circuits, observables, parameter_values):
            self._validate_parameter_length(param_val, circ_ind)
            parameter = self._parameters[circ_ind]
            for basis_ind in obs_maps[circ_ind][obs_ind]:
                if (
                    circ_ind in exp_map
                    and basis_ind in exp_map[circ_ind]
                    and len(self._parameters[circ_ind]) > 0
                ):
                    param_vals = exp_map[circ_ind][basis_ind][1]
                    if param_val not in param_vals:
                        param_vals.append(param_val)
                else:
                    exp_map[circ_ind][basis_ind] = (parameter, [param_val])

        return exp_map

    @staticmethod
    def _flatten(experiments_dict, exp_map):
        experiments_unsorted = {}
        experiments_list = []
        parameter_binds = []
        for circ_ind in experiments_dict:
            indices = []
            for i, (parameters, param_vals) in exp_map[circ_ind].items():
                for param_val in param_vals:
                    experiments_unsorted[i] = experiments_dict[circ_ind][i]
                    indices.append(i)
                    parameter_binds.extend(
                        [
                            {
                                param: param_val[i]
                                for i, param in enumerate(parameters)
                            }
                        ]
                )
            for i in sorted(indices):
                experiments_list.append(experiments_unsorted[i])
                
        return experiments_list, parameter_binds

    def _create_meas_circuit(self, basis: Pauli, circuit_index: int):
        qargs = np.arange(basis.num_qubits)[basis.z | basis.x]
        meas_circuit = QuantumCircuit(basis.num_qubits, len(qargs))
        for clbit, qarg in enumerate(qargs):
            if basis.x[qarg]:
                if basis.z[qarg]:
                    meas_circuit.sdg(qarg)
                meas_circuit.h(qarg)
            meas_circuit.measure(qarg, clbit)
        meas_circuit.metadata = {"basis": basis}

        if self._skip_transpilation:
            return meas_circuit

        layout = self._layouts[circuit_index]
        passmanager = PassManager([SetLayout(layout)])
        opt1q = Optimize1qGatesDecomposition(target=convert_to_target(self._backend.configuration()))
        passmanager.append(opt1q)
        if isinstance(self._backend.coupling_map, CouplingMap):
            coupling_map = self._backend.coupling_map
            passmanager.append(FullAncillaAllocation(coupling_map))
            passmanager.append(EnlargeWithAncilla())
        passmanager.append(ApplyLayout())
        return passmanager.run(meas_circuit)

    @staticmethod
    def _combine_circs(circuit: QuantumCircuit, meas_circuits: list[QuantumCircuit]):
        circs = []
        for meas_circuit in meas_circuits:
            new_circ = circuit.copy()
            for creg in meas_circuit.cregs:
                new_circ.add_register(creg)
            new_circ.compose(meas_circuit, inplace=True)
            _update_metadata(new_circ, meas_circuit.metadata)
            circs.append(new_circ)
        return circs

    @staticmethod
    def _calculate_result_index(circ_ind, obs_ind, term_ind, param_val, obs_maps, exp_map) -> int:
        basis_ind = obs_maps[circ_ind][obs_ind][term_ind]

        result_index = 0
        for _circ_ind, basis_map in exp_map.items():
            for _basis_ind, (_, (_, param_vals)) in enumerate(basis_map.items()):
                if circ_ind == _circ_ind and basis_ind == _basis_ind:
                    result_index += param_vals.index(param_val)
                    return result_index
                result_index += len(param_vals)

    def _create_post_processing(
        self, circuits, observables, parameter_values, obs_maps, exp_map, approximation=False, approximation_shots=None, seed=None
    ) -> list[_PostProcessing]:
        post_processings = []
        for circ_ind, obs_ind, param_val in zip(circuits, observables, parameter_values):
            result_indices: list[int | None] = []
            paulis = []
            coeffs = []
            observable = self._observables[obs_ind]
            for term_ind, (pauli, coeff) in enumerate(zip(observable.paulis, observable.coeffs)):
                # Identity
                if not pauli.x.any() and not pauli.z.any():
                    result_indices.append(None)
                    paulis.append(PauliList(pauli))
                    coeffs.append([coeff])
                    continue

                result_index = self._calculate_result_index(
                    circ_ind, obs_ind, term_ind, param_val, obs_maps, exp_map
                )
                if result_index in result_indices:
                    i = result_indices.index(result_index)
                    paulis[i] += pauli
                    coeffs[i].append(coeff)
                else:
                    result_indices.append(result_index)
                    paulis.append(PauliList(pauli))
                    coeffs.append([coeff])
            post_processings.append(_PostProcessing(result_indices, paulis, coeffs, approximation, approximation_shots, seed))
        return post_processings

    def _validate_parameter_length(self, parameter, circuit_index):
        if len(parameter) != len(self._parameters[circuit_index]):
            raise ValueError(
                f"The number of values ({len(parameter)}) does not match "
                f"the number of parameters ({len(self._parameters[circuit_index])})."
            )

    def _transpile(self, circuits):
        if self._skip_transpilation:
            return circuits
        return transpile(circuits, self._backend, **self._transpile_options.__dict__)

    def _transpile_circuits(self, circuits):
        if self._skip_transpilation:
            for i in set(circuits):
                num_qubits = self._circuits[i].num_qubits
                self._layouts[i] = list(range(num_qubits))
            return
        for i in set(circuits):
            if i not in self._transpiled_circuits:
                circuit = self._circuits[i].copy()
                circuit.measure_all()
                num_qubits = circuit.num_qubits
                circuit = self._transpile(circuit)
                bit_map = {bit: index for index, bit in enumerate(circuit.qubits)}
                layout = [bit_map[qr[0]] for _, qr, _ in circuit[-num_qubits:]]
                circuit.remove_final_measurements()
                self._transpiled_circuits[i] = circuit
                self._layouts[i] = layout

class _PostProcessing:
    def __init__(
        self,
        result_indices: list[int],
        paulis: list[PauliList],
        coeffs: list[list[float]],
        approximation: bool = False,
        approximation_shots: int | None = None,
        seed: int = 0,
    ):
        self._result_indices = result_indices
        self._paulis = paulis
        self._coeffs = [np.array(c) for c in coeffs]
        self._approximation = approximation
        self._approximation_shots = approximation_shots
        self._seed = seed

    def run(self, results: list[ExperimentResult]) -> tuple[float, dict]:
        """Coumpute expectation value.

        Args:
            results: list of ExperimentResult.

        Returns:
            tuple of an expectation value and metadata.
        """
        combined_expval = 0.0
        combined_var = 0.0
        simulator_metadata = []
        if self._approximation:
            for c_i, paulis, coeffs in zip(self._result_indices, self._paulis, self._coeffs):
                if c_i is None:
                    # Observable is identity
                    expvals, variances = np.array([1], dtype=complex), np.array([0], dtype=complex)
                    shots = 0
                else:
                    result = results[c_i]
                    probabilities = result.data.probabilities
                    quasi_dist = QuasiDistribution(probabilities)
                    basis = result.header.metadata["basis"]
                    indices = np.where(basis.z | basis.x)[0]
                    measured_paulis = PauliList.from_symplectic(
                        paulis.z[:, indices], paulis.x[:, indices], 0
                    )
                    expvals, variances = _pauli_expval_with_variance_from_dist(quasi_dist, measured_paulis)
                    shots = results[c_i].shots
                combined_expval += np.dot(expvals, coeffs)
                combined_var += np.dot(variances, coeffs**2)
            if self._approximation_shots:
                rng = np.random.default_rng(self._seed)
                standard_error = np.sqrt(combined_var / self._approximation_shots)
                combined_expval = rng.normal(combined_expval, standard_error)
                metadata = {
                    "shots": self._approximation_shots,
                    "variance": combined_var,
                    "simulator_metadata": simulator_metadata,
                }
            else:
                metadata = {
                    "shots": shots,
                    "simulator_metadata": simulator_metadata,
                }
        else:
            for c_i, paulis, coeffs in zip(self._result_indices, self._paulis, self._coeffs):
                if c_i is None:
                    # Observable is identity
                    expvals, variances = np.array([1], dtype=complex), np.array([0], dtype=complex)
                    shots = None
                else:
                    result = results[c_i]
                    count = result.data.counts
                    shots = sum(count.values())
                    basis = result.header.metadata["basis"]
                    indices = np.where(basis.z | basis.x)[0]
                    measured_paulis = PauliList.from_symplectic(
                        paulis.z[:, indices], paulis.x[:, indices], 0
                    )
                    expvals, variances = _pauli_expval_with_variance(count, measured_paulis)
                    simulator_metadata.append(result._metadata)
                combined_expval += np.dot(expvals, coeffs)
                combined_var += np.dot(variances, coeffs**2)
            metadata = {
                "shots": shots,
                "variance": np.real_if_close(combined_var).item(),
                "simulator_metadata": simulator_metadata,
            }
        return combined_expval, metadata


def _update_metadata(circuit: QuantumCircuit, metadata: dict) -> QuantumCircuit:
    if circuit.metadata:
        circuit.metadata.update(metadata)
    else:
        circuit.metadata = metadata
    return circuit


def _pauli_expval_with_variance(counts: dict, paulis: PauliList) -> tuple[np.ndarray, np.ndarray]:
    # Diag indices
    size = len(paulis)
    diag_inds = _paulis2inds(paulis)

    expvals = np.zeros(size, dtype=float)
    denom = 0  # Total shots for counts dict
    for bin_outcome, freq in counts.items():
        outcome = int(bin_outcome, 2)
        denom += freq
        for k in range(size):
            coeff = (-1) ** _parity(diag_inds[k] & outcome)
            expvals[k] += freq * coeff

    # Divide by total shots
    expvals /= denom

    variances = 1 - expvals**2
    return expvals, variances

def _pauli_expval_with_variance_from_dist(quasi_dist: QuasiDistribution, paulis: PauliList) -> tuple[np.ndarray, np.ndarray]:
    # Diag indices
    size = len(paulis)
    diag_inds = _paulis2inds(paulis)

    expvals = np.zeros(size, dtype=float)
    for outcome, freq in quasi_dist.items():
        for k in range(size):
            coeff = (-1) ** _parity(diag_inds[k] & outcome)
            expvals[k] += freq * coeff


    variances = 1 - expvals**2
    return expvals, variances


def _paulis2inds(paulis: PauliList) -> list[int]:
    nonid = paulis.z | paulis.x
    packed_vals = np.packbits(nonid, axis=1, bitorder="little").astype(  # pylint:disable=no-member
        object
    )
    power_uint8 = 1 << (8 * np.arange(packed_vals.shape[1], dtype=object))
    inds = packed_vals @ power_uint8
    return inds.tolist()


def _parity(integer: int) -> int:
    """Return the parity of an integer"""
    return bin(integer).count("1") % 2


def _paulis2basis(paulis: PauliList) -> Pauli:
    return Pauli(
        (
            np.logical_or.reduce(paulis.z),  # pylint:disable=no-member
            np.logical_or.reduce(paulis.x),  # pylint:disable=no-member
        )
    )


class _ExperimentManager:
    def __init__(self):
        self.keys: list[tuple[int, int]] = []
        self.experiment_circuits: list[QuantumCircuit] = []
        self.parameter_binds: list[dict[ParameterExpression, list[float]]] = []
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
        if key in self.keys and parameter_bind:
            key_index = self.keys.index(key)
            for k, vs in self.parameter_binds[key_index].items():
                vs.append(parameter_bind[k])
            self._input_indices[key_index].append(self._num_experiment)
        else:
            self.experiment_circuits.append(experiment_circuit)
            self.keys.append(key)
            self.parameter_binds.append({k: [v] for k, v in parameter_bind.items()})
            self._input_indices.append([self._num_experiment])

        self._num_experiment += 1

    def get_observable_key(self, index):
        """return key of observables"""
        for i, inputs in enumerate(self._input_indices):
            if index in inputs:
                return self.keys[i][1]
        raise RuntimeError("Unexpected behavior.")
