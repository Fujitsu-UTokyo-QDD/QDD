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
from collections.abc import Iterable
import math
from warnings import warn

import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.compiler import transpile
from qiskit.primitives import BaseEstimatorV2
from qiskit.primitives.containers import (
    DataBin,
    EstimatorPubLike,
    PrimitiveResult,
    PubResult,
)
from qiskit.primitives.containers.estimator_pub import EstimatorPub
from qiskit.primitives.primitive_job import PrimitiveJob
from qiskit.primitives.utils import _circuit_key, _observable_key
from qiskit.providers import Options
from qiskit.quantum_info import Pauli, PauliList, SparsePauliOp
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


class Estimator(BaseEstimatorV2):
    """
    QDD implmentation of Estimator.
    """

    def __init__(
        self,
        *,
        backend_options: dict | None = None,
        transpile_options: dict | None = None,
        run_options: dict | None = None,
        default_precision: float = 0.0,
        skip_transpilation: bool = False,
        abelian_grouping: bool = True,
    ):
        """
        Args:
            backend_options: Options passed to QDD backend
            transpile_options: Options passed to transpile.
            run_options: Options passed to run.
            default_precision: Default precision for the estimator. If precision is not
                specified in the run method, this value is used.
            skip_transpilation: If True, transpilation is skipped.
            abelian_grouping: Whether the observable should be grouped into commuting.
                If approximation is True, this parameter is ignored and assumed to be False.
        """
        self._circuits = []
        self._parameters = []
        self._observables = []

        backend_options = {} if backend_options is None else backend_options
        self._backend = QddProvider().get_backend()
        self._backend.set_options(**backend_options)
        self._transpile_options = (
            Options(**transpile_options) if transpile_options else Options()
        )
        self._run_options = Options(**run_options) if run_options else Options()
        self._default_precision = default_precision
        self._skip_transpilation = skip_transpilation
        self._cache: dict[tuple[tuple[int], tuple[int], bool], tuple[dict, dict]] = {}
        self._transpiled_circuits: dict[int, QuantumCircuit] = {}
        self._layouts: dict[int, list[int]] = {}
        self._circuit_ids: dict[tuple, int] = {}
        self._observable_ids: dict[tuple, int] = {}
        self._abelian_grouping = abelian_grouping

    def run(
        self,
        pubs: Iterable[EstimatorPubLike],
        *,
        precision: float | None = None,
    ) -> PrimitiveJob:

        self._precision = (
            precision if precision is not None else self._default_precision
        )

        if self._precision < 0:
            raise ValueError("Precision must be non-negative.")
        elif self._precision == 0:
            self._run_options.update_options(shots=None)
        else:
            shots = int(math.ceil(1 / self._precision**2))
            if shots > self._backend._MAX_SHOTS:
                warn(
                    f"The number of shots ({shots}) is too large. Becaus the precision is too small, "
                    f"It sets the number of shots to the maximum value ({self._backend._MAX_SHOTS})."
                )
                shots = self._backend._MAX_SHOTS
            self._run_options.update_options(
                shots=int(math.ceil(1 / self._precision**2))
            )
        print(self._run_options)

        coerced_pubs = [EstimatorPub.coerce(pub, precision) for pub in pubs]

        job = PrimitiveJob(
            self._run,
            coerced_pubs,
        )
        job._submit()
        return job

    def _run(self, pubs):
        return PrimitiveResult([self._run_pub(pub) for pub in pubs])

    def _run_pub(self, pub):
        circuit_indices = []
        observable_indices = []
        parameter_values = []
        circuit = pub.circuit
        observables = pub.observables.ravel()
        precision = pub.precision
        for observable in observables:
            circuit_index = self._circuit_ids.get(_circuit_key(circuit))
            if circuit_index is not None:
                circuit_indices.append(circuit_index)
            else:
                circuit_indices.append(len(self._circuits))
                self._circuit_ids[_circuit_key(circuit)] = len(self._circuits)
                self._circuits.append(circuit)
                self._parameters.append(circuit.parameters)

            basis = []
            coeffs = []
            for k, v in observable.items():
                basis.append(k)
                coeffs.append(v)
            observable = SparsePauliOp(data=basis, coeffs=coeffs)
            observable_index = self._observable_ids.get(_observable_key(observable))
            if observable_index is not None:
                observable_indices.append(observable_index)
            else:
                observable_indices.append(len(self._observables))
                self._observable_ids[_observable_key(observable)] = len(
                    self._observables
                )
                self._observables.append(observable)

            parameter_values.append(tuple(pub.parameter_values.as_array()))

        expected_values, stds, metadata_list = self._compute(
            circuit_indices, observable_indices, parameter_values, precision
        )
        if len(expected_values) == 1:
            data = DataBin(
                evs=expected_values[0], stds=stds[0], metadata=metadata_list[0]
            )
        else:
            data = DataBin(evs=expected_values, stds=stds, metadata=metadata_list)

        return PubResult(data=data, metadata=None)

    def _compute(self, circuits, observables, parameter_values, precision):
        # Key for cache
        key = (tuple(circuits), tuple(observables), self._precision)

        # Create expectation value experiments.
        if key in self._cache:  # Use a cache
            experiments_dict, obs_maps = self._cache[key]
            exp_map = self._pre_process_params(
                circuits, observables, parameter_values, obs_maps
            )
            experiments, parameter_binds = self._flatten(experiments_dict, exp_map)
            post_processings = self._create_post_processing(
                circuits,
                observables,
                parameter_values,
                obs_maps,
                exp_map,
                self._precision,
                self._run_options.get("shots"),
            )
        else:
            self._transpile_circuits(circuits)
            circ_obs_map = defaultdict(list)
            # Aggregate observables
            for circ_ind, obs_ind in zip(circuits, observables):
                circ_obs_map[circ_ind].append(obs_ind)
            experiments_dict = {}
            obs_maps = (
                {}
            )  # circ_ind => obs_ind => term_ind (Original Pauli) => basis_ind
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
                if (
                    len(bases) == 1 and not bases[0].x.any() and not bases[0].z.any()
                ):  # identity
                    break
                meas_circuits = [
                    self._create_meas_circuit(basis, circ_ind) for basis in bases
                ]
                circuit = (
                    self._circuits[circ_ind]
                    if self._skip_transpilation
                    else self._transpiled_circuits[circ_ind]
                )
                experiments_dict[circ_ind] = self._combine_circs(circuit, meas_circuits)
            self._cache[key] = experiments_dict, obs_maps

            exp_map = self._pre_process_params(
                circuits, observables, parameter_values, obs_maps
            )

            # Flatten
            experiments, parameter_binds = self._flatten(experiments_dict, exp_map)

            # Create PostProcessing
            post_processings = self._create_post_processing(
                circuits,
                observables,
                parameter_values,
                obs_maps,
                exp_map,
                self._precision,
                self._run_options.get("shots"),
            )

        # Run experiments
        if experiments:
            results = (
                self._backend.run(
                    circuits=experiments,
                    parameter_binds=parameter_binds if any(parameter_binds) else None,
                    **self._run_options.__dict__,
                )
                .result()
                .results
            )
        else:
            results = []

        expected_values = []
        stds = []
        metadata_list = []
        for post_processing in post_processings:
            expectation_value, metadata = post_processing.run(results)
            if self._run_options.get("shots") is not None:
                variance = metadata["variance"]
                std = np.sqrt(variance / self._run_options["shots"])
            else:
                std = 0
            expected_values.append(expectation_value)
            stds.append(std)
            metadata_list.append(metadata)
        return expected_values, stds, metadata_list

    def _pre_process_params(self, circuits, observables, parameter_values, obs_maps):
        exp_map = defaultdict(
            dict
        )  # circ_ind => basis_ind => (parameter, parameter_values)
        for circ_ind, obs_ind, param_val in zip(
            circuits, observables, parameter_values
        ):
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
                        [{param: param_val[i] for i, param in enumerate(parameters)}]
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
        opt1q = Optimize1qGatesDecomposition(target=self._backend.target)
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
    def _calculate_result_index(
        circ_ind, obs_ind, term_ind, param_val, obs_maps, exp_map
    ) -> int:
        basis_ind = obs_maps[circ_ind][obs_ind][term_ind]

        result_index = 0
        for _circ_ind, basis_map in exp_map.items():
            for _basis_ind, (_, (_, param_vals)) in enumerate(basis_map.items()):
                if circ_ind == _circ_ind and basis_ind == _basis_ind:
                    result_index += param_vals.index(param_val)
                    return result_index
                result_index += len(param_vals)

    def _create_post_processing(
        self,
        circuits,
        observables,
        parameter_values,
        obs_maps,
        exp_map,
        precision: float,
        shots=None,
        seed=None,
    ) -> list[_PostProcessing]:
        post_processings = []
        for circ_ind, obs_ind, param_val in zip(
            circuits, observables, parameter_values
        ):
            result_indices: list[int | None] = []
            paulis = []
            coeffs = []
            observable = self._observables[obs_ind]
            for term_ind, (pauli, coeff) in enumerate(
                zip(observable.paulis, observable.coeffs)
            ):
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
            post_processings.append(
                _PostProcessing(
                    result_indices,
                    paulis,
                    coeffs,
                    precision,
                    shots,
                    seed,
                )
            )
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
                layout = [bit_map[ci.qubits[0]] for ci in circuit[-num_qubits:]]
                circuit.remove_final_measurements()
                self._transpiled_circuits[i] = circuit
                self._layouts[i] = layout


class _PostProcessing:
    def __init__(
        self,
        result_indices: list[int],
        paulis: list[PauliList],
        coeffs: list[list[float]],
        precision: float,
        shots: int | None = None,
        seed: int = 0,
    ):
        self._result_indices = result_indices
        self._paulis = paulis
        self._coeffs = [np.array(c) for c in coeffs]
        self._precision = precision
        self._shots = shots
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
        if self._shots is not None:
            for c_i, paulis, coeffs in zip(
                self._result_indices, self._paulis, self._coeffs
            ):
                if c_i is None:
                    # Observable is identity
                    expvals, variances = np.array([1], dtype=complex), np.array(
                        [0], dtype=complex
                    )
                else:
                    result = results[c_i]
                    count = result.data.counts
                    basis = result.header.metadata["basis"]
                    indices = np.where(basis.z | basis.x)[0]
                    measured_paulis = PauliList.from_symplectic(
                        paulis.z[:, indices], paulis.x[:, indices], 0
                    )
                    expvals, variances = _pauli_expval_with_variance(
                        count, measured_paulis
                    )
                    simulator_metadata.append(result._metadata)
                combined_expval += np.dot(expvals, coeffs)
                combined_var += np.dot(variances, coeffs**2)

        else:
            for c_i, paulis, coeffs in zip(
                self._result_indices, self._paulis, self._coeffs
            ):
                if c_i is None:
                    # Observable is identity
                    expvals, variances = np.array([1], dtype=complex), np.array(
                        [0], dtype=complex
                    )
                else:
                    result = results[c_i]
                    probabilities = result.data.probabilities
                    quasi_dist = QuasiDistribution(probabilities)
                    basis = result.header.metadata["basis"]
                    indices = np.where(basis.z | basis.x)[0]
                    measured_paulis = PauliList.from_symplectic(
                        paulis.z[:, indices], paulis.x[:, indices], 0
                    )
                    expvals, variances = _pauli_expval_with_variance_from_dist(
                        quasi_dist, measured_paulis
                    )
                    simulator_metadata.append(result._metadata)
                combined_expval += np.dot(expvals, coeffs)
                combined_var += np.dot(variances, coeffs**2)

        combined_expval = np.real_if_close(combined_expval).item()
        metadata = {
            "shots": self._shots,
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


def _pauli_expval_with_variance(
    counts: dict, paulis: PauliList
) -> tuple[np.ndarray, np.ndarray]:
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


def _pauli_expval_with_variance_from_dist(
    quasi_dist: QuasiDistribution, paulis: PauliList
) -> tuple[np.ndarray, np.ndarray]:
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
    packed_vals = np.packbits(
        nonid, axis=1, bitorder="little"
    ).astype(  # pylint:disable=no-member
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
