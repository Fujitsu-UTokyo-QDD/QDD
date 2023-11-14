from typing import Any, Dict, List, Optional, Tuple, Union
import traceback
import os
import sys
import uuid
import dataclasses
from collections import Counter
from warnings import warn
import time

from qiskit.providers import BackendV1, JobV1, Options, Provider
from qiskit.providers.models import BackendConfiguration
from qiskit import QuantumCircuit as QiskitCircuit
from qiskit.result import Result
from qiskit.circuit import Barrier, Clbit, Instruction, Measure, ParameterExpression, Qubit, Reset
import qiskit.circuit.library.standard_gates as qiskit_gates
from qiskit.extensions import Initialize

from qdd import __version__
from qdd.qdd_failed_job import QddFailedJob
from qdd.qdd_job import QddJob
from .circuit_property import CircuitProperty

from qdd import pyQDD
from mpi4py import MPI

_qiskit_gates_1q: Dict = {
    qiskit_gates.HGate: "H",
    qiskit_gates.IGate: "I",
    qiskit_gates.SdgGate: "Sdag",
    qiskit_gates.SGate: "S",
    qiskit_gates.SXdgGate: "SXdag",
    qiskit_gates.SXGate: "SX",
    qiskit_gates.TdgGate: "Tdag",
    qiskit_gates.TGate: "T",
    qiskit_gates.XGate: "X",
    qiskit_gates.YGate: "Y",
    qiskit_gates.ZGate: "Z"
}

_qiskit_rotations_1q: Dict = {
    qiskit_gates.RXGate: pyQDD.rxmat,
    qiskit_gates.RYGate: pyQDD.rymat,
    qiskit_gates.RZGate: pyQDD.rzmat,
    qiskit_gates.U1Gate: pyQDD.u1,
    qiskit_gates.U2Gate: pyQDD.u2,
    qiskit_gates.U3Gate: pyQDD.u3,
    qiskit_gates.UGate: pyQDD.u,
    qiskit_gates.PhaseGate: pyQDD.p,
    qiskit_gates.RGate: pyQDD.r,
}

_qiskit_gates_2q: Dict = {
    qiskit_gates.CXGate: pyQDD.CX,
    qiskit_gates.SwapGate: pyQDD.SWAP
}

_qiskit_1q_control: Dict = {
    qiskit_gates.CYGate: "Y",
    qiskit_gates.CZGate: "Z",
    qiskit_gates.CSXGate: "SX",
    qiskit_gates.CCXGate: "X",
    qiskit_gates.C3XGate: "X",
    qiskit_gates.C4XGate: "X",
    qiskit_gates.MCXGate: "X",
}

_supported_qiskit_gates: Dict = {
    **_qiskit_gates_1q,
    **_qiskit_rotations_1q,
    **_qiskit_gates_2q,
    **_qiskit_1q_control,
}

@dataclasses.dataclass
class QddExperiments:
    circs: List[QiskitCircuit]
    circuit_props: List[CircuitProperty]
    options: dict

class QddBackend(BackendV1):
    """A backend used for evaluating circuits with QDD simulator."""

    _save_SV = False

    _DEFAULT_CONFIG: Dict[str, Any] = {
        'backend_name': 'qasm_simulator',
        'backend_version': __version__,
        'n_qubits': 100,  # inclusive
        'basis_gates': sorted([
            'x', 'y', 'z', 'h', 's', 'sdg', 't', 'tdg',
            'id', 'sx', 'sxdg',
            'cx', 'swap',
            'rx', 'ry', 'rz',
            'u1','u2','u3','u','p','r',
            'cy','cz','csx','ccx',"mcx",
            #"cswap", "mcswap", 
            #'cu','cp','cu1','cu2','cu3', # rotation + 1 control
            #"mcu1", "mcu2", "mcu3","mcu","mcp","mcphase", "mcrx", "mcry", "mcrz", "mcr", # rotation + multi controls
            #"rxx","ryy", "rzz", "rzx", # Two-qubits parameterized gates # HOW TO ???
            #"unitary", # HOW TO ???
        ]),
        'gates': [],

        'local': True,
        'simulator': True,
        'conditional': False,
        'open_pulse': False,
        'memory': False,
        'max_shots': int(1e6),  # same as AerSimulator
        'coupling_map': None,
        'description': 'Backend for QDD',
    }

    _DEFAULT_SHOTS = 1024

    # If user specify runtime options that are not included in the default option list, the runtime options are ignored
    # and ignorance warnings will be emitted.
    # This dict is used for suppressing the ignorance warnings.
    _OPTIONS_IGNORED_WITHOUT_WARN = {
        # Entries of "the option key": "whether the option value can be ignored without warning"

        # 'parameter_binds' is a very exceptional option. The option is always not ignored even though it is not in the
        # default option list; so, ignorance warning should not be emitted.
        'parameter_binds': lambda v: True,
        'max_credits': lambda v: True,  # it is obvious to users that max_credits has no meaning in the Qdd simulator
    }

    def __init__(self, provider: Provider, configuration=None):
        if configuration==None:
            configuration = BackendConfiguration.from_dict(QddBackend._DEFAULT_CONFIG)
        super().__init__(
            configuration=configuration,
            provider=provider)

    @classmethod
    def _default_options(cls) -> Options:
        # Note: regarding the 'parameter_binds' option, QddBackend does not include it in the default option list
        # below because AerSimulator also does not.
        # Normally, user-specified runtime options are filtered out in execute(...) if they are not listed below.
        # However, 'parameter_binds' is an exceptional one; it is not excluded regardless of whether to be listed below.
        return Options(
            shots=QddBackend._DEFAULT_SHOTS,
            memory=False,
            seed_simulator=None,
            use_mpi=False,
        )
    
    @staticmethod
    def _create_experiment_header(qiskit_circ: QiskitCircuit) -> dict:
        clbit_labels = []
        creg_sizes = []
        memory_slots = 0
        for creg in qiskit_circ.cregs:
            for i in range(creg.size):
                clbit_labels.append([creg.name, i])
            creg_sizes.append([creg.name, creg.size])
            memory_slots += creg.size

        qubit_labels = []
        qreg_sizes = []
        num_qubits = 0
        for qreg in qiskit_circ.qregs:  # 'qregs' includes ancilla registers
            for i in range(qreg.size):
                qubit_labels.append([qreg.name, i])
            qreg_sizes.append([qreg.name, qreg.size])
            num_qubits += qreg.size

        header = {
            'clbit_labels': clbit_labels,
            'creg_sizes': creg_sizes,
            'global_phase': float(qiskit_circ.global_phase),
            'memory_slots': memory_slots,
            'metadata': qiskit_circ.metadata,
            'n_qubits': num_qubits,
            'name': qiskit_circ.name,
            'qreg_sizes': qreg_sizes,
            'qubit_labels': qubit_labels,
        }

        return header

    def run(self, circuits: Union[QiskitCircuit, List[QiskitCircuit]], **run_options) -> JobV1:
        qiskit_circs: List[QiskitCircuit] = []

        try:
            # convert the given 'circuits' argument into a list of QiskitCircuit
            if isinstance(circuits, QiskitCircuit):
                qiskit_circs = [circuits]
            elif isinstance(circuits, list) and all(isinstance(circ, QiskitCircuit) for circ in circuits):
                qiskit_circs = circuits
            else:
                # the argument type is invalid; abort circuit evaluation.
                if isinstance(circuits, list):
                    arg_circ_type = str(list(map(type, circuits)))
                else:
                    arg_circ_type = str(type(circuits))
                raise RuntimeError(
                    f'{QddBackend.__name__}.run(...) accepts one or more'
                    f' qiskit.{QiskitCircuit.__qualname__} objects only,'
                    f' but the following arguments were specified.{os.linesep}'
                    f'    type={arg_circ_type},{os.linesep}'
                    f'    value={circuits}')

            # evaluate the given circuits
            return self._run_internal(qiskit_circs, **run_options)
        
        except Exception as e:
            # print the stack trace and error messages so that the user can notice the circuit evaluation has failed.
            traceback.print_exc()
            circ_names = [circ.name for circ in qiskit_circs]

            print(f'Failed to evaluate the given circuits. {circ_names}.'
                  f'{os.linesep}Please see error messages.', file=sys.stderr)

            job_id = 'N/A'
            result = Result.from_dict({
                'results': [],
                'backend_name': self.configuration().backend_name,
                'backend_version': self.configuration().backend_version,
                'job_id': job_id,
                'qobj_id': 'N/A',
                'success': False,
                'status': str(e),
            })
            return QddFailedJob(self, job_id, result)

    def _run_internal(self, qiskit_circs: List[QiskitCircuit], **run_options) -> JobV1:
        self._validate_run_options(run_options)
        actual_options = {
            'shots': run_options.get('shots', self.options.shots),
            'memory': run_options.get('memory', self.options.memory),
            'seed_simulator': run_options.get('seed_simulator', self.options.seed_simulator),
            'use_mpi': run_options.get('use_mpi', self.options.use_mpi),
        }

        if ('parameter_binds' in run_options) and (run_options['parameter_binds'] is not None):
            param_bound_qiskit_circs = []
            for qiskit_circ in qiskit_circs:
                for bind in run_options['parameter_binds']:  # the type of run_options['parameter_binds'] is List[Dict]
                    param_bound_qiskit_circs.append(qiskit_circ.bind_parameters(bind))
        else:
            # no parameter bindings are specified
            param_bound_qiskit_circs = qiskit_circs
        
        circ_props = QddBackend._validate_and_get_circuit_properties(param_bound_qiskit_circs, self._save_SV==False)
        experiments = QddExperiments(circs=param_bound_qiskit_circs, circuit_props=circ_props, options=actual_options)

        # run circuits via issuing a job
        job_id = str(uuid.uuid4())
        job = QddJob(backend=self, job_id=job_id, run_exp_fn=self._run_experiment, qobj=experiments)
        job.submit()

        return job

    def _run_experiment(self, experiments, job_id) -> Result:
        """Runs the given experiments"""

        results = [self._evaluate_circuit(circ, circ_prop, experiments.options)
                   for circ, circ_prop
                   in zip(experiments.circs, experiments.circuit_props)]

        result = Result.from_dict({
            'results': results,
            'backend_name': self.configuration().backend_name,
            'backend_version': self.configuration().backend_version,
            'job_id': job_id,
            'qobj_id': 'N/A',
            'success': True,
        })
        return result

    def _create_qubitmap(self, circ: QiskitCircuit):
        qubits = circ.qubits
        mapdata = {}
        id = 0
        for qubit in qubits:
            mapdata[qubit] = id
            id+=1
        self.qubitmap = mapdata

    def get_qID(self, qubit):
        return self.qubitmap[qubit]
    
    def _create_cbitmap(self, circ: QiskitCircuit):
        cbits = circ.clbits
        mapdata = {}
        id = 0
        for cbit in cbits:
            mapdata[cbit] = id
            id+=1
        self.cbitmap = mapdata

    def get_cID(self, cbit):
        return self.cbitmap[cbit]
    
    def _evaluate_circuit(self, circ: QiskitCircuit, circ_prop: CircuitProperty, options: dict):
        use_mpi = options['use_mpi']
        if use_mpi:
            circ = MPI.COMM_WORLD.bcast(circ, root=0)
            circ_prop = MPI.COMM_WORLD.bcast(circ_prop, root=0)
            if pow(2, circ.num_qubits) <= MPI.COMM_WORLD.Get_size():
                print("ERROR: Too many nodes for MPI")
                assert(pow(2, circ.num_qubits) > MPI.COMM_WORLD.Get_size())
            print("MPI enabled")

        start = time.time()
        n_qubit = circ.num_qubits
        n_cbit = circ.num_clbits
        self._create_qubitmap(circ)
        self._create_cbitmap(circ)
        sampled_values = [None] * options['shots']
        print(len(circ.data), " gates")
        if circ_prop.stable_final_state:
            current = pyQDD.makeZeroState(n_qubit) if use_mpi ==False else pyQDD.makeZeroStateMPI(n_qubit)
            for i, qargs, cargs in circ.data:
                qiskit_gate_type = type(i)

                # filter out special cases first
                if qiskit_gate_type == Barrier:
                    continue
                if (qiskit_gate_type == Measure):
                    continue
                assert(len(cargs) == 0)

                if qiskit_gate_type in _supported_qiskit_gates:
                    if qiskit_gate_type in _qiskit_gates_1q:
                        gate = pyQDD.makeGate(n_qubit, _qiskit_gates_1q[qiskit_gate_type], self.get_qID(qargs[0]))
                        current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                    elif qiskit_gate_type in _qiskit_rotations_1q:
                        if qiskit_gate_type == qiskit_gates.U3Gate or qiskit_gate_type == qiskit_gates.UGate:
                            matrix = _qiskit_rotations_1q[qiskit_gate_type](i.params[0],i.params[1],i.params[2])
                        elif qiskit_gate_type == qiskit_gates.U2Gate or qiskit_gate_type == qiskit_gates.RGate:
                            matrix = _qiskit_rotations_1q[qiskit_gate_type](i.params[0],i.params[1])
                        else:
                            matrix = _qiskit_rotations_1q[qiskit_gate_type](i.params[0])
                        gate = pyQDD.makeGate(n_qubit, matrix, self.get_qID(qargs[0]))
                        current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                    elif qiskit_gate_type in _qiskit_gates_2q:
                        gate = _qiskit_gates_2q[qiskit_gate_type](n_qubit, self.get_qID(qargs[1]), self.get_qID(qargs[0]))
                        current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                    elif qiskit_gate_type in _qiskit_1q_control:
                        controls = []
                        for idx in range(len(qargs)-1):
                            controls.append(self.get_qID(qargs[idx]))
                        gate = pyQDD.makeControlGate(n_qubit, _qiskit_1q_control[qiskit_gate_type], self.get_qID(qargs[-1]), controls)
                        current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                    else:
                        raise RuntimeError(f'Unsupported gate or instruction:'
                                       f' type={qiskit_gate_type.__name__}, name={i.name}.'
                                       f' It needs to transpile the circuit before evaluating it.')
                else:
                    if qiskit_gate_type == Measure:
                        continue
                    # We assume the given Qiskit circuit has already been transpiled into a circuit of basis gates only.
                    raise RuntimeError(f'Unsupported gate or instruction:'
                                       f' type={qiskit_gate_type.__name__}, name={i.name}.'
                                       f' It needs to transpile the circuit before evaluating it.')
                current = pyQDD.gc(current);
            
            for i in range(options['shots']):
                _, result_tmp = pyQDD.measureAll(current, False) if use_mpi==False else pyQDD.measureAllMPI(current, False)
                result_final_tmp = ['0'] * n_cbit
                mapping: Dict[Clbit, Qubit] = circ_prop.clbit_final_values
                for cbit in mapping:
                    result_final_tmp[self.get_cID(cbit)] = result_tmp[len(result_tmp)-1-self.get_qID(mapping[cbit])]                
                sampled_values[i] = ''.join(reversed(result_final_tmp))

        else:
            if use_mpi==True:
                print("ERROR: Intermediate measurement is not supported in MPI mode.")
                assert(0)
            for shot in range(options['shots']):
                current = pyQDD.makeZeroState(n_qubit)  if use_mpi ==False else pyQDD.makeZeroStateMPI(n_qubit)
                val_cbit = ['0'] * n_cbit
                for i, qargs, cargs in circ.data:
                    qiskit_gate_type = type(i)

                    # filter out special cases first
                    if qiskit_gate_type == Barrier:
                        continue

                    if qiskit_gate_type in _supported_qiskit_gates:
                        
                        skipGate = False
                        for c_idx in cargs:
                            if val_cbit[self.get_cID(c_idx)] == '0':
                                skipGate = True
                                break
                        if skipGate:
                            continue

                        if qiskit_gate_type in _qiskit_gates_1q:
                            gate = pyQDD.makeGate(n_qubit, _qiskit_gates_1q[qiskit_gate_type], self.get_qID(qargs[0]))
                            current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                        elif qiskit_gate_type in _qiskit_rotations_1q:
                            if qiskit_gate_type == qiskit_gates.U3Gate or qiskit_gate_type == qiskit_gates.UGate:
                                matrix = _qiskit_rotations_1q[qiskit_gate_type](i.params[0],i.params[1],i.params[2])
                            elif qiskit_gate_type == qiskit_gates.U2Gate or qiskit_gate_type == qiskit_gates.RGate:
                                matrix = _qiskit_rotations_1q[qiskit_gate_type](i.params[0],i.params[1])
                            else:
                                matrix = _qiskit_rotations_1q[qiskit_gate_type](i.params[0])
                            gate = pyQDD.makeGate(n_qubit, matrix, self.get_qID(qargs[0]))
                            current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                        elif qiskit_gate_type in _qiskit_gates_2q:
                            gate = _qiskit_gates_2q[qiskit_gate_type](n_qubit, self.get_qID(qargs[1]), self.get_qID(qargs[0]))# target, control
                            current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                        elif qiskit_gate_type in _qiskit_1q_control:
                            controls = []
                            for idx in range(len(qargs)-1):
                                controls.append(self.get_qID(qargs[idx]))
                            gate = pyQDD.makeControlGate(n_qubit, _qiskit_1q_control[qiskit_gate_type], self.get_qID(qargs[-1]), controls)
                            current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiply_MPI(gate, current)
                        else:
                            raise NotImplementedError
                    else:
                        if qiskit_gate_type == Measure:
                            current, val_cbit[self.get_cID(cargs[0])] = pyQDD.measureOneCollapsing(current, self.get_qID(qargs[0]))
                        elif qiskit_gate_type == Reset:
                            current,_meas_result = pyQDD.measureOneCollapsing(current, self.get_qID(qargs[0]))
                            if _meas_result == '1':
                                gate = pyQDD.makeGate(n_qubit, "X", self.get_qID(qargs[0]))
                                current = pyQDD.mv_multiply(gate, current) if use_mpi ==False else pyQDD.mv_multiplyMPI(gate, current)
                        else:
                            # We assume the given Qiskit circuit has already been transpiled into a circuit of basis gates only.
                            raise RuntimeError(f'Unsupported gate or instruction:'
                                       f' type={qiskit_gate_type.__name__}, name={i.name}.'
                                       f' It needs to transpile the circuit before evaluating it.')
                    current = pyQDD.gc(current);
                sampled_values[shot] = ''.join(reversed(val_cbit))

        hex_sampled_counts = Counter(sampled_values)
        result_data: Dict[str, Any] = {'counts': hex_sampled_counts}
        if options['memory']:
            result_data['memory'] = sampled_values
        if self._save_SV:
            result_data["statevector"] = pyQDD.getVector(current) if use_mpi == False else pyQDD.getVectorMPI(current)
        header = QddBackend._create_experiment_header(circ)
        result = {
            'success': True,
            'shots': options['shots'],
            'data': result_data,
            # Note: header information is used by several Qiskit functions; it must be contained in every result object.
            # E.g., header['memory_slots'] is used in Result#get_counts() for formatting sampled counts.
            'header': header,
        }

        print("nQubit", n_qubit, "nGates", len(circ.data), "nNodes", pyQDD.get_nNodes(current))
        return result

    

    @staticmethod
    def _validate_and_get_circuit_properties(qiskit_circs: List[QiskitCircuit], require_measurement=True) -> List[CircuitProperty]:
        circ_props = [QddBackend._validate_and_get_circuit_property(circ, require_measurement) for circ in qiskit_circs]
        return circ_props

    @staticmethod
    def _validate_and_get_circuit_property(circ: QiskitCircuit, require_measurement=True) -> CircuitProperty:
        """Checks whether the given circuit is valid to run (e.g., #qubits <= #max-qubits),
         and, if valid, returns a CircuitProperty instance that will be used for simulation."""

        # #qubit must be lower than #max-qubits
        max_qubits: int = QddBackend._DEFAULT_CONFIG['n_qubits']
        if circ.num_qubits > max_qubits:  # num_qubits includes ancilla registers
            raise RuntimeError(f'Circuit "{circ.name}" has {circ.num_qubits} qubits,'
                               f' but #qubits must be <= {max_qubits}.')

        # Check whether the final state computed by evaluating the circuit is stable over shots,
        # which affects the simulation strategy (see QddBackend._evaluate_circuit(...) for the details).
        # If either of the followings holds, the final state is regarded as unstable.
        # - There are conditional gates.
        # - For the same qubit, there are instructions (except for measurements) after a measurement.
        # - There are probabilistic instructions (e.g., reset, kraus).
        #   Note: Reset is decomposed to two operations: Measure(q, c) and circuit.x(q).c_if(c, 1);
        #         so, Reset is regarded as a stochastic operation.
        # - There are 'initialize' instructions that are not located at the first step in the circuit
        #   and do not operate on all qubits (i.e., non-full-width initialization at or after the second step).
        #   Note: 'Initialize' performs reset instructions at first, thereby regarded as stochastic.
        #         Exceptions are full-width initializations and an initialization at the first step,
        #         which is deterministic.
        #
        # Note: consecutive measurements for the same qubit do not result in an unstable final state
        # because they are equal to a single measurement for the qubit.
        # Note: a measure gate with broadcast parameters, e.g., Measure(0, [0, 1]), is decomposed into
        # two measure gates: Measure(0, 0) and Measure(0, 1).
        # Note: all the quantum channel instructions (e.g., Choi, SuperOp, Kraus, ...) have the name of 'kraus'.
        stable_final_state: bool = True
        measured_qubits = set()
        clbit_final_values: Dict[Clbit, Qubit] = {}  # for each clbit, this holds the qubit last assigned to the clbit.
        for i, (inst, qargs, cargs) in enumerate(circ.data):
            if type(inst) == Measure:
                # For a Measure instruction, both qargs and cargs always have a size of 1.
                # (Multi-target measurements (Measure([...], [...]) is decomposed to single-target measurement gates.)
                clbit_final_values[cargs[0]] = qargs[0]
                measured_qubits |= set(qargs)
            elif inst.condition is not None \
                    or any(qubit in measured_qubits for qubit in qargs) \
                    or type(inst) == Reset \
                    or inst.name == 'kraus' \
                    or inst.name == 'superop':
                # Note: the condition value of 'superop' might have no meaning because SuperOp instances have
                # the name of 'kraus'; however, the implementation of the Qiskit aer simulator explicitly tests
                # the name against 'superop' to check whether the circuit has probabilistic instructions.
                # https://github.com/Qiskit/qiskit-aer/blob/0.9.1/src/controllers/aer_controller.hpp#L1621
                # So, we do the same check here just in case.
                stable_final_state = False
            elif type(inst) == Initialize:
                if i != 0 and len(qargs) < circ.num_qubits:
                    stable_final_state = False

            # if all the properties this loop try to detect have already been identified,
            # there's no need to check further
            if not stable_final_state and len(measured_qubits) > 0:
                break

        # Notice: here, 'measured_qubits' may not contain all measured qubits due to 'break' above

        # The circuit must have measurement gates
        if not measured_qubits:
            if require_measurement:
                raise RuntimeError(f'Circuit "{circ.name}" has no measurement gates.'
                               f' Every circuit must have measurement gates when using {QddBackend.__name__}.')

        return CircuitProperty(stable_final_state=stable_final_state, clbit_final_values=clbit_final_values)
    
    def _validate_run_options(self, run_options):
        """Checks whether the given options are valid to run with (e.g., shots < max_shots)."""

        # Check if invalid options (i.e., causing errors) are specified
        if 'shots' in run_options:
            shots = run_options['shots']
            max_shots = QddBackend._DEFAULT_CONFIG['max_shots']
            if not (0 < shots <= max_shots):
                raise RuntimeError(f'Shots={shots} is specified, but #shots must be positive and <= {max_shots}'
                                   f' when using {QddBackend.__name__}.')

        # Warn if unsupported options are specified
        # Note: we cannot raise an error here because some Qiskit functions (e.g., execute(...)) adds some options to
        # the user-specified options (i.e., user-unspecified options sometimes comes to here).
        for run_opt in run_options:
            # if the runtime option is not contained in the default option list, output a warning.
            if not hasattr(self.options, run_opt):
                # some options can be ignored without warnings (e.g., an option value of no effect)
                can_ignore = run_opt in QddBackend._OPTIONS_IGNORED_WITHOUT_WARN \
                    and QddBackend._OPTIONS_IGNORED_WITHOUT_WARN[run_opt](run_options[run_opt])
                if not can_ignore:
                    warn(f'Option {run_opt}={run_options[run_opt]} is not used by {QddBackend.__name__}.',
                         UserWarning, stacklevel=2)