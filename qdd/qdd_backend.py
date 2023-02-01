from typing import Any, Dict, List, Optional, Tuple, Union
import traceback
import os
import sys
import uuid
import dataclasses

from qiskit.providers import BackendV1, JobV1, Options, Provider
from qiskit.providers.models import BackendConfiguration
from qiskit import QuantumCircuit as QiskitCircuit
from qiskit.result import Result
from qiskit.circuit import Barrier, Clbit, Instruction, Measure, ParameterExpression, Qubit, Reset
import qiskit.circuit.library.standard_gates as qiskit_gates

from qdd import __version__
from qdd.qdd_failed_job import QddFailedJob
from qdd.qdd_job import QddJob

import pyQDD

_qiskit_gates_1q: Dict = {
    qiskit_gates.HGate: pyQDD.H,
    qiskit_gates.IGate: pyQDD.I,
    qiskit_gates.SdgGate: pyQDD.Sdag,
    qiskit_gates.SGate: pyQDD.S,
    qiskit_gates.SXdgGate: pyQDD.SXdag,
    qiskit_gates.SXGate: pyQDD.SX,
    qiskit_gates.TdgGate: pyQDD.Tdag,
    qiskit_gates.TGate: pyQDD.T,
    qiskit_gates.XGate: pyQDD.X,
    qiskit_gates.YGate: pyQDD.Y,
    qiskit_gates.ZGate: pyQDD.Z
}

_qiskit_rotations_1q: Dict = {
    qiskit_gates.RXGate: pyQDD.RX,
    qiskit_gates.RYGate: pyQDD.RY,
    qiskit_gates.RZGate: pyQDD.RZ,
}

_qiskit_gates_2q: Dict = {
    qiskit_gates.CXGate: pyQDD.CX,
    qiskit_gates.CYGate: pyQDD.CZ,
    qiskit_gates.CZGate: pyQDD.CZ,
    qiskit_gates.SwapGate: pyQDD.SWAP
}

_supported_qiskit_gates: Dict = {
    **_qiskit_gates_1q,
    **_qiskit_rotations_1q,
    **_qiskit_gates_2q,
    Measure: None
}

@dataclasses.dataclass
class QddExperiments:
    circs: List[QiskitCircuit]
    options: dict

class QddBackend(BackendV1):
    """A backend used for evaluating circuits with QDD simulator."""

    _DEFAULT_CONFIG: Dict[str, Any] = {
        'backend_name': 'qdd_backend',
        'backend_version': __version__,
        'n_qubits': 100,  # inclusive
        'basis_gates': sorted([
            'x', 'y', 'z', 'h', 's', 'sdg', 't', 'tdg',
            'id', 'sx', 'sxdg',
            'cx', 'cy', 'cz', 'swap',
            'rx', 'ry', 'rz',
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
        'max_credits': lambda v: True,  # it is obvious to users that max_credits has no meaning in the Qulacs simulator
    }

    def __init__(self, provider: Provider):
        super().__init__(
            configuration=BackendConfiguration.from_dict(QddBackend._DEFAULT_CONFIG),
            provider=provider)

    @classmethod
    def _default_options(cls) -> Options:
        # Note: regarding the 'parameter_binds' option, QulacsBackend does not include it in the default option list
        # below because AerSimulator also does not.
        # Normally, user-specified runtime options are filtered out in execute(...) if they are not listed below.
        # However, 'parameter_binds' is an exceptional one; it is not excluded regardless of whether to be listed below.
        return Options(
            shots=QddBackend._DEFAULT_SHOTS,
            memory=False,
            seed_simulator=None,
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
        actual_options = {
            'shots': run_options.get('shots', self.options.shots),
            #'memory': run_options.get('memory', self.options.memory),
            'seed_simulator': run_options.get('seed_simulator', self.options.seed_simulator),
        }
        
        experiments = QddExperiments(circs=qiskit_circs, options=actual_options)

        # run circuits via issuing a job
        job_id = str(uuid.uuid4())
        job = QddJob(backend=self, job_id=job_id, run_exp_fn=self._run_experiment, qobj=experiments)
        job.submit()

        return job

    def _run_experiment(self, experiments, job_id) -> Result:
        """Runs the given experiments"""

        results = [self._evaluate_circuit(circuit, experiments.options) for circuit in experiments.circs]
        result = Result.from_dict({
            'results': results,
            'backend_name': self.configuration().backend_name,
            'backend_version': self.configuration().backend_version,
            'job_id': job_id,
            'qobj_id': 'N/A',
            'success': True,
        })
        return result
    
    def add_gate(self, i, qargs, cargs, ddcirc):
        gatetype = type(i)
        fn = _supported_qiskit_gates[type(i)]
        if gatetype in _qiskit_gates_1q:
            fn(self.get_qID(qargs[0]), ddcirc)
        elif gatetype in _qiskit_gates_2q:
            fn(self.get_qID(qargs[0]), self.get_qID(qargs[1]), ddcirc)
        elif gatetype in _qiskit_rotations_1q:
            fn(i.params[0], self.get_qID(qargs[0]), ddcirc)
        else:
            print("Measure ignored")

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

    def _evaluate_circuit(self, circ: QiskitCircuit, options: dict):
        n_qubit = circ.num_qubits
        n_cbit = circ.num_clbits
        self._create_qubitmap(circ)
        if n_cbit>0:
            raise RuntimeError('Classic bits not supported')
        
        inputstate = pyQDD.makeZeroState(n_qubit)
        ddcirc = pyQDD.QuantumCircuit(n_qubit,8,4,inputstate)

        for i, qargs, cargs in circ.data:
            qiskit_gate_type = type(i)

            # filter out special cases first
            if qiskit_gate_type == Barrier:
                continue
            if (qiskit_gate_type == Measure) and (self.ignore_measure is True):
                continue
            assert(len(cargs)==0)

            if qiskit_gate_type in _supported_qiskit_gates:
                self.add_gate(i, qargs, cargs, ddcirc)
            else:
                # We assume the given Qiskit circuit has already been transpiled into a circuit of basis gates only.
                raise RuntimeError(f'Unsupported gate or instruction:'
                                   f' type={qiskit_gate_type.__name__}, name={i.name}.'
                                   f' It needs to transpile the circuit before evaluating it.')
            
        pyQDD.simulate(ddcirc)

        # Qiskit Backends (e.g., StateVectorSimulator) returns sampled counts in hex format.
        # E.g., {'0x001': 2, '0x100': 3}
        # On the other hand, Qulacs returns sampled counts as a list of decimal values.
        # E.g., [1, 4, 4, 1, 4]
        # So, we need to convert the format of Qulacs's sampled counts into the Qiskit's format.
        
        #hex_memory = list(map(hex, sampled_values))
        #hex_sampled_counts = Counter(hex_memory)
        result_data: Dict[str, Any] = {'counts': {'0x001':2, '0x100':3}} # TODO: FIX dummy data
        #if options['memory']:
        #    result_data['memory'] = hex_memory

        header = self._create_experiment_header(circ)
        result = {
            'success': True,
            'shots': options['shots'],
            'data': result_data,
            # Note: header information is used by several Qiskit functions; it must be contained in every result object.
            # E.g., header['memory_slots'] is used in Result#get_counts() for formatting sampled counts.
            'header': header,
        }
        return result
