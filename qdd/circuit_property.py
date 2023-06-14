import dataclasses
from typing import Dict

from qiskit.circuit import Clbit, Qubit


@dataclasses.dataclass
class CircuitProperty:
    """This class holds circuit properties required for simulation.
    Attributes:
        stable_final_state (bool): whether the final state vector computed by evaluating this circuit is stable
            (unchanged) over multiple shots.
        clbit_final_values (Dict[Clbit, Qubit]): a dict from clbit to qubit. For each clbit, this dict holds the qubit
            last assigned to the clbit.
            For example, if the value of qubit1 is assigned to bit1 and then the value of qubit2 is assigned to bit1,
            this dict contains an entry of <clbit1, qubit2>.
    """

    stable_final_state: bool
    clbit_final_values: Dict[Clbit, Qubit]
