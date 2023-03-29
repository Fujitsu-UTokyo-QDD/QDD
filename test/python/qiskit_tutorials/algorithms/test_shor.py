import math
import time

import pytest
from qiskit import Aer
from qiskit.algorithms import Shor
from qiskit.utils import QuantumInstance

from qdd import QddProvider


@pytest.mark.slow
def test_shor():
    """Tests the simulation behavior of the Shor's factorization algorithm"""

    num_to_factor = 15
    print(f'Factoring {num_to_factor}')
    print(f'#qubits required is {4 * math.ceil(math.log(num_to_factor, 2)) + 2}.')

    # use Qdd backend
    real_time, cpu_time = time.perf_counter(), time.process_time()
    backend = QddProvider().get_backend()
    qi = QuantumInstance(backend, shots=1024, seed_transpiler=50, seed_simulator=80)
    shor = Shor(quantum_instance=qi)
    result_qdd = shor.factor(num_to_factor)
    real_time_spent, cpu_time_spent = time.perf_counter() - real_time, time.process_time() - cpu_time
    print(f'{result_qdd.factors}'
          f' (solved via Qdd in real={real_time_spent:.3f}(sec), cpu={cpu_time_spent:.3f}(sec))')

    # use Aer backend
    real_time, cpu_time = time.perf_counter(), time.process_time()
    backend = Aer.get_backend('aer_simulator')
    qi = QuantumInstance(backend, shots=1024, seed_transpiler=50, seed_simulator=80)
    shor = Shor(quantum_instance=qi)
    result_aer = shor.factor(num_to_factor)
    real_time_spent, cpu_time_spent = time.perf_counter() - real_time, time.process_time() - cpu_time
    print(f'{result_aer.factors} (solved via Aer in real={real_time_spent:.3f}(sec), cpu={cpu_time_spent:.3f}(sec))')

    # validate results
    assert sorted(map(sorted, result_qdd.factors)) == sorted(map(sorted, result_aer.factors))
