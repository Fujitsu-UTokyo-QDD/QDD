from qiskit import QuantumCircuit
from qdd import QddProvider
import pytest

@pytest.mark.mpi
def test_mpi():

    backend = QddProvider().get_backend()
    backend.set_options(use_mpi=True)

    circ = QuantumCircuit(3)
    circ.h(0)
    circ.cx(0,1)
    circ.measure_all()

    job = backend.run(circ)
    print(job.result().results[0].data.counts)

if __name__=="__main__":
    test_mpi()
