from qiskit.providers import Backend, JobStatus, JobV1
from qiskit.result import Result


class QddFailedJob(JobV1):
    """A job representing a failure of circuit evaluation."""

    def __init__(self, backend: Backend, job_id: str, result: Result):
        super().__init__(backend, job_id)
        self._result = result

    def submit(self):
        raise RuntimeError('Submit operation is not supported.')

    def result(self) -> Result:
        return self._result

    def cancel(self):
        raise RuntimeError('Cancel operation is not supported.')

    def status(self) -> JobStatus:
        return JobStatus.ERROR
