from qiskit.providers import Backend, JobStatus, JobV1
from qiskit.result import Result


class QddJob(JobV1):
    """A job issued by Qdd backend."""

    _async = False

    def __init__(self, backend, job_id, run_exp_fn, qobj):
        super().__init__(backend, job_id)
        self._run_exp_fn = run_exp_fn
        self._qobj = qobj
        self._result = None

    def cancel(self):
        raise RuntimeError('Cancel operation is not supported.')

    def submit(self):
        self._result = self._run_exp_fn(self._qobj, self.job_id())
        return self._result

    def result(self) -> Result:
        return self._result

    def status(self) -> JobStatus:
        if self._result is None:
            return JobStatus.RUNNING

        return JobStatus.DONE
