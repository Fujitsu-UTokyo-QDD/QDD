from qiskit.providers import Backend, JobStatus, JobV1
from qiskit.result import Result


class QddJob(JobV1):
    """Synchronous job returned by ``QddBackend.run``.

    QDD executes the supplied callable when the job is submitted and stores the
    resulting Qiskit ``Result`` object.
    """

    _async = False

    def __init__(self, backend, job_id, run_exp_fn, qobj):
        """Create a QDD job.

        Args:
            backend: Backend that created the job.
            job_id: Identifier reported through Qiskit's job API.
            run_exp_fn: Callable that performs the actual experiment run.
            qobj: Internal experiment payload passed to ``run_exp_fn``.
        """
        super().__init__(backend, job_id)
        self._run_exp_fn = run_exp_fn
        self._qobj = qobj
        self._result = None

    def cancel(self):
        raise RuntimeError("Cancel operation is not supported.")

    def submit(self):
        self._result = self._run_exp_fn(self._qobj, self.job_id())
        return self._result

    def result(self) -> Result:
        """Return the result produced by ``submit``."""
        return self._result

    def status(self) -> JobStatus:
        """Return ``RUNNING`` until a result exists, then ``DONE``."""
        if self._result is None:
            return JobStatus.RUNNING

        return JobStatus.DONE
