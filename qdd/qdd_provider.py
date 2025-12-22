from qiskit.providers.providerutils import filter_backends

from .qdd_backend import QddBackend
from .qdd_sv_backend import QddSVBackend


class QddProvider:

    def __init__(self, token=None):
        self.token = token
        self._backends = [QddBackend(provider=self), QddSVBackend(provider=self)]

    def backends(self, name=None, **kwargs):
        if name:
            backend_candidates = [
                backend for backend in self._backends if backend.name == name
            ]
        else:
            backend_candidates = self._backends.copy()
        return filter_backends(backend_candidates, **kwargs)

    def get_backend(self, name=None, **kwargs) -> QddBackend:
        if name == None:
            name = "qasm_simulator"
        backends = self.backends(name, **kwargs)
        if len(backends) == 0:
            raise RuntimeError(f"Backend {name} not found")
        if len(backends) > 1:
            raise RuntimeError(f"Multiple backends found for {name}")
        return self.backends(name, **kwargs)[0]
