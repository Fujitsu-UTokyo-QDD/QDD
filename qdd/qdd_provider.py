from qiskit.providers import ProviderV1
from qiskit.providers.providerutils import filter_backends
from qdd.qdd_backend import QddBackend

class QddProvider(ProviderV1):
    
    def __init__(self):
        super().__init__()
        self._backends = [QddBackend(provider=self)]

    def backends(self, name=None, **kwargs):
        if name:
            backend_candidates = [backend for backend in self._backends if backend.name() == name]
        else:
            backend_candidates = self._backends.copy()
        return filter_backends(backend_candidates, **kwargs)

    def get_backend(self, name=None, **kwargs):
        return super().get_backend(name, **kwargs)