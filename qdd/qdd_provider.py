from qiskit.providers import ProviderV1 as Provider
from qiskit.providers.providerutils import filter_backends

from .qdd_backend import QddBackend

class QddProvider(Provider):

    def __init__(self, token=None):
        super().__init__()
        self.token = token
        self.backends = [QddBackend(provider=self)]

    def backends(self, name=None, **kwargs):
        if name:
            backends = [
                backend for backend in backends if backend.name() == name]
        return filter_backends(backends, **kwargs)
