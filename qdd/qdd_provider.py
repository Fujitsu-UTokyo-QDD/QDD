from qiskit.providers.providerutils import filter_backends

from .qdd_backend import QddBackend
from .qdd_sv_backend import QddSVBackend


class QddProvider:
    """Provider for QDD simulator backends.

    The provider exposes the QDD qasm simulator and statevector simulator using
    the standard Qiskit provider-style lookup methods.
    """

    def __init__(self, token=None):
        """Create a provider instance.

        Args:
            token: Reserved for compatibility with provider-style APIs. QDD
                does not currently use an authentication token.
        """
        self.token = token
        self._backends = [QddBackend(provider=self), QddSVBackend(provider=self)]

    def backends(self, name=None, **kwargs):
        """Return available backends matching the requested filters.

        Args:
            name: Optional backend name. Supported names are
                ``"qasm_simulator"`` and ``"statevector_simulator"``.
            **kwargs: Additional filters accepted by Qiskit's
                ``filter_backends`` helper.
        """
        if name:
            backend_candidates = [
                backend for backend in self._backends if backend.name == name
            ]
        else:
            backend_candidates = self._backends.copy()
        return filter_backends(backend_candidates, **kwargs)

    def get_backend(self, name=None, **kwargs) -> QddBackend:
        """Return a single backend by name.

        Args:
            name: Backend name. If omitted, ``"qasm_simulator"`` is used.
            **kwargs: Additional filters accepted by :meth:`backends`.

        Raises:
            RuntimeError: If no backend or multiple backends match the query.
        """
        if name == None:
            name = "qasm_simulator"
        backends = self.backends(name, **kwargs)
        if len(backends) == 0:
            raise RuntimeError(f"Backend {name} not found")
        if len(backends) > 1:
            raise RuntimeError(f"Multiple backends found for {name}")
        return self.backends(name, **kwargs)[0]
