from qiskit.providers import BackendV1, JobV1, Options, Provider
from qiskit.providers.models import BackendConfiguration
from qdd.qdd_backend import QddBackend

class QddSVBackend(QddBackend):
    def __init__(self, provider: Provider, configuration=None):
        self._save_SV = True
        conf_dict = self._DEFAULT_CONFIG.copy()
        conf_dict['backend_name'] = 'statevector_simulator'
        conf_dict['n_qubits'] = 20
        if configuration == None:
            configuration = BackendConfiguration.from_dict(conf_dict)
        
        super().__init__(
            configuration=configuration,
            provider=provider)