import pytest

from qdd import QddProvider


def test_backend_config():
    backend = QddProvider().get_backend()
    config = backend.configuration()

    print(f'{config.backend_name=}')
    assert isinstance(config.backend_name, str)
    print(f'{config.backend_version=}')
    assert isinstance(config.backend_version, str)
    print(f'{config.n_qubits=}')
    assert isinstance(config.n_qubits, int)
    print(f'{config.open_pulse=}')
    assert isinstance(config.open_pulse, bool)
    print(f'{config.basis_gates=}')
    assert isinstance(config.basis_gates, list)

    with pytest.raises(Exception):
        print(config.dt)
    with pytest.raises(Exception):
        print(config.meas_levels)
    with pytest.raises(Exception):
        print(config.dtm)
    with pytest.raises(Exception):
        print(config.meas_map)
    with pytest.raises(Exception):
        print(config.drive(0))
    with pytest.raises(Exception):
        print(config.measure(0))
    with pytest.raises(Exception):
        print(config.acquire(0))


def test_backend_properties():
    backend = QddProvider().get_backend()
    props = backend.properties()
    assert props is None


def test_backend_defaults():
    backend = QddProvider().get_backend()
    with pytest.raises(Exception):
        backend.defaults()
