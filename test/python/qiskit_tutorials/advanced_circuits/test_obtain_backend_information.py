import pytest

from qdd import QddProvider


def test_backend_config():
    backend = QddProvider().get_backend()
    target = backend.target

    print(f"{backend.name=}")
    assert isinstance(backend.name, str)
    print(f"{backend.backend_version=}")
    assert isinstance(backend.backend_version, str)
    print(f"{backend.num_qubits=}")
    assert isinstance(backend.num_qubits, int)
    print(f"{target.instructions=}")
    assert isinstance(target.instructions, list)

    assert backend.dt is None
    with pytest.raises(Exception):
        print(backend.dtm)
    with pytest.raises(Exception):
        print(backend.meas_map)
    with pytest.raises(Exception):
        print(backend.drive_channel(0))
    with pytest.raises(Exception):
        print(backend.measure_channel(0))
    with pytest.raises(Exception):
        print(backend.acquire_channel(0))


def test_backend_properties():
    backend = QddProvider().get_backend()
    with pytest.raises(Exception):
        backend.properties()


def test_backend_defaults():
    backend = QddProvider().get_backend()
    with pytest.raises(Exception):
        backend.defaults()
