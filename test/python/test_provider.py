import pytest
from qdd import QddProvider, QddBackend

def test_get_backend():
    backend = QddProvider().get_backend()
    #print(type(backend))
    assert type(backend) == QddBackend
