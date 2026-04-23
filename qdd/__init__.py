try:
    from ._version import version as __version__
except ImportError:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:
        __version__ = "0.0.0+unknown"
    else:
        try:
            __version__ = version("qdd")
        except PackageNotFoundError:
            __version__ = "0.0.0+unknown"

from .qdd_backend import QddBackend
from .qdd_job import QddJob
from .qdd_provider import QddProvider

__all__ = [QddBackend.__name__, QddProvider.__name__, QddJob.__name__]
