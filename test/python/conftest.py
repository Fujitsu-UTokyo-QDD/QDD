import os


if os.environ.get("QDD_TEST_USE_MPI") == "1":
    from qdd.qdd_backend import QddBackend

    _default_options = QddBackend._default_options.__func__

    def _default_options_with_mpi(cls):
        options = _default_options(cls)
        options.update_options(use_mpi=True)
        return options

    QddBackend._default_options = classmethod(_default_options_with_mpi)
