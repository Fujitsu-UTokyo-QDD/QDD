# QDD Documentation

QDD is a decision diagram based quantum computing simulator that works as a
Qiskit backend. It can reduce memory usage compared with typical state-vector
simulators for circuits where the decision diagram representation remains
compact.

Use this documentation for installation, common usage, source builds, testing,
and API references.

```{toctree}
:maxdepth: 2
:caption: User Guide

installation
quickstart
usage
performance-characteristics
mpi
```

```{toctree}
:maxdepth: 2
:caption: Development

build-from-source
testing
development
internals/decision-diagram
internals/mpi
fx1000
```

```{toctree}
:maxdepth: 2
:caption: API Reference

api/python
api/cpp
```
