# QDD

## About

QDD is a decision diagram based quantum computing simulator.
Memory usage can be significantly reduced compared to typical state vector based simulators.

It works as a Qiskit backend.
At this point, we do not distribute wheel file, so you need to build it locally.

## License

BSD 3-Clause Clear License

## Install

Supported environment

- Linux x86_64 with Python3.8 ~ 3.10

Command

```sh
$ pip install qdd
```

## Install from source

### Build

You must prepare the following software.

- cmake (>=3.25)
  - You can download the executable from [Official GitHub](https://github.com/Kitware/CMake/releases).
  - If you use Ubuntu, [the Kitware repository](https://apt.kitware.com/) is easier to install.
- build
  - Just install it via pip.

You can build QDD as follows.

```sh
$ cmake . -DCMAKE_BUILD_TYPE=Release
$ cmake --build . -j
$ python3 -m pip install build
$ python3 -m build
```

BUILD_TYPE can be either `Release`, `Debug` or `RelWithDebInfo`.
You can find an installable wheel file in 'dist' directory.

Or, you can just execute 'scripts/local_build_and_test.sh build' in the project root to build the package automatically.

### Wheel Instllation

Please move to your working directory and install the wheel file created in the build phase.

```sh
$ pip install {QDD_DIR}/dist/qdd-XXX.whl
```

## Usage

QDD works as a Qiskit backend.

```py
from qiskit import QuantumCircuit
from qiskit.primitives import BackendSampler

from qdd import QddProvider

backend = QddProvider().get_backend()
circ = QuantumCircuit(3)
circ.h(0)
circ.cx(0,1)
circ.measure_all()
sampler = BackendSampler(backend=backend)
qdd_job = sampler.run(circuits=circ)
print(qdd_job.result())
```

QDD has its own implementation of the Sampler class

```py
from qdd.qdd_sampler import Sampler
sampler = Sampler()
qdd_job = sampler.run(circuits=circ)
print(qdd_job.result())
```

If you need statevector, create the backend as follows.

```py
backend = QddProvider().get_backend('statevector_simulator')
```

# Limitation of Liability

In no event and under no legal theory, whether in tort (including negligence), contract, or otherwise, unless required by applicable law (such as deliberate and grossly negligent acts) or agreed to in writing, shall any Contributor be liable to You for damages, including any direct, indirect, special, incidental, or consequential damages of any character arising as a result of this License or out of the use or inability to use the Work (including but not limited to damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses), even if such Contributor has been advised of the possibility of such damages.
