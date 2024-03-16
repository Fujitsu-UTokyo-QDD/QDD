# -*- coding: utf-8 -*-
import os
from pathlib import Path
from setuptools import Extension, setup

packages = \
['qdd']

package_data = \
{'': ['*'], 'qdd': ['CMakeFiles/*', 'CMakeFiles/pyQDD.dir/*']}

install_requires = \
['qiskit-aer>=0.13.3']

class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())

setup_kwargs = {
    'name': 'qdd',
    'version': '0.1.0',
    'description': 'Qiskit Provider for QDD',
    'long_description': "# QDD\n## About\nQDD is a decision diagram based quantum computing simulator.\nMemory usage can be significantly reduced compared to typical state vector based simulators.\n\nIt works as a Qiskit backend.\nAt this point, we do not distribute wheel file, so you need to build it locally.\n\n## License\nBSD 3-Clause Clear License\n\n## Build\nYou must prepare the following software.\n\n* cmake (>=3.25)\n  * You can download the executable from [Official GitHub](https://github.com/Kitware/CMake/releases).\n  * If you use Ubuntu, [the Kitware repository](https://apt.kitware.com/) is easier to install.\n* poetry\n  * See [the installation manual](https://python-poetry.org/docs/#installing-with-the-official-installer).\n\nYou can build QDD as follows.\n```sh\n$ cmake . -DCMAKE_BUILD_TYPE=Release\n$ cmake --build . -j\n$ poetry build\n```\nBUILD_TYPE can be either `Release`, `Debug` or `RelWithDebInfo`.\nYou can find an installable wheel file in 'dist' directory.\n\n## Instllation\nPlease move to your working directory and install the wheel file created in the build phase.\n```sh\n$ pip install {QDD_DIR}/dist/qdd-XXX.whl\n```\n\n## Usage\nQDD works as a Qiskit backend.\n\n```py\nfrom qiskit import QuantumCircuit, execute\nfrom qdd import QddBackend, QddProvider\n\nbackend = QddProvider().get_backend()\ncirc = QuantumCircuit(3)\ncirc.h(0)\ncirc.cx(0,1)\ncirc.measure_all()\nqdd_job = execute(circ, backend=backend, shots=1000)\nprint(qdd_job.result())\n```\n\nIf you need statevector, create the backend as follows.\n```py\nbackend = QddProvider().get_backend('statevector_simulator')\n```\n\n# Limitation of Liability\nIn no event and under no legal theory, whether in tort (including negligence), contract, or otherwise, unless required by applicable law (such as deliberate and grossly negligent acts) or agreed to in writing, shall any Contributor be liable to You for damages, including any direct, indirect, special, incidental, or consequential damages of any character arising as a result of this License or out of the use or inability to use the Work (including but not limited to damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses), even if such Contributor has been advised of the possibility of such damages.\n",
    'author': 'Yusuke Kimura',
    'author_email': 'yusuke-kimura@fujitsu.com',
    'maintainer': 'None',
    'maintainer_email': 'None',
    'url': 'None',
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'python_requires': '>=3.8,<3.11',
    'ext_modules':[CMakeExtension("qdd/pyQDD", "build/qdd")],
}

setup(**setup_kwargs)
