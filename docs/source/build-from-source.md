# Build From Source

Source builds are useful for development, local testing, and platforms where a
prebuilt wheel is not available.

## Requirements

- CMake 3.25 or newer
- C++17 compatible compiler
- Python build frontend

```sh
python -m pip install build
```

## Build C++ Targets and Wheel

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
python -m build
```

`CMAKE_BUILD_TYPE` may be `Release`, `Debug`, or `RelWithDebInfo`. The Python
wheel is written to the `dist` directory.

## Install the Built Wheel

Move to the environment where QDD should be installed and install the generated
wheel.

```sh
pip install /path/to/QDD/dist/qdd-*.whl
```

## Helper Script

The repository includes a helper script for local builds and tests:

```sh
./scripts/local_build_and_test_new.sh build
```

If no virtual environment is active, the script creates one under the project
root.
