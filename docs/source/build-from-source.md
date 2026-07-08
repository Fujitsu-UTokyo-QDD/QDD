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

## Windows Build with clang-cl

Windows builds use clang-cl (LLVM 18 or later) to target the MSVC ABI on
x86_64 and ARM64.
`cl.exe` is not used.

### Install LLVM

```powershell
choco install llvm
```

Add `C:\Program Files\LLVM\bin` to `PATH`.

### Install Ninja

```sh
pip install ninja
```

### Configure and Build

```powershell
cmake -S . -B build -G Ninja `
  -DCMAKE_C_COMPILER=clang-cl `
  -DCMAKE_CXX_COMPILER=clang-cl `
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Run C++ Tests

```sh
ctest --test-dir build --output-on-failure -C Release
```

### Build and Install Python Wheel

```powershell
$env:CMAKE_ARGS = "-G Ninja -DCMAKE_C_COMPILER=clang-cl -DCMAKE_CXX_COMPILER=clang-cl"
pip install ".[test]"
pytest test/python -m "not slow and not mpi"
```

> **Note:** MPI is not supported on Windows.

## macOS Build

macOS builds use AppleClang. Both Intel (`x86_64`) and Apple Silicon (`arm64`)
are built as native, architecture-specific wheels rather than Universal2
wheels.

Release wheels target macOS 10.13 or newer on Intel and macOS 11.0 or newer on
Apple Silicon.

### Install Build Tools

```sh
brew install cmake
python -m pip install --upgrade pip build
```

### Configure and Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Run C++ Tests

```sh
ctest --test-dir build --output-on-failure
```

### Build and Test Python Package

```sh
pip install ".[test]"
pytest test/python -m "not slow and not mpi"
```

> **Note:** MPI is not supported on macOS.
