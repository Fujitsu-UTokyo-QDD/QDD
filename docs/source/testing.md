# Testing

QDD includes Python tests and C++ tests.

## Python Tests

Install the package with test dependencies:

```sh
pip install -e .[test]
pytest test/python
```

The default pytest configuration skips slow and MPI tests.

```sh
pytest -m "not slow and not mpi"
```

Avoid installing an editable checkout into a long-lived user environment unless
you intend to develop QDD. Editable installs make Python import from the working
tree.

## C++ Tests

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBINDINGS=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Helper Script

```sh
./scripts/local_build_and_test_new.sh test
```

This command builds the C++ extension and runs the configured test suite.
