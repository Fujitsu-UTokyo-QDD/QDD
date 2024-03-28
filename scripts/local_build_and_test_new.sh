#!/bin/bash
# "./scripts/local_build_and_test.sh [arg]"

# The argument must be either "build" or "test".
# Either way, a shared library (.so file) will be built first. If the executed environment does not use virtualenv, it will create one named "venv-qdd".
# After this script, it is deactivated, but not deleted entirely.
# If "build" is specified, it will install build package and then build sdist and bdist.
# If "test" is specified, it will install toml package and test-requires packages in pyproject.toml and undergo pytest.

if [[ "$1" != "cmake" ]] && [[ "$1" != "build" ]] && [[ "$1" != "test" ]]; then
    echo "You need to specify among 'cmake', 'build', or 'test' as the first argument"
    exit 1
fi && \

# create a virtualenv if you are not in it
if [[ "$VIRTUAL_ENV" == "" ]]; then
    echo "You are not in a virtualenv. Creating a new virtualenv"
    python3 -m venv venv_qdd && source venv_qdd/bin/activate
    VENV_CREATED=1
else
    VENV_CREATED=0
fi && \

# .so build
if [ -d "qdd" ]; then # if qdd does not exist, it is not likely to be the project root, aborting
    cmake . -DCMAKE_BUILD_TYPE=Release && cmake --build . -j
else
    echo "Make sure to execute this file from the project root directory"
    exit 1
fi


if [[ "$1" == "build" ]]; then
    ( # build (sdist and wheel)
        if ! pip show build &> /dev/null; then
            python3 -m pip install build
        fi
    ) && python3 -m build
fi && \

if [[ "$1" == "test" ]]; then
    # Install toml, if not installed
    (
        if ! pip show toml &> /dev/null; then
            python3 -m pip install toml
        fi
    ) && \

    # Install requirements and execute pytest
    #python3 -m pip install .[test] && \
    python3 -m pip install $(python3 scripts/get_test_reqs.py) && \
    python3 -m pytest test && \
    ./test/qdd_test
fi

# Deactivate venv if created
echo $VENV_CREATED
if [[ $VENV_CREATED == 1 ]]; then
    echo "Deactivating virtualenv"
    deactivate
fi
