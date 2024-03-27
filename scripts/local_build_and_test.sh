#!/bin/bash

if [[ "$1" != "build" ]] && [[ "$1" != "test" ]]; then
    echo "You need to specify whether 'build' or 'test' as the first argument"
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
    if ! compgen -G "qdd/*.so" &> /dev/null; then
        echo ".so file does not exist. Building cmake"
        cmake . -DCMAKE_BUILD_TYPE=Release && cmake --build . -j
    fi
else
    echo "Make sure to execute this file from the project root directory"
    exit 1
fi && \

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
