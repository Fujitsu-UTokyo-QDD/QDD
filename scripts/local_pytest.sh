#!/bin/bash

# create a virtualenv if you are not in it
(
    if [[ "$VIRTUAL_ENV" != ""]]
    then
        echo "You are not in a virtualenv. Creating a new virtualenv"
        python3 -m venv venv_qdd && source venv_qdd/bin/activate
    fi
) && \

# create .so file, if any of them does not exist
(
    if [ -d "qdd" ] # if qdd does not exist, it is not likely to be the project root, aborting
    then
        if ! compgen -G "qdd/*.so" &> /dev/null
        then
            echo ".so file does not exist. Building cmake."
            cmake . -DCMAKE_BUILD_TYPE=Release && cmake --build . -j
        fi
    else
        echo "Make sure to execute this file from the project root directory"
        exit 1
    fi
) && \

# Install toml, if not installed
(if ! pip show toml &> /dev/null; then
    python3 -m pip install toml
fi) && \

# Install requirements and execute pytest
python3 -m pip install $(python3 scripts/get_test_reqs.py) && \
python3 -m pytest test && \
./test/qdd_test
