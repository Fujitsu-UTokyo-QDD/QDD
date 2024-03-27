#cmake . -DCMAKE_BUILD_TYPE=Release
#cmake --build . -j
#python3 -m pip install build
#python3 -m build
python3 -m pip install toml; python3 -m pip install $(python3 get_test_reqs.py)
python3 -m pytest test
./test/qdd_test
