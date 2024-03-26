# A python script to obtain requirements for test from pyproject.toml (and install them all)
# python3 -m pip install toml; python3 -m pip install $(python3 get_test_reqs.py)
import toml

data = toml.load(open("pyproject.toml"))

testlibs = data["tool"]["cibuildwheel"]["test-requires"]

# requirements_txt = open("requirements.txt", "w")
# print("\n".join(testlibs), file=requirements_txt)
# requirements_txt.close()
print("\n".join(testlibs))
