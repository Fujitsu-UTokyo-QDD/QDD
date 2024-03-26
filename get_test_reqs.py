import toml

data = toml.load(open("pyproject.toml"))

testlibs = data["tool"]["cibuildwheel"]["test-requires"]

# requirements_txt = open("requirements.txt", "w")
# print("\n".join(testlibs), file=requirements_txt)
# requirements_txt.close()
print("\n".join(testlibs))
