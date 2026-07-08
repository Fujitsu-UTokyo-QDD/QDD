cd /tmp
rm -rf test-qdd
mkdir test-qdd
cd test-qdd

python -m venv .venv
source .venv/bin/activate
pip install -e ~/QDD[test]

cp -r ~/QDD/test test
cp ~/QDD/pyproject.toml .

pytest

