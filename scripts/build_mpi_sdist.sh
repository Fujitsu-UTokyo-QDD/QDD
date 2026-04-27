#!/usr/bin/env bash
set -euo pipefail

# Save original pyproject.toml and ensure restoration on exit
cp pyproject.toml pyproject.toml.bak
trap 'mv pyproject.toml.bak pyproject.toml' EXIT

# Verify the target line exists before attempting replacement
if ! grep -q '^name = "qdd"$' pyproject.toml; then
  echo 'ERROR: could not find `name = "qdd"` in pyproject.toml' >&2
  exit 1
fi

# Verify no conflicting section exists
if grep -q '^\[tool\.scikit-build\.cmake\.define\]$' pyproject.toml; then
  echo 'ERROR: [tool.scikit-build.cmake.define] already exists; aborting' >&2
  exit 1
fi

# Rewrite: rename project and enable MPI
sed -i 's/^name = "qdd"$/name = "qdd-mpi"/' pyproject.toml
cat >> pyproject.toml <<'EOF'

[tool.scikit-build.cmake.define]
isMPI = "ON"
EOF

# Build sdist
rm -rf dist-mpi
pip install build
python -m build --sdist --outdir dist-mpi/

echo ""
echo "qdd-mpi sdist built in dist-mpi/"
ls -la dist-mpi/
