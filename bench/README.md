# What is this repository?
This repository is for benchmarking ddsim (decision diagram quantum simulator) compared with Qiskit Aer.

# How to use?
## Dependency setup
At first, please install poetry.
You will need to add the installed directry to $PATH. Follow the instruction of the command.
```
curl -sSL https://install.python-poetry.org | python3 -
```

## Quantum Volume, VQE
If you want to try 5 qubits Quantum Volume circuit with Qiskit Aer, run the following command.
```
poetry run python QuantumVolume.py 5 Aer
```
If you want to try 3 qubits VQE with ddsim, run the following command.
```
poetry run python VQE.py 3 ddsim
```

## Shor
The input variable of `Shor.py` is not the number of qubits, but the value to be factored.
The following command try to factorize 15. (The result should be 3&5.) 
```
poetry run python Shor.py 15 ddsim
```

# Need help?
Send email to `yusuke-kimura@fujitsu.com`