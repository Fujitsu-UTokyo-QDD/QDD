import pennylane as qml
from pennylane import numpy as np
import time
import sys

# Define the number of qubits and the number of layers
num_qubits = int(sys.argv[1])
num_layers = int(sys.argv[2])

# Define the device to be used for simulation
dev = qml.device("default.qubit", wires=num_qubits)

# Define the quantum circuit
@qml.qnode(dev)
def circuit():
    for layer in range(num_layers):
        # Apply a random gate to a random qubit in each layer
        gate = np.random.choice([qml.RX, qml.RY, qml.RZ, qml.CNOT])
        if gate == qml.CNOT:
            control, target = np.random.choice(num_qubits, size=2, replace=False)
            gate(wires=[control, target])
        else:
            qubit = np.random.choice(num_qubits)
            angle = np.random.uniform(low=-np.pi, high=np.pi)
            gate(angle, wires=qubit)
    
    return qml.probs(wires=range(num_qubits))

if __name__ == '__main__':
    print("qubits: "+ str(num_qubits)+", gates: "+str(num_layers))
    start_time = time.time()
    circuit()
    # print(circuit())
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Elapsed time: ", elapsed_time, "seconds")
