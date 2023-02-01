#pragma once
#include "graph.hpp"

void I(int target, QuantumCircuit &qc){ qc.emplace_back(Imat, target); }
void H(int target, QuantumCircuit &qc){ qc.emplace_back(Hmat, target); }
void X(int target, QuantumCircuit &qc){ qc.emplace_back(Xmat, target); }
void Y(int target, QuantumCircuit &qc){ qc.emplace_back(Ymat, target); }
void Z(int target, QuantumCircuit &qc){ qc.emplace_back(Zmat, target); }
void S(int target, QuantumCircuit &qc){ qc.emplace_back(Smat, target); }
void Sdag(int target, QuantumCircuit &qc){ qc.emplace_back(Sdagmat, target); }
void T(int target, QuantumCircuit &qc){ qc.emplace_back(Tmat, target); }
void Tdag(int target, QuantumCircuit &qc){ qc.emplace_back(Tdagmat, target); }
void SX(int target, QuantumCircuit &qc){ qc.emplace_back(SXmat, target); }
void SXdag(int target, QuantumCircuit &qc){ qc.emplace_back(SXdagmat, target); }
void RX(int target, float angle, QuantumCircuit &qc) {
    std::complex<float> i1 = {std::cos(angle / 2), 0};
    std::complex<float> i2 = {0, -std::sin(angle / 2)};
    qc.emplace_back(GateMatrix{i1,i2,i2,i1}, target);
}
void RY(int target, float angle, QuantumCircuit &qc) {
    std::complex<float> i1 = {std::cos(angle / 2), 0};
    std::complex<float> i2 = {-std::sin(angle / 2), 0};
    std::complex<float> i3 = {std::sin(angle / 2), 0};
    qc.emplace_back(GateMatrix{i1,i2,i3,i1}, target);
}
void RZ(int target, float angle, QuantumCircuit &qc) {
    std::complex<float> i1 = {std::cos(angle / 2), -std::sin(angle / 2)};
    std::complex<float> i2 = {std::cos(angle / 2), std::sin(angle / 2)};
    qc.emplace_back(GateMatrix{i1,cf_zero,cf_zero,i1}, target);
}
void CX(int target, int control, QuantumCircuit &qc){
    Controls controls;
    controls.emplace(Control{control, Control::Type::pos});
    qc.emplace_back(Xmat, target, controls);
}
void CY(int target, int control, QuantumCircuit &qc){
    Controls controls;
    controls.emplace(Control{control, Control::Type::pos});
    qc.emplace_back(Ymat, target, controls);
}
void CZ(int target, int control, QuantumCircuit &qc){
    Controls controls;
    controls.emplace(Control{control, Control::Type::pos});
    qc.emplace_back(Zmat, target, controls);
}
void SWAP(int t0, int t1, QuantumCircuit &qc) { qc.emplace_gate(makeSwap(qc.getQubits(), t0, t1)); }

void simulate(QuantumCircuit &circuit){
    circuit.buildCircuit();
    vEdge result = circuit.wait().vectorResult();
    result.printVector();
    return;
}
