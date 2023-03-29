#include <pybind11/pybind11.h>
#include <map>
#include "dd.h"
#include "common.h"

namespace py = pybind11;

std::map<std::string, GateMatrix> gateMap{
    {"I", Imat},
    {"H", Hmat},
    {"X", Xmat},
    {"Y", Ymat},
    {"Z", Zmat},
    {"S", Smat},
    {"Sdag", Sdagmat},
    {"T", Tmat},
    {"Tdag", Tdagmat},
    {"SX", SXmat},
    {"SXdag", SXdagmat},
    {"V", Vmat},
    {"Vdag", Vdagmat}
};

std::mt19937_64 mt;

mEdge makeGate(QubitCount q, std::string name, Qubit target){
    return makeGate(q, gateMap[name], target);
}
mEdge makeGate(QubitCount q, std::string name, Qubit target, const Controls &c){
    return makeGate(q, gateMap[name], target, c);
}

std::string _measureAll(vEdge &rootEdge, bool collapse){
    return measureAll(rootEdge, collapse, mt);
}

char _measureOneCollapsing(vEdge &rootEdge, const Qubit index, const bool assumeProbabilityNormalization){
    return measureOneCollapsing(rootEdge, index, assumeProbabilityNormalization, mt);
}

PYBIND11_MODULE(pyQDD, m){
    py::class_<vEdge>(m, "vEdge").def("printVector",&vEdge::printVector).def("printVector_sparse",&vEdge::printVector_sparse);
    py::class_<mEdge>(m, "mEdge").def("printMatrix",&mEdge::printMatrix);
    m.def("makeZeroState", makeZeroState);
    m.def("mv_multiply", mv_multiply).def("mm_multiply", mm_multiply);

    // Gates
    m.def("makeGate", py::overload_cast<QubitCount, GateMatrix, Qubit>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, GateMatrix, Qubit, const Controls &>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, std::string, Qubit>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, std::string, Qubit, const Controls &>(&makeGate));
    m.def("RX", RX).def("RY", RY).def("RZ", RZ).def("CX", CX).def("SWAP", makeSwap);

    // Measure
    m.def("measureAll", _measureAll)
     .def("measureOneCollapsing", _measureOneCollapsing);
}