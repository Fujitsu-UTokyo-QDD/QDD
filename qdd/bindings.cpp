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

mEdge makeGate(QubitCount q, std::string name, Qubit target){
    return makeGate(q, gateMap[name], target);
}
mEdge makeGate(QubitCount q, std::string name, Qubit target, const Controls &c){
    return makeGate(q, gateMap[name], target, c);
}

PYBIND11_MODULE(pyQDD, m){
    py::class_<vEdge>(m, "vEdge");
    py::class_<mEdge>(m, "mEdge");
    m.def("makeZeroState", makeZeroState);
    m.def("mv_multiply", mv_multiply).def("mm_multiply", mm_multiply);

    // Gates
    m.def("makeGate", py::overload_cast<QubitCount, GateMatrix, Qubit>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, GateMatrix, Qubit, const Controls &>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, std::string, Qubit>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, std::string, Qubit, const Controls &>(&makeGate));
    m.def("RX", RX).def("RY", RY).def("RZ", RZ).def("CX", CX).def("SWAP", makeSwap);
    
}