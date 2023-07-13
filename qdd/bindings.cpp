#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
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

std::pair<vEdge, std::string> _measureAll(vEdge &rootEdge, bool collapse){
    std::string result = measureAll(rootEdge, collapse, mt);
    return std::pair<vEdge, std::string>(rootEdge, result);
}



std::pair<vEdge, char> _measureOneCollapsing(vEdge &rootEdge, const Qubit index){
    char result = measureOneCollapsing(rootEdge, index, mt);
    return std::pair<vEdge, char>(rootEdge, result);
}



std::vector<std::complex<double>> _getVector(vEdge &rootEdge){
    size_t dim;
    std_complex *vec = rootEdge.getVector(&dim);
    std::vector<std::complex<double>> result;
    for (int i = 0; i < dim;i++){
        std::complex<double> tmp(vec[i].r, vec[i].i);
        result.push_back(tmp);
    }
    return result;
}

Controls get_controls(std::vector<Qubit> qs){
    Controls controls;
    for(Qubit q: qs){
        controls.emplace(Control{q, Control::Type::pos});
    }
    return controls;
}

mEdge makeControlGate(QubitCount q, std::string name, Qubit target, const std::vector<Qubit> controls){
    Controls c = get_controls(controls);
    return makeGate(q, name, target, c);
}

PYBIND11_MODULE(pyQDD, m){
    py::class_<vEdge>(m, "vEdge").def("printVector",&vEdge::printVector).def("printVector_sparse",&vEdge::printVector_sparse);
    py::class_<mEdge>(m, "mEdge").def("printMatrix",&mEdge::printMatrix).def("getEigenMatrix", &mEdge::getEigenMatrix);
    m.def("makeZeroState", makeZeroState);
    m.def("mv_multiply", mv_multiply).def("mm_multiply", mm_multiply);

    // Gates
    m.def("makeGate", py::overload_cast<QubitCount, GateMatrix, Qubit>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, GateMatrix, Qubit, const Controls &>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, std::string, Qubit>(&makeGate))
     .def("makeGate", py::overload_cast<QubitCount, std::string, Qubit, const Controls &>(&makeGate))
     .def("makeControlGate", makeControlGate);
    m.def("RX", RX).def("RY", RY).def("RZ", RZ).def("CX", CX).def("SWAP", makeSwap);
    m.def("rxmat", rx).def("rymat", ry).def("rzmat", rz).def("u1", u1).def("u2", u2).def("u3", u3).def("u", u).def("p", p).def("r", r);

    // Measure
    m.def("measureAll", _measureAll)
     .def("measureOneCollapsing", _measureOneCollapsing);
    m.def("getVector", _getVector);
}