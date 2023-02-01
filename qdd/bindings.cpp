#include <pybind11/pybind11.h>
#include "python_helper.hpp"

namespace py = pybind11;
PYBIND11_MODULE(pyQDD, m){
    py::class_<vEdge>(m, "vEdge");
    py::class_<QuantumCircuit>(m, "QuantumCircuit")
        .def(py::init<uint32_t, int, int, vEdge>(), "Constructor");
    m.def("makeZeroState", makeZeroState);
    m.def("simulate", simulate);
    m.def("I", I);
    m.def("H", H);
    m.def("X", X);
    m.def("Y", Y);
    m.def("Z", Z);
    m.def("S", S);
    m.def("Sdag", Sdag);
    m.def("T", T);
    m.def("Tdag", Tdag);
    m.def("SX", SX);
    m.def("SXdag", SXdag);
    m.def("RX", RX);
    m.def("RY", RY);
    m.def("RZ", RZ);
    m.def("CX", CX);
    m.def("CY", CY);
    m.def("CZ", CZ);
    m.def("SWAP", SWAP);
}
