#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#ifdef isMPI
#include <boost/mpi/collectives.hpp>
#endif
#include "common.h"
#include "dd.h"
#ifdef isMT
#include "task.h"
#endif
namespace py = pybind11;

#ifdef isMT
bool mt_initialized = false;
int n_threads = 1;
Scheduler *s;
#endif

std::map<std::string, GateMatrix> gateMap{
    {"I", Imat},       {"H", Hmat},   {"X", Xmat},         {"Y", Ymat},
    {"Z", Zmat},       {"S", Smat},   {"Sdag", Sdagmat},   {"T", Tmat},
    {"Tdag", Tdagmat}, {"SX", SXmat}, {"SXdag", SXdagmat}, {"V", Vmat},
    {"Vdag", Vdagmat}};

std::random_device rd;
std::mt19937_64 mt(rd());

std::pair<vEdge, std::string> _measureAll(vEdge &rootEdge, bool collapse) {
    std::string result = measureAll(rootEdge, collapse, mt);
    return std::pair<vEdge, std::string>(rootEdge, result);
}

std::pair<vEdge, char> _measureOneCollapsing(vEdge &rootEdge,
                                             const Qubit index) {
    char result = measureOneCollapsing(rootEdge, index, mt);
    return std::pair<vEdge, char>(rootEdge, result);
}

std::pair<vEdge, double> _measureOne(vEdge &rootEdge, const Qubit index) {
    double result = measureOne(rootEdge, index, mt);
    return std::pair<vEdge, double>(rootEdge, result);
}

std::vector<std::complex<double>> _getVector(vEdge &rootEdge) {
    size_t dim;
    std_complex *vec = rootEdge.getVector(&dim);
    std::vector<std::complex<double>> result;
    for (int i = 0; i < dim; i++) {
        std::complex<double> tmp(vec[i].r, vec[i].i);
        result.push_back(tmp);
    }
    return result;
}

Controls get_controls(std::vector<Qubit> qs) {
    Controls controls;
    for (Qubit q : qs) {
        controls.emplace(Control{q, Control::Type::pos});
    }
    return controls;
}

mEdge makeGate(QubitCount q, std::string name, Qubit target) {
    return makeGate(q, gateMap[name], target);
}

mEdge makeGate(QubitCount q, GateMatrix m, Qubit target,
               const std::vector<Qubit> controls) {
    Controls c = get_controls(controls);
    return makeGate(q, m, target, c);
}

mEdge makeGate(QubitCount q, std::string name, Qubit target,
               const std::vector<Qubit> controls) {
    return makeGate(q, gateMap[name], target, controls);
}

mEdge makeTwoQubitGate(QubitCount q, TwoQubitGateMatrix m, Qubit target0,
                       Qubit target1, const std::vector<Qubit> controls) {
    Controls c = get_controls(controls);
    return makeTwoQubitGate(q, m, target0, target1, c);
}

mEdge applyGlobal(mEdge org, double angle){
    mEdge result = org;
    Complex e = {std::cos(angle), std::sin(angle)};
    result.w *= e;
    return result;
}

vEdge applyGlobal(vEdge org, double angle){
    vEdge result = org;
    Complex e = {std::cos(angle), std::sin(angle)};
    result.w *= e;
    return result;
}


#ifdef isMPI

boost::mpi::communicator _world;

std::vector<double> _probabilitiesMPI(const vEdge &v) {
    return probabilitiesMPI(_world, v);
}

vEdge _makeZeroStateMPI(QubitCount q) { return makeZeroStateMPI(q, _world); }
vEdge _makeOneStateMPI(QubitCount q) { return makeOneStateMPI(q, _world); }
void _printVectorMPI(vEdge &v) {
    _world.barrier();
    v.printVectorMPI(_world);
}
vEdge _mv_multiply_MPI(mEdge lhs, vEdge rhs, std::size_t total_qubits,
                       std::size_t largest_qubit) {
    return mv_multiply_MPI(lhs, rhs, _world, total_qubits, largest_qubit);
}
vEdge _mv_multiply_MPI_bcast(mEdge lhs, vEdge rhs, std::size_t total_qubits,
                             std::size_t largest_qubit) {
    return mv_multiply_MPI_bcast3(lhs, rhs, _world, total_qubits,
                                  largest_qubit);
}
std::pair<vEdge, std::string> _measureAllMPI(vEdge &rootEdge,
                                             const bool collapse) {
    std::string result = measureAllMPI(_world, rootEdge, collapse, mt);
    return std::pair<vEdge, std::string>(rootEdge, result);
}

std::pair<vEdge, char> _measureOneCollapsingMPI(vEdge &rootEdge,
                                                const Qubit index,
                                                const Qubit n_qubits) {
    char result =
        measureOneCollapsingMPI(_world, rootEdge, index, n_qubits, mt, true);
    return std::pair<vEdge, char>(rootEdge, result);
}

std::pair<vEdge, double> _measureOneMPI(vEdge &rootEdge, const Qubit index,
                                        const Qubit n_qubits) {
    double result = measureOneMPI(_world, rootEdge, index, n_qubits, mt);
    return std::pair<vEdge, double>(rootEdge, result);
}

void dump_mpi() {
    std::cout << _world.rank() << "/" << _world.size() << std::endl;
}

std::vector<std::complex<double>> _getVectorMPI(vEdge &edge) {
    size_t dim;
    std_complex *vec = edge.getVector(&dim);
    std::vector<std_complex> result;
    for (int i = 0; i < dim; i++) {
        result.push_back(vec[i]);
    }

    std::vector<std::vector<std_complex>> all_results;
    bmpi::all_gather(_world, result, all_results);

    std::vector<std::complex<double>> final_result;
    for (int i = 0; i < all_results.size(); i++) {
        for (int j = 0; j < all_results[i].size(); j++) {
            final_result.push_back(
                std::complex(all_results[i][j].r, all_results[i][j].i));
        }
    }
    return final_result;
}
#endif

#ifdef isMT
void terminate_mt() {
    if (mt_initialized) {
        mt_initialized = false;
        n_threads = 1;
        delete s;
    }
    return;
}

int init_mt(int n = 8) {
    if (n != n_threads) {
        terminate_mt();
        if (n > 1) {
            mt_initialized = true;
            n_threads = n;
            s = new Scheduler(n - 1);
        }
    }
    return n_threads;
}
#endif

std::vector<double> _probabilities(const vEdge &rootEdge) {
    std::vector<double> result = probabilities(rootEdge);
    return result;
}

PYBIND11_MODULE(pyQDD, m) {
    py::class_<vEdge>(m, "vEdge")
        .def("getEigenVector", &vEdge::getEigenVector)
        .def("printVector", &vEdge::printVector)
        .def("printVector_sparse", &vEdge::printVector_sparse)
#ifdef isMPI
        .def("printVectorMPI", _printVectorMPI)
#endif
        ;
    py::class_<mEdge>(m, "mEdge")
        .def("printMatrix", &mEdge::printMatrix)
        .def("getEigenMatrix", &mEdge::getEigenMatrix);
    m.def("makeZeroState", makeZeroState);
    m.def("mv_multiply", mv_multiply).def("mm_multiply", mm_multiply);
    m.def("get_nNodes", get_nNodes)
        .def("gc", gc)
        .def("gc_mat", gc_mat)
        .def("set_gc_thr", set_gc_thr);
    m.def("applyGlobal", py::overload_cast<mEdge, double>(&applyGlobal))
     .def("applyGlobal", py::overload_cast<vEdge, double>(&applyGlobal));

    // Gates
    m.def("makeGate",
          py::overload_cast<QubitCount, GateMatrix, Qubit>(&makeGate))
        .def("makeGate", py::overload_cast<QubitCount, GateMatrix, Qubit,
                                           const std::vector<Qubit>>(&makeGate))
        .def("makeGate",
             py::overload_cast<QubitCount, std::string, Qubit>(&makeGate))
        .def("makeGate",
             py::overload_cast<QubitCount, std::string, Qubit,
                               const std::vector<Qubit>>(&makeGate));
    m.def("RX", RX).def("RY", RY).def("RZ", RZ).def("CX", CX);
    m.def("rxmat", rx)
        .def("rymat", ry)
        .def("rzmat", rz)
        .def("u1", u1)
        .def("u2", u2)
        .def("u3", u3)
        .def("u", u)
        .def("p", p)
        .def("r", r);

    m.def("makeTwoQubitGate",
          py::overload_cast<QubitCount, TwoQubitGateMatrix, Qubit, Qubit>(
              &makeTwoQubitGate))
        .def("makeTwoQubitGate",
             py::overload_cast<QubitCount, TwoQubitGateMatrix, Qubit, Qubit,
                               const std::vector<Qubit>>(&makeTwoQubitGate));
    m.def("RXX", RXX)
        .def("RYY", RYY)
        .def("RZZ", RZZ)
        .def("RZX", RZX)
        .def("SWAP", SWAP)
        .def("ISWAP", ISWAP)
        .def("CSWAP", CSWAP);
    m.def("rxxmat", rxx_matrix)
        .def("ryymat", ryy_matrix)
        .def("rzzmat", rzz_matrix)
        .def("rzxmat", rzx_matrix)
        .def("swapmat", swap_matrix)
        .def("iswapmat", iswap_matrix);

    m.def("unitary",
          py::overload_cast<QubitCount, ComplexMatrix &>(makeLargeGate))
        .def("unitary",
             py::overload_cast<QubitCount, ComplexMatrix &,
                               const std::vector<Qubit> &>(makeLargeGate));

    // Measure
    m.def("measureAll", _measureAll)
        .def("measureOneCollapsing", _measureOneCollapsing);
    m.def("getVector", _getVector);
    m.def("probabilities", _probabilities);
    m.def("measureOne", _measureOne);

    // Others
    m.def("genDot", py::overload_cast<vEdge &>(&genDot))
        .def("genDot", py::overload_cast<mEdge &>(&genDot));

#ifdef isMPI
    // MPI
    m.def("dump_mpi", dump_mpi)
        .def("probabilitiesMPI", _probabilitiesMPI)
        .def("makeZeroStateMPI", _makeZeroStateMPI)
        .def("makeOneStateMPI", _makeOneStateMPI)
        .def("mv_multiply_MPI", _mv_multiply_MPI)
        .def("mv_multiply_MPI_bcast", _mv_multiply_MPI_bcast)
        .def("measureAllMPI", _measureAllMPI)
        .def("getVectorMPI", _getVectorMPI)
        .def("save_binary", save_binary)
        .def("load_binary", load_binary)
        .def("measureOneCollapsingMPI", _measureOneCollapsingMPI)
        .def("measureOneMPI", _measureOneMPI);
#endif
#ifdef isMT
    m.def("initMT", init_mt).def("terminateMT", terminate_mt);
#endif
}
