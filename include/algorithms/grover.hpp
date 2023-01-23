#pragma once
#include "common.h"
#include <cstdlib>
#include <bitset>
#include "graph.hpp"


class Grover {
    public:
        explicit Grover(QubitCount q, int workers, std::size_t seed = 0);
        void full_grover();
    private:
        std::size_t iterations{1};
        std::size_t seed{0};
        std::string oracle;
        QubitCount n_qubits{0};
        QubitCount n_anciallae{1};
        QuantumCircuit qc;
        std::mt19937_64 mt;

};

