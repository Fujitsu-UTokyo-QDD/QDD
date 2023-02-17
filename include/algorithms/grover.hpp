#pragma once
#include "common.h"
#include <cstdlib>
#include <bitset>
#include "graph.hpp"
#include "task.h"


class Grover {
    public:
        explicit Grover(QubitCount q, int workers, int reduce,  std::size_t seed = 0);
        void full_grover();
    private:
        
        mEdge makeFullIteration();

        std::size_t iterations{1};
        std::size_t seed{0};
        std::string oracle;
        QubitCount n_qubits{0};
        QubitCount n_anciallae{1};
        std::mt19937_64 mt;
        int _nworkers;

        int _reduce;

};

vEdge grover(QubitCount n_qubits);
vEdge groverFiber(Scheduler& s ,QubitCount n_qubits);
