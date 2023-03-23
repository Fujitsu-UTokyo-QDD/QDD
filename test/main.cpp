#include "common.h"
#include "dd.h"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>

using namespace std::chrono_literals;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

#define M_PI 3.14159265358979323846

void benchmark(QubitCount q, Qubit threshold, int ngates) {
    typedef mEdge (*R)(QubitCount, int, float, Qubit);

    vEdge v = makeHybridZeroState(q, threshold);
    std::vector<mEdge> gates;
    std::vector<R> func{RX_Hybrid, RY_Hybrid, RZ_Hybrid};

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<> dist_gate(0, 2);
    std::uniform_int_distribution<> dist_target(0, q - 1);
    std::uniform_real_distribution<> dist_angle(-M_PI, M_PI);
    for (auto i = 0; i < ngates; i++) {
        gates.push_back(func[dist_gate(gen)](q, dist_target(gen),
                                             dist_angle(gen), threshold));
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < ngates; i++) {
        v = mv_multiply(gates[i], v);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << ms.count() << " milliseconds" << std::endl;
}

int main(int argc, char *argv[]) {
    QubitCount q = std::atoi(argv[1]);
    Qubit threshold = std::atoi(argv[2]);
    int ngates = std::atoi(argv[3]);
    std::cout << q << " qubits, " << threshold << " threshold," << ngates
              << " gates" << std::endl;
    benchmark(q, threshold, ngates);
}
