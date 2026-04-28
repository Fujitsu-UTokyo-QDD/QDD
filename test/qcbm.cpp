#include <iostream>
#include <random>

#include "cache.hpp"
#include "table.hpp"

std::random_device seed_gen;
std::mt19937_64 mt(seed_gen());
std::uniform_real_distribution<double> dist(0.0, 1.0L);

double get_random() { return dist(mt); }

vEdge exec(int nQubits) {
    int n_SingleGates = nQubits * 28;
    int n_CXGates = nQubits * 9;
    std::cout << "nQubits=" << nQubits
              << " Total Gates=" << n_SingleGates + n_CXGates << std::endl;

    vEdge v = makeZeroState(nQubits);

    // first_rotation
    double angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RX(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RZ(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    std::cout << "First rotation" << std::endl;

    // entangler
    for (int i = 0; i < nQubits; i++) {
        int control = i;
        int target = (i + 1) % nQubits;
        auto g = CX(nQubits, target, control);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    std::cout << "Entangler" << std::endl;

    for (int k = 0; k < 8; k++) {
        // mid rotation
        angle = get_random();
        for (int target = 0; target < nQubits; target++) {
            auto g = RZ(nQubits, target, angle);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "mid rotation (1) " << k << std::endl;
        angle = get_random();
        for (int target = 0; target < nQubits; target++) {
            auto g = RX(nQubits, target, angle);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "mid rotation (2) " << k << std::endl;
        angle = get_random();
        for (int target = 0; target < nQubits; target++) {
            auto g = RZ(nQubits, target, angle);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "mid rotation (3) " << k << std::endl;
        // entangler
        for (int i = 0; i < nQubits; i++) {
            int control = i;
            int target = (i + 1) % nQubits;
            auto g = CX(nQubits, target, control);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "entangler " << k << std::endl;
    }
    // last rotation
    angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RZ(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RX(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    return v;
}

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();
    auto v = exec(std::atoi(argv[1]));

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> ms = end - start;
    std::cout << "nQubits " << std::atoi(argv[1]) << " nNodes " << get_nNodes(v)
              << " " << ms.count() / 1000000 << " sec" << std::endl;
}