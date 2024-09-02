#include <iostream>
#include <random>

#include "cache.hpp"
#include "table.hpp"

std::random_device seed_gen;
std::mt19937_64 mt(seed_gen());
std::uniform_real_distribution<double> dist(0.0, 1.0L);

// int GC_SIZE = 131072*2*2*2*2;

double get_random() { return dist(mt); }

vEdge exec(int nQubits, bmpi::communicator& world) {
    int n_SingleGates = nQubits * 28;
    int n_CXGates = nQubits * 9;
    if (world.rank() == 0)
        std::cout << "nQubits=" << nQubits
                  << " Total Gates=" << n_SingleGates + n_CXGates << std::endl;

    vEdge v = makeZeroStateMPI(nQubits, world);

    // first_rotation
    double angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RX(nQubits, target, angle);
        v = mv_multiply_MPI(g, v, world, nQubits, target);
        v = gc(v);
    }
    angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RZ(nQubits, target, angle);
        v = mv_multiply_MPI(g, v, world, nQubits, target);
        v = gc(v);
    }

    // entangler
    for (int i = 0; i < nQubits; i++) {
        int control = i;
        int target = (i + 1) % nQubits;
        auto g = CX(nQubits, target, control);
        v = mv_multiply_MPI(g, v, world, nQubits,
                            target > control ? target : control);
        v = gc(v);
    }

    for (int k = 0; k < 8; k++) {
        // mid rotation
        angle = get_random();
        for (int target = 0; target < nQubits; target++) {
            auto g = RZ(nQubits, target, angle);
            v = mv_multiply_MPI(g, v, world, nQubits, target);
            v = gc(v);
        }
        angle = get_random();
        for (int target = 0; target < nQubits; target++) {
            auto g = RX(nQubits, target, angle);
            v = mv_multiply_MPI(g, v, world, nQubits, target);
            v = gc(v);
        }
        angle = get_random();
        for (int target = 0; target < nQubits; target++) {
            auto g = RZ(nQubits, target, angle);
            v = mv_multiply_MPI(g, v, world, nQubits, target);
            v = gc(v);
        }
        // entangler
        for (int i = 0; i < nQubits; i++) {
            int control = i;
            int target = (i + 1) % nQubits;
            auto g = CX(nQubits, target, control);
            v = mv_multiply_MPI(g, v, world, nQubits,
                                target > control ? target : control);
            v = gc(v);
        }
        if (world.rank() == 0)
            std::cout << k + 1 << " th iteration" << std::endl;
    }
    // last rotation
    angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RZ(nQubits, target, angle);
        v = mv_multiply_MPI(g, v, world, nQubits, target);
        v = gc(v);
    }
    angle = get_random();
    for (int target = 0; target < nQubits; target++) {
        auto g = RX(nQubits, target, angle);
        v = mv_multiply_MPI(g, v, world, nQubits, target);
        v = gc(v);
    }
    return v;
}

int main(int argc, char** argv) {
    bmpi::environment env(argc, argv);
    bmpi::communicator world;

    auto start = std::chrono::high_resolution_clock::now();

    auto v = exec(std::atoi(argv[1]), world);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> ms = end - start;
    if (world.rank() == 0)
        std::cout << "MPIsize " << world.size() << " nQubits "
                  << std::atoi(argv[1]) << " nNodes " << get_nNodes(v) << " "
                  << ms.count() / 1000000 << " sec" << std::endl;
}