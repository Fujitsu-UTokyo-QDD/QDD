#include "common.h"
#include "table.hpp"
#include <iostream>
#include "engine.h"
#include "dd.h"
#include <random>
#include <chrono>


using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

auto benchmark(Engine* eng){
    GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};

    std::vector<mEdge> gate_queue;
    const std::size_t NGATES = 50;
    const uint64_t NQUBITS = 8;

    std::mt19937_64 rng;
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);

    std::uniform_int_distribution<int> qubit_dist(0, NQUBITS-1);
    std::uniform_int_distribution<int> gate_dist(0, gates->size() - 1 );

    std::vector<Job*> jobs;
    Controls  c = {};
    auto t1 = std::chrono::high_resolution_clock::now();
    for(std::size_t i = 0; i < NGATES; i++){
        auto g = gates[gate_dist(rng)];
        Job* j = eng->submit(makeGate, g, NQUBITS, qubit_dist(rng),Controls{});
        jobs.push_back(j);
    }
    for(auto i = 0; i < NGATES; i++){
        assert(jobs[i]->valid());
        jobs[i]->getResult();
    }
    
    //jobs[10]->getResult().printMatrix(); 
    std::cout<<"all results should be available"<<std::endl; 
    auto result = eng->addReduce(jobs, NGATES/eng->worker_number());
    auto t2 = std::chrono::high_resolution_clock::now();
    
    duration<double, std::milli> ms = t2 - t1;
    std::cout<<ms.count()<<" ms"<<std::endl;
    return;
}

int main()
{
    Engine eng(8, 20);
    benchmark(&eng);
    eng.terminate();
}

