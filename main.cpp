#include "common.h"
#include "table.hpp"
#include <iostream>
#include "engine.h"
#include "dd.h"
#include <random>
#include <chrono>
#include "lockfree_hashmap.hpp"
#include <bitset>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb.h>



using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


using namespace oneapi::tbb;


auto benchmark_dense(Engine* eng){
    GateMatrix gates[] = {Hmat, SXmat, SXdagmat, Vmat, Vdagmat};

    std::vector<mEdge> gate_queue;
    const std::size_t NGATES = 50;
    const uint64_t NQUBITS = 10;

    std::mt19937_64 rng;
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);

    std::uniform_int_distribution<int> qubit_dist(0, NQUBITS-1);
    std::uniform_int_distribution<int> gate_dist(0, gates->size() - 1 );

    std::vector<Job*> jobs;
    Controls  c = {};
    for(std::size_t i = 0; i < NGATES; i++){
        std::vector<mEdge> single_gate;
        for(auto j = 0; j < NQUBITS; j++ ){
            auto g = gates[gate_dist(rng)];
            Job* job = eng->submit(makeGate, g, 1,0, Controls{});
            single_gate.push_back(job->getResult());
        }
        mEdge e = single_gate[0];
        for(auto j = 1; j < NQUBITS -1; j++) e = eng->submit(kronecker, e, single_gate[j])->getResult();
        Job* j = eng->submit(kronecker, e, single_gate[NQUBITS-1]);
        jobs.push_back(j);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto result = eng->mulReduce(jobs, 5);
    auto t2 = std::chrono::high_resolution_clock::now();
    
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" macro sec"<<std::endl;
    return;
}

auto benchmark(Engine* eng){
    GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};

    std::vector<mEdge> gate_queue;
    const std::size_t NGATES = 50;
    const uint64_t NQUBITS = 3;

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
    /*
    for(auto j: jobs){
        j->getResult();
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout<<"makeGate: "<< duration<double, std::micro>(t2-t1).count()<<std::endl;
    */
    auto result = eng->mulReduce(jobs, 5);
    auto t2 = std::chrono::high_resolution_clock::now();
    
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" macro sec"<<std::endl;
    return;
}


struct Bar {
    Qubit v;
    double b;
    Bar* next;
    Bar* self;
    inline bool operator==(const Bar& other) const noexcept {
        return b == other.b && v == other.v;
    }
};

template<>
struct std::hash<Bar>{
    std::size_t operator()(const Bar& b) const noexcept {
        return std::hash<double>()(b.b);
    }
};


void test_ht(){
    CHashTable<Bar> ht(10); 

    std::mt19937_64 rng;
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};

    rng.seed(ss);

    std::uniform_real_distribution<> val_dist(0.0, 5.0);
    std::uniform_int_distribution<> qubit_dist(0, 9);


    std::vector<int> qubits;
    std::vector<double> vals;
    
    std::size_t TESTS = 10000;
    for(auto i = 0 ; i< TESTS; i++){
        qubits.emplace_back(qubit_dist(rng));
        vals.emplace_back(val_dist(rng));
    }

    assert(qubits.size() == TESTS && vals.size() == TESTS);
   parallel_for(blocked_range<std::size_t>(0, TESTS), [&](blocked_range<std::size_t> r){
        for(auto i = r.begin(); i != r.end(); i++){
           Bar* b = ht.getNode();
           b->v = qubits[i];
           b->b = vals[i];
           Bar* ret = ht.lookup(b);
           ret->self = ret;
        }
    }); 


   parallel_for(blocked_range<std::size_t>(0, TESTS), [&](blocked_range<std::size_t> r){
        for(auto i = r.begin(); i != r.end(); i++){
           Bar* b = ht.getNode();
           b->v = qubits[i];
           b->b = vals[i];
           Bar* ret = ht.lookup(b);
           assert(ret->self == ret);
        }
    }); 

    std::cout<<"ht success"<<std::endl;

}

void test_basic(){

    //GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};
    GateMatrix gates[] = {Hmat, SXmat, SXdagmat, Vmat, Vdagmat};

    std::vector<mEdge> gate_queue;
    const std::size_t NGATES = 100;
    const uint64_t NQUBITS = 10;

    std::mt19937_64 rng;
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(0);

    std::uniform_int_distribution<int> qubit_dist(0, NQUBITS-1);
    std::uniform_int_distribution<int> gate_dist(0, gates->size() - 1 );

    Controls  c = {};
    for(std::size_t i = 0; i < NGATES; i++){
        std::cout<<i<<std::endl;
        auto g = gates[gate_dist(rng)];
        auto e = makeGate(nullptr, g, NQUBITS, qubit_dist(rng), c);
        //e.printMatrix();
        gate_queue.push_back(e);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto result = multiply(nullptr, gate_queue[0], gate_queue[1]);
    for(auto i = 2; i < NGATES; i++){
        result = multiply(nullptr, result, gate_queue[i]);
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" micro sec"<<std::endl;
    return;


}


int main(){

    test_basic();

 
    return 0;
}


