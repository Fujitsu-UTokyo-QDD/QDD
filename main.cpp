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
#include <future>
#include <cstdlib>



using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


using namespace oneapi::tbb;




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
        auto g = gates[gate_dist(rng)];
        auto e = makeGate(nullptr, g, NQUBITS, qubit_dist(rng), c);
        //e.printMatrix();
        gate_queue.push_back(e);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto r1 = std::async(std::launch::async, [&](){
        auto result = multiply(nullptr, gate_queue[0], gate_queue[1]);
        for(auto i = 2; i < NGATES/2; i++){
            result = multiply(nullptr, result, gate_queue[i]);
        }
        return result;
        });

    auto r2 = std::async(std::launch::async, [&](){
        auto result = multiply(nullptr, gate_queue[NGATES/2], gate_queue[NGATES/2+1]);
        for(auto i = NGATES/2+2; i < NGATES; i++){
            result = multiply(nullptr, result, gate_queue[i]);
        }
        return result;
            });

    auto result = multiply(nullptr, r1.get(), r2.get());
    auto t2 = std::chrono::high_resolution_clock::now();

    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" micro sec"<<std::endl;
    return;


}


void test_kronecker(){

    //GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};
    GateMatrix gates[] = {Hmat, SXmat, SXdagmat, Vmat, Vdagmat};

    mEdge e1 = makeGate(nullptr, gates[1], 1, 0, Controls{});
    mEdge e2 = makeGate(nullptr, gates[2], 1, 0, Controls{});
    mEdge e3 = makeGate(nullptr, gates[3], 1, 0, Controls{});
    mEdge result = kronecker(nullptr, e1, e2);
    result = kronecker(nullptr, result, e3);
    result.printMatrix();


    return;


}

static mEdge make_dense(QubitCount q, int seed = 0){
    
    GateMatrix gates[] = {Hmat, SXmat, SXdagmat, Vmat, Vdagmat};
    std::mt19937_64 rng;
    rng.seed(seed);

    std::uniform_int_distribution<int> gate_dist(0, gates->size() - 1 );
    
    std::vector<mEdge> gate_queue;
    for(auto i = 0; i < q; i++){
       gate_queue.push_back(makeGate(nullptr, gates[gate_dist(rng)], 1, 0, Controls{})) ;
    }
    mEdge result = kronecker(nullptr, gate_queue[0], gate_queue[1]);

    for(auto i = 2; i < q; i++){
        result = kronecker(nullptr, result, gate_queue[i]);
    }

    return result;

}

auto benchmark_dense(){

    std::vector<mEdge> gate_queue;
    const std::size_t NGATES = 100;
    const uint64_t NQUBITS = 10;


    Controls  c = {};
    for(auto i = 0; i < NGATES; i++){
        gate_queue.emplace_back(make_dense(NQUBITS, i));
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto r1 = std::async(std::launch::async, [&](){
        auto result = multiply(nullptr, gate_queue[0], gate_queue[1]);
        for(auto i = 2; i < NGATES/2; i++){
            result = multiply(nullptr, result, gate_queue[i]);
        }
        return result;
        });

    auto result = multiply(nullptr, gate_queue[NGATES/2], gate_queue[NGATES/2+1]);
    for(auto i = NGATES/2+2; i < NGATES; i++){
        result = multiply(nullptr, result, gate_queue[i]);
    }

    result = multiply(nullptr, result, r1.get());
    auto t2 = std::chrono::high_resolution_clock::now();
    
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" macro sec"<<std::endl;
    return;
}

void test_vec(){
    auto v1 = makeZeroState(nullptr, 4);
    v1.printVector();
    std::cout << std::endl;
    auto v2 = makeOneState(nullptr, 2);
    v2.printVector();
    return;
}

struct F {
    int a;
    F(){ std::cout<<"default"<<std::endl;}
    F(int aa):a(aa) {std::cout<<"user"<<std::endl;}
    F(const F& f):a(f.a) {std::cout<<"copy"<<std::endl;}
    F(F&& f):a(f.a) {std::cout<<"move"<<std::endl;}
};

template<typename F, typename... Args>
void foo(F&& f, Args&&... args){
    std::cout<<"enter foo"<<std::endl; 
    auto t = std::bind(f, std::forward<Args>(args)...);
    t();
}

int main(int argc, char* argv[]){
    int nthreads;
    if(argc > 1) nthreads = std::atoi(argv[1]);
    else nthreads = 4;
    std::cout<<"Use "<<nthreads<<" threads"<<std::endl;
    Engine eng(nthreads, 20);
    mEdge e1 = make_dense(10, 0);
    mEdge e2 = make_dense(10, 1);
    Job* j = eng.submit(add, e1, e2);
    auto t1 = std::chrono::high_resolution_clock::now();
    mEdge result = j->getResult();
    auto t2 = std::chrono::high_resolution_clock::now();
    
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" macro sec"<<std::endl;
    eng.terminate();
}


