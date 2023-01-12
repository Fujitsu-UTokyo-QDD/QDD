#include "common.h"
#include "table.hpp"
#include <iostream>
#include "engine.h"
#include "dd.h"
#include <random>
#include <chrono>
#include "lockfree_hashmap.hpp"
#include <bitset>


using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

auto benchmark(Engine* eng){
    GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};

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
    auto t1 = std::chrono::high_resolution_clock::now();
    for(std::size_t i = 0; i < NGATES; i++){
        auto g = gates[gate_dist(rng)];
        Job* j = eng->submit(makeGate, g, NQUBITS, qubit_dist(rng),Controls{});
        jobs.push_back(j);
    }

    auto result = eng->mulReduce(jobs, NGATES/eng->worker_number());
    auto t2 = std::chrono::high_resolution_clock::now();
    
    duration<double, std::milli> ms = t2 - t1;
    std::cout<<ms.count()<<" ms"<<std::endl;
    return;
}
struct Bar {
    double i;

    inline bool operator==(const Bar& other) const noexcept {
        return i == other.i;
    }
};

int test()
{


    
    LockFreeMap<Bar> map;




    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0,1);

    constexpr std::size_t num = 100000;
    double vals[num];
    Bar* ans[num];
    for(auto i = 0; i < num; i++){
        vals[i]=dist(e2);
    }

     
    std::vector<std::thread> threads1;
    std::vector<std::thread> threads2;
    int nthread = 10;
    std::size_t step = num/nthread;
    for(int i = 0; i < nthread; i++){
        threads1.emplace_back([ &](int ii){
            for(auto t = ii * step; t < (ii+1)*step; t++){  
                Bar* b = map.find_or_insert(double_to_bits(vals[t]), { .i = vals[t]});
                ans[t] = b;
            }
        }, i);
    }

    for(int i = 0; i < nthread; i++){
        threads1[i].join();
    }
    for(int i = 0; i < nthread; i++){
        threads2.emplace_back([&](int ii){
            for(auto t = ii * step; t < (ii+1)*step; t++){  
                Bar* b = map.find_or_insert(double_to_bits(vals[t]), { .i = vals[t]});
                assert(b == ans[t]);
                assert(b->i == vals[t]);
            }
        }, i);
    }
    for(int i = 0; i < nthread; i++){
        threads2[i].join();
    }

    std::cout<<"success"<<std::endl;

    return 0;    
}
int main(){
    Engine eng(8,20);
    benchmark(&eng);
    eng.terminate();
    return 0;
}


