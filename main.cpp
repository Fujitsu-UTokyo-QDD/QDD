#include <iostream>
#include "dd.h"
#include <random>
#include <chrono>
#include <bitset>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb.h>
#include <future>
#include <cstdlib>
#include "graph.hpp"



using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


using namespace oneapi::tbb;




int main(int argc, char* argv[]){
    /*
    QuantumCircuit qc(2);
    qc.emplace_back(Hmat, 1);
    qc.emplace_back(Xmat, 0);
    qc.emplace_back(Xmat, 0);
    qc.emplace_back(Xmat, 0);
    qc.emplace_back(Xmat, 0);
    qc.emplace_back(Xmat, 0);
    qc.emplace_back(Xmat, 0);
    qc.buildCircuit();
    qc.dump_task_graph();
    */

    Executor e(10);
    e.spawn();
}


