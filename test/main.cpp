#include <iostream>
#include "graph.hpp"
#include "dd.h"
#include <chrono>
#include "algorithms/grover.hpp"
#include "algorithms/shor.hpp"
#include <cstdlib>
#include <cmath>



using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

vEdge run(int r){
    vEdge input = makeZeroState(10);
    QuantumCircuit qc(10,8, r,input);
    qc.emplace_back(Hmat, 1);
    qc.emplace_back(Xmat, 2);
    qc.emplace_back(Ymat, 3);
    qc.emplace_back(Vdagmat, 4);
    qc.emplace_back(Hmat, 5);
    qc.emplace_back(Sdagmat, 6);
    qc.emplace_back(Hmat, 7);
    qc.emplace_back(Xmat, 8);
    qc.emplace_back(Ymat, 9);
    qc.emplace_back(Vdagmat, 1);
    qc.emplace_back(Hmat, 2);
    qc.emplace_back(Sdagmat, 3);
    qc.emplace_back(Hmat, 1);
    qc.emplace_back(Xmat, 2);
    qc.emplace_back(Ymat, 3);
    qc.emplace_back(Vdagmat, 4);
    qc.emplace_back(Hmat, 5);
    qc.emplace_back(Sdagmat, 6);
    qc.emplace_back(Hmat, 7);
    qc.emplace_back(Xmat, 8);
    qc.emplace_back(Ymat, 9);
    qc.emplace_back(Vdagmat, 1);
    qc.emplace_back(Hmat, 2);
    qc.emplace_back(Sdagmat, 3);
    qc.buildCircuit();
    std::cout<<"Task graph: "<<std::endl;
    qc.dump_task_graph();
    auto t1 = std::chrono::high_resolution_clock::now();
    vEdge result = qc.wait().vectorResult();
    auto t2 = std::chrono::high_resolution_clock::now();
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" micros"<<std::endl;
    return result;

}



//[start, end)
static void qft_rotations(QuantumCircuit& qc, Qubit start, Qubit end){
    if(start == end) 
        return;

    end -= 1;
    qc.emplace_back(Hmat, end);

    for(Qubit q = start; q < end; q++){
        Controls c{Control{q, Control::Type::pos}};
        float r = std::cos(std::numbers::pi/(std::pow(2, end - q))); 
        float i = std::sin(std::numbers::pi/(std::pow(2, end - q))); 
        GateMatrix g{cf_one, cf_zero, cf_zero, {r,i}}; 
        qc.emplace_back(g, end, c);
    }

    qft_rotations(qc, start, end);
    
}


//[start, end)
static void qft_swap(QuantumCircuit& qc, Qubit start, Qubit end){
    QubitCount total_qubits = qc.getQubits();
    Qubit e = end-1;
    for(Qubit q = start; q < (start+end)/2; q++){
        qc.emplace_gate(makeSwap(total_qubits, q, e-- ));
    }
}

static void full_qft(QuantumCircuit& qc, Qubit start, Qubit end){
    qft_rotations(qc, start, end);
    qft_swap(qc, start, end);
}

int main(int argc, char* argv[]){


    int q = std::atoi(argv[1]);
    int workers = std::atoi(argv[2]);
    int r = std::atoi(argv[3]);

        Grover g(q, workers, r) ;   
        g.full_grover();
        return 0;
    /*
    

    int n = std::atoi(argv[1]);
    int w = std::atoi(argv[2]);
    int r = std::atoi(argv[3]);
    Shor s(n, w, r, true);    
    s.run();

    */
    return 0;
    
}
