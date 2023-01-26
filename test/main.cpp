#include <iostream>
#include "graph.hpp"
#include "dd.h"
#include <chrono>
#include "algorithms/grover.hpp"
#include "algorithms/shor.hpp"
#include <cstdlib>



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




int main(int argc, char* argv[]){

/*     
    int q = std::atoi(argv[1]);
    int workers = std::atoi(argv[2]);


    for(auto r = 20; r<= 120 ; r += 20){
        std::cout<<r<<" reduce"<<std::endl; 
        Grover g(q, workers, r) ;   
        g.full_grover();
    }
    
*/
    int n = std::atoi(argv[1]);
    int w = std::atoi(argv[2]);
    int r = std::atoi(argv[3]);
    Shor s(n, w, r, true);    
    s.run();


    return 0;
}
