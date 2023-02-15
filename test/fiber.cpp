#include <iostream>
#include "task.h"
#include <functional>
#include <chrono>
#include <boost/fiber/future/future.hpp>
#include <boost/fiber/future/packaged_task.hpp>
using namespace std::chrono_literals;



int bar(int i){
    for(int j = 0; j < 10000; j++){}
    //boost::this_fiber::sleep_for(std::chrono::seconds(i));
    return i;
}

void foo(int until) {
    std::vector<boost::fibers::future<int>> results;
    for(int i = 0; i < until; i++){
        boost::fibers::packaged_task<int()> pt(std::bind(bar, i)); 
        results.emplace_back(pt.get_future());
        boost::fibers::fiber(std::move(pt)).detach();
    }

    for(int i = 0; i < until; i++){
        std::cout<<i<<": "<<results[i].get()<<std::endl;
    }

}

int main(int argc, char* argv[]){

    foo(4);
    return 0;
    Scheduler s(4);
    s.spawn();

    s.addGate(RZ(15,0,1.7952706710012407));
    s.addGate(RY(15, 0, 1.0056905557557458));
    s.addGate(RZ(15,0,-2.860782987649066));
    s.addGate(RY(15,1,2.9482444835445483));
    s.buildCircuit();
    return 0;
}
