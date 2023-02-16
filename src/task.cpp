
#if defined(__linux__)
#include <sched.h>
#endif


#include "task.h"
#include <iostream>
#include <chrono>



using namespace std::chrono_literals;





void Scheduler::spawn(){





    for(int i = 0; i < _nworkers; i++){
        _workers[i]._id = i;
        _workers[i]._sched = this;
 

        _workers[i]._thread = new std::thread([this](int id){

            boost::fibers::use_scheduling_algorithm<boost::fibers::algo::work_stealing>(this->_nworkers + 1);
            {
                std::unique_lock<boost::fibers::mutex> lk(this->mtx_stop);
                this->cond_stop.wait(lk);
            }

        }, i);
        
                
    }

#if defined(__linux__)
    for(int i = 0; i < _nworkers; i++){
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        if(pthread_setaffinity_np(_workers[i]._thread->native_handle(), sizeof(cpu_set_t), &cpuset)){
            std::cout<<"pthread_setaffinity_np failed"<<std::endl;
            exit(1);
        }
    }
#endif

    boost::fibers::use_scheduling_algorithm<boost::fibers::algo::work_stealing>(this->_nworkers + 1);
}

Scheduler::Scheduler(int n ): _nworkers(n){

    for(int i = 0; i < _nworkers; i++){
        _workers.emplace_back(WorkerThread());
    }
    this->spawn();

}

Scheduler::~Scheduler(){
    cond_stop.notify_all();
    
    for(WorkerThread& w: _workers){
        w._thread->join();
    }

}

void Scheduler::addGate(const mEdge& e){
    _gates.emplace_back(e);
}


vEdge Scheduler::buildCircuit(vEdge input){

    vEdge v = input;

    for(auto i = 0; i < _gates.size(); i++){
        v = mv_multiply_fiber(_gates[i], v);
    }



    return v;

}
