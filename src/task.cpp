
#if defined(__linux__)
#include <sched.h>
#endif


#include "task.h"
#include <iostream>
#include <chrono>



using namespace std::chrono_literals;



Scheduler::Scheduler(int n ): _stop(false), _ready(false), _nworkers(n){

    for(int i = 0; i < _nworkers; i++){
        _workers.emplace_back(WorkerThread());
    }

}

Scheduler::~Scheduler(){
    
    _stop = true;
    for(WorkerThread& w: _workers){
        w._thread->join();
    }

}


void Scheduler::spawn(){




    std::mutex mtx;
    std::condition_variable cond;
    int spawned = 0;

    for(int i = 0; i < _nworkers; i++){
        _workers[i]._id = i;
        _workers[i]._sched = this;
 

        _workers[i]._thread = new std::thread([this](int id, std::mutex& mtx, std::condition_variable& cond, int& spawned){
            
            {
                std::unique_lock<std::mutex> lock(mtx);    
                if(++spawned == _nworkers){
                    cond.notify_one();
                }
            }


            boost::fibers::use_scheduling_algorithm<boost::fibers::algo::work_stealing>(this->_nworkers + 1);

            std::exception_ptr ptr{nullptr};
            try {
                while(!this->_stop.load()) {
                    boost::this_fiber::sleep_for(500ms);
                    std::cout<< id <<std::endl;
                }
            } catch(...) {
                ptr = std::current_exception();
            }

        }, i, std::ref(mtx), std::ref(cond), std::ref(spawned));
        
                
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

    std::unique_lock<std::mutex> lock{mtx};
    cond.wait(lock, [&]() { return spawned == _nworkers; });
    boost::fibers::use_scheduling_algorithm<boost::fibers::algo::work_stealing>(this->_nworkers + 1);
}


void Scheduler::addGate(const mEdge& e){
    _gates.emplace_back(e);
}


mEdge Scheduler::buildCircuit(){

    mEdge lhs, rhs;
    if(_gates.size() == 0){
        return mEdge();
    }else if(_gates.size() == 1){
        return _gates[0]; 
    }

    lhs = _gates[0];
    std::size_t next = 1;
    do{
        lhs = mm_multiply_fiber(lhs, _gates[next]);
        next++;
    
    }while(next < _gates.size());
    std::this_thread::sleep_for(3s);
    return lhs;

}
