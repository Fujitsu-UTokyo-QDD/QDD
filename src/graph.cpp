#include <sched.h>
#include "graph.hpp"
#include <mutex>
#include <exception>
#include <condition_variable>
#include <optional>
#include <pthread.h>



void Executor::spawn(){

    std::mutex mtx;
    std::condition_variable cond;
    int spawned = 0;

    for(int i = 0; i < _nworkers; i++){
        _workers[i]->_id = i;
        _workers[i]->_executor = this;
        _workers[i]->_dist = std::uniform_int_distribution<int>(0, this->_nworkers - 1);
 

        _workers[i]->_thread = new std::thread([this](int id, std::mutex& mtx, std::condition_variable& cond, int& spawned){
            
            {
                std::unique_lock<std::mutex> lock(mtx);    
                if(++spawned == _nworkers){
                    cond.notify_one();
                }
            }

            std::exception_ptr ptr{nullptr};
            try {
                while(!(*this->_stop)) {
                    try_execute_self(_workers[id]);
                    try_execute_else(_workers[id]);
                }
            } catch(...) {
                ptr = std::current_exception();
            }

        }, i, std::ref(mtx), std::ref(cond), std::ref(spawned));
        
                
    }

    for(int i = 0; i < _nworkers; i++){
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        if(pthread_setaffinity_np(_workers[i]->_thread->native_handle(), sizeof(cpu_set_t), &cpuset)){
            std::cout<<"pthread_setaffinity_np failed"<<std::endl;
            exit(1);
        }
    }

    std::unique_lock<std::mutex> lock{mtx};
    cond.wait(lock, [&]() { return spawned == _nworkers; });
}



void Executor::try_execute_self(Worker* w){
    
    while(auto job = w->_wsq.pop()){
        job.value()->execute(w);
    }
   
}


void Executor::try_execute_else(Worker* w){

    if(auto job = this->_wsq.steal()){
        job.value()->execute(w);
    }
    else{
        int victim = w->_dist(w->_rdgen);
        if(auto job = this->_workers[victim]->_wsq.steal()){
            job.value()->execute(w);
        }
        
    }


}



