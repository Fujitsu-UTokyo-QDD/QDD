#include "engine.h"
#include <cassert>
#include <thread>
#include "table.hpp"






void Worker::run() {
    
    _thread = std::thread([this](){
            while(!(*_stop) || _local_jobs.size() > 0){
                if(auto job = _local_jobs.pop()){
                    //auto t1 = std::chrono::high_resolution_clock::now(); 
                    job.value()->execute(this); 
                    //auto t2 = std::chrono::high_resolution_clock::now(); 
                    //auto d = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
                    //timer += d;
                    executed++;
                }else if(auto job = _eng->steal(worker_dist(rng))){
                    job.value()->execute(this); 
                    executed++;
                }else{
                    std::this_thread::yield();
                }
                
            }
    
            return;

     });

}



std::optional<Job*> Engine::steal(std::size_t target) {
    if(auto j = _local_jobs.steal()){
        return j;
    }

    if(_workers[target]->_local_jobs.size() > 0)
        return _workers[target]->_local_jobs.steal();
    else
        return std::nullopt;

    
}





