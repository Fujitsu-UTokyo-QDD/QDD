#include "engine.h"
#include <cassert>
#include "table.hpp"





Index Worker::uniquefy(const mNode& n){

   //return m_uniqueTable.find_or_insert(this, n); 
}





void Worker::run() {
    
    _thread = std::thread([this](){
            while(!(*_stop)){
                if(auto job = _eng->request()){
                    auto t1 = std::chrono::high_resolution_clock::now(); 
                    job.value()->execute(this); 
                    auto t2 = std::chrono::high_resolution_clock::now(); 
                    auto d = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
                    timer += d;
                    executed++;
                }
                
            }
    
            return;

     });

}

/*
void Worker::run() {
    
    _thread = std::thread([this](){
            while(!(*_stop)){
                Job* j;
                if(_queue.pop(j)){
                    auto t1 = std::chrono::high_resolution_clock::now(); 
                    j->execute(this);
                    auto t2 = std::chrono::high_resolution_clock::now(); 
                    auto d = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
                    timer += d;
                    //if(_id == 2) std::cout<<"task"<<executed<<":  "<< d.count()<<std::endl; 
                    executed++;
                }
                
            }
    
            return;

     });

}
*/
Job* Engine::steal() {
    return nullptr;
}

std::size_t Engine::next_worker() {
    auto w = _current_worker.load();
    _current_worker.store((_current_worker+1)%_total_worker);
    return w;
}



