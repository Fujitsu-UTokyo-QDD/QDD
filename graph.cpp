#include "graph.hpp"
#include <mutex>
#include <exception>
#include <optional>



void Executor::spawn(){

    std::mutex mtx;
    std::condition_variable cond;
    int spawned = 0;

    for(int i = 0; i < _nworkers; i++){
        _workers[i]._id = i;
        _workers[i]._executor = this;
 

        _workers[i]._thread = new std::thread([this](int id, std::mutex& mtx, std::condition_variable& cond, int& spawned){
            
            {
                std::unique_lock<std::mutex> lock(mtx);    
                if(++spawned == _nworkers){
                    cond.notify_one();
                }
            }

            std::exception_ptr ptr{nullptr};
            try {
                while(1) {
                    try_execute_self(&_workers[id]);
                    if(!try_execute_else(&_workers[id])) break;
                }
            } catch(...) {
                ptr = std::current_exception();
            }

        }, i, std::ref(mtx), std::ref(cond), std::ref(spawned));
        
                
    }

    std::unique_lock<std::mutex> lock{mtx};
    cond.wait(lock, [&]() { return spawned == _nworkers; });
    std::cout<<"finish"<<std::endl;
}



void Executor::try_execute_self(Worker* w){
    
    while(auto job = w->_wsq.pop()){
           
    }
   

}


bool Executor::try_execute_else(Worker* w){

}


void Executor::invoke(Worker* w, Node* n){

}

