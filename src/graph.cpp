#include <sched.h>
#include "graph.hpp"
#include <mutex>
#include <exception>
#include <condition_variable>
#include <optional>
#include <pthread.h>

oneapi::tbb::enumerable_thread_specific<std::chrono::duration<double, std::micro>> steal_timer{std::chrono::duration<double, std::micro>::zero()};


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
            auto t1 = std::chrono::high_resolution_clock::now();

            _workers[id]->execute();
            if(id == 0) _workers[id]->collect(_nworkers);
            auto t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::micro> ms = t2 - t1;
            std::cout<<id<<" finished in "<<ms.count()<<" micro s"<<std::endl;


        }, i, std::ref(mtx), std::ref(cond), std::ref(spawned));
        
                
    }

#if defined(__linux__)
    for(int i = 0; i < _nworkers; i++){
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        if(pthread_setaffinity_np(_workers[i]->_thread->native_handle(), sizeof(cpu_set_t), &cpuset)){
            std::cout<<"pthread_setaffinity_np failed"<<std::endl;
            exit(1);
        }
    }
#endif

    std::unique_lock<std::mutex> lock{mtx};
    cond.wait(lock, [&]() { return spawned == _nworkers; });
}



void Executor::try_execute_self(Worker* w){
    
    while(auto job = w->_wsq.pop()){
        job.value()->execute(w);
    }
   
}


void Executor::try_execute_else(Worker* w){
    /*
    if(auto job = this->_wsq.steal()){
        job.value()->execute(w);
    }
    else{
        int victim = w->_dist(w->_rdgen);
        if(auto job = this->_workers[victim]->_wsq.steal()){
            job.value()->execute(w);
        }
        
    }
    */
    std::size_t num_steals = 0;

    do{
        int victim = w->_dist(w->_rdgen);
        if(auto job = this->_workers[victim]->_wsq.steal()){
            job.value()->execute(w);
        }else{
            num_steals++;
            if(num_steals == MAX_STEALS){
                std::this_thread::yield();
                break;
            }
        }
    }while(true);


}


void Worker::execute(){

    while(this_round.size() > 1){
        for(Node* n: this_round){
            n->execute(this);
        }
        this_round.clear();
        this_round.swap(next_round);
    }

    assert(this_round.size() == 1);
    _executor->_total_queue[_id] = this_round[0];
    
   this->_executor->_sem.release(); 

}


void Worker::collect(int N){
    do{
        _executor->_sem.acquire();
    }while(--N > 0);

    this_round.clear();
    next_round.clear();

    this_round = std::move(_executor->_total_queue);


    while(!this_round.empty()){
        for(Node* n: this_round){
            n->execute(this);
        }
        this_round.clear();
        this_round.swap(next_round);
    }

}
