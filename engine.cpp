#include "engine.h"
#include <cassert>


NodeTable uniqueTable(32*NodeTable::region_size, 64*NodeTable::region_size);

/*
void JobQueue::push(Job* job){
    long long bottom = _bottom.load(std::memory_order_acquire);
    _jobs[bottom&MASK] = job;
    _bottom.store(bottom+1, std::memory_order_release); 
}
*/
void JobQueue::push(Job* job){
    while(true){
        long bot = _bottom;
        long expected = bot + 1;
        Job* j = _jobs[bot];
        _jobs[bot] = job;
        if(__sync_bool_compare_and_swap(&_bottom, bot, expected)){
            return;
        }else{
            _jobs[bot] = j;
            continue;
        }
    }
}
/*
Job* JobQueue::pop(){
    
    long long bottom = _bottom.load(std::memory_order_acquire);
    bottom = std::max(0LL, bottom - 1);
    _bottom.store(bottom, std::memory_order_release);

    long long top = _top.load();

    if(top <= bottom){
    
        Job* job = _jobs[bottom&MASK];
        if(top != bottom)
        {
            return job;
            
        }else{
            //there is only one job(or none) left
            //so we have to compete with other steals' action on top
            long long expectedTop = top;
            long long desiredTop= top + 1LL;
            if(!_top.compare_exchange_weak(expectedTop, desiredTop, std::memory_order_acq_rel)){
                job = nullptr;
            }

            _bottom.store(top + 1LL, std::memory_order_release);
            return job;
        }
    }else{
        _bottom.store(top, std::memory_order_release);
        return nullptr;
    }

}
*/
Job* JobQueue::pop(){
   
    long bot = _bottom;
    long expected = bot - 1;
    if(bot == 0) return nullptr;
    Job* j = _jobs[expected];

    if(__sync_bool_compare_and_swap(&_bottom, bot, expected)){
        return j;
    }else{
        return nullptr;
    }

}

Job* JobQueue::steal(){

    return nullptr; 

}


std::size_t JobQueue::size() const {
    long long bottom = _bottom;
    long long top = _top;
    return bottom - top;
}

bool JobQueue::empty() const {
    return size() == 0;
}


Complex Worker::getComplexFromCache(const double_pair& p){
    return ccache.getCached_v(p);
}

void Worker::returnComplexToCache(Complex & c){
    ccache.returnCached(c);
}

Complex Worker::getComplexFromTable(const double_pair& p){
    double r = p.first * p.first + p.second * p.second;
    r = std::sqrt(r);
    assert(r <= 1);
    return ctable.find_or_insert(p);
}

Index Worker::uniquefy(const mNode& n){
    assert(std::all_of(n.children.begin(), n.children.end(), [](const mEdge& e){
        return e.w.src == Complex::Table;
    }));
   return uniqueTable.find_or_insert(n); 
}


void Worker::submit(Job* j){
    while(!_queue.push(j));
}


void Worker::run() {
    
    _thread = std::thread([this](){

            while(!(*_stop)){
                Job* j;
                if(_queue.pop(j)){
                    j->execute(this);
                }else{
                    //need to steal
                    j = this->_eng->steal();
                    if(j != nullptr){
                        j->execute(this);
                    }else {
                        std::this_thread::yield();
                    }
                }
                
            }
    
            if(_queue.size() != 0) assert(false);
            return;

     });

}

Job* Engine::steal() {
    return nullptr;
}

std::size_t Engine::next_worker() {
    auto w = _current_worker.load();
    _current_worker.store((_current_worker+1)%_total_worker);
    return w;
}



