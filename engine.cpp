#include "engine.h"


NodeTable uniqueTable(8*NodeTable::region_size, 16*NodeTable::region_size);

void JobQueue::push(Job* job){
    long long bottom = _bottom.load(std::memory_order_acquire);
    _jobs[bottom&MASK] = job;
    _bottom.store(bottom+1, std::memory_order_release); 
}


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


Job* JobQueue::steal(){
    long long top = _top.load(std::memory_order_seq_cst);
    long long bottom = _bottom.load(std::memory_order_seq_cst);

    if(top < bottom){
        Job* job = _jobs[top&MASK];
        
        long long expectedTop = top;
        long long desiredTop= top + 1LL;
        if(!_top.compare_exchange_weak(expectedTop, desiredTop, std::memory_order_acq_rel)){
            job = nullptr;
        }
        return job;
    }else{
        return nullptr; 
    }
}


std::size_t JobQueue::size() const {
    long long bottom = _bottom.load();
    long long top = _top.load();
    return bottom - top;
}

bool JobQueue::empty() const {
    return size() == 0;
}


Complex Worker::getComplexFromCache(const double_pair& p){
    return ccache.getCached_v(p);
}

void Worker::returnComplexToCache(Complex && c){
    ccache.returnCached(std::move(c));
}

Complex Worker::getComplexFromTable(const double_pair& p){
    return ctable.find_or_insert(p);
}

Index Worker::uniquefy(const mNode& n){
    assert(std::all_of(n.children.begin(), n.children.end(), [](const mEdge& e){
        return e.w.src == Complex::Table;
    }));
   return uniqueTable.find_or_insert(n); 
}
