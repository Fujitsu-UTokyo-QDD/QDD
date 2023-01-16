#pragma once
#include <vector>
#include <type_traits>
#include "dd.h"
#include <future>
#include <functional>
#include <thread>
#include <unordered_map>
#include <atomic>
#include "queue.hpp"
#include <chrono>
#include <iostream>
#include "lockfree_hashmap.hpp"
#include "chase_lev_queue.hpp"



/**
 *  lhs, rhs, op, result
 *
 * */
class Worker;

class Job{
    public:

        template<typename F, typename... Args>
        Job(std::size_t p, F&& f, Args&&... args){
            using namespace std::placeholders;
            using return_type = std::invoke_result_t<F, Worker*, Args...>;
            static_assert(std::is_same_v<return_type, mEdge>);

            auto pt = std::packaged_task<mEdge(Worker*)>(std::bind(std::forward<F>(f), _1, std::forward<Args>(args)...));
            task = std::move(pt);
            result = task.get_future();
            
            _parent = p;
        }

        void execute( Worker* w) {
           task(w);
        }

        bool valid(){
            return result.valid();
        }

        mEdge getResult() {
            return result.get();
        }
    private:
        std::packaged_task<mEdge(Worker*)> task;
        std::shared_future<mEdge> result;

        std::size_t _parent;

};



class Engine;

class Worker{
    friend class Engine;
    public:

        Worker(Engine* eng, std::size_t id,  bool* stop): _eng(eng), _id(id),  _stop(stop), _queue(1024), timer(std::chrono::microseconds::zero()), _mregion(-1), _vregion(-1){};
    //Worker(Engine* eng, std::size_t id,  bool* stop): _eng(eng), _id(id),  _stop(stop), timer(std::chrono::microseconds::zero()), _mregion(-1), _vregion(-1){};

        void run();


        template<typename T>
        int64_t& get_region();


        template<typename F, typename... Args>
            Job* submit(F&& f, Args&&... args){

            }

        Index uniquefy(const mNode& n);

        Engine* _eng;
        bool* _stop;
        std::thread _thread;
        std::size_t _id;

        int64_t _mregion;
        int64_t _vregion;

        int64_t addCacheHit{0};
        int64_t mulCacheHit{0};

        // worker local
        LockFreeQueue<Job*> _queue;
        //SemQueue _queue;

        std::chrono::microseconds timer;
        int executed{0};
        

};

template<>
inline int64_t& Worker::get_region<mNode>(){
    return _mregion;
}
template<>
inline int64_t& Worker::get_region<vNode>(){
    return _vregion;
}

struct Query{
    mEdge lhs;
    mEdge rhs;
    int32_t current_var;
    bool available{false};
    mEdge result;

    bool removed{false};

    Query(const mEdge& l, const mEdge& r, int32_t var):lhs(l),rhs(r), current_var(var){}

    Query(const Query& other): lhs(other.lhs), rhs(other.rhs), current_var(other.current_var){}

    ~Query(){}
    inline bool operator==(const Query& other) const noexcept {
        return ((lhs == other.lhs && rhs == other.rhs) || (lhs == other.rhs && rhs == other.lhs)) && (current_var && other.current_var);
    }

    void set_result(Worker* w ,const mEdge& r){

       result = {r.w, r.n};
       __atomic_store_n(&available, true, __ATOMIC_RELEASE);
    }

    bool load_result(Worker* w,mEdge& r){
        bool a =  __atomic_load_n(&available, __ATOMIC_ACQUIRE);
        if(a){
            r.w = result.w;
            r.n = result.n;
        }
        return a;
    }

    bool is_removed() const {
        bool r = __atomic_load_n(&removed, __ATOMIC_ACQUIRE);
        return r;
    }
};

template<>
struct std::hash<Query>{
    std::size_t operator()(const Query& q) const noexcept {
        std::size_t h1 = std::hash<mEdge>()(q.lhs);
        std::size_t h2 = std::hash<mEdge>()(q.rhs);
        std::size_t h3 = std::hash<int32_t>()(q.current_var);
        std::size_t h = h1 + h2 + h3; 
        return h;
    }
};





class Engine {
    public:
    
        Engine(std::size_t workers, QubitCount q)
            :_total_worker(workers),  _current_worker(0), _stop(false){ 
            identityTable.resize(q);
            for(auto i = 0; i < _total_worker; i++) {
                Worker* w = new Worker(this, i+1, &_stop);
                _workers.push_back(w);
                w->run();
            }
        }
        
        template<typename F, typename... Args>
            Job* submit(F&& f, Args&&... args){
               Job* j = new Job( 0, std::forward<F>(f), std::forward<Args>(args)...); 
                
                std::unique_lock lk(_job_queue_lock);
                _all_jobs.emplace(j);

               return j;
            }

            void submit(Job* j){
                std::unique_lock lk(_job_queue_lock);
                _all_jobs.emplace(j);
                return;
            }

            std::optional<Job*> request() {
                return _all_jobs.steal();
            }
        
        mEdge addReduce(std::vector<Job*>& jobs, std::size_t grain_size){

            std::vector<Job*> next_round;

            std::size_t i = 0;
            std::size_t njobs = jobs.size();

            while(njobs >= grain_size){
                njobs = 0;
                auto total_size = jobs.size();
                for(; i < total_size; i += grain_size){
                    next_round.push_back(this->submit(addSerial, jobs, i+0, std::min(total_size,i+grain_size)));
                    njobs++;
                }

                i = std::min(i, jobs.size());
                
                jobs.insert(jobs.end(), next_round.begin(), next_round.end());
                next_round.clear();
            }


            Job* j   = this->submit(addSerial, jobs, i, jobs.size());
            return j -> getResult();

            
        }
        mEdge mulReduce(std::vector<Job*>& jobs, std::size_t grain_size){
            std::vector<Job*> next_round;

            std::size_t i = 0;
            std::size_t njobs = jobs.size();

            while(njobs >= grain_size){
                njobs = 0;
                auto total_size = jobs.size();
                for(; i < total_size; i += grain_size){
                    next_round.push_back(this->submit(mulSerial, jobs, i+0, std::min(total_size,i+grain_size)));
                    njobs++;
                }

                i = std::min(i, jobs.size());
                
                jobs.insert(jobs.end(), next_round.begin(), next_round.end());
                next_round.clear();
            }

            Job* j   = this->submit(mulSerial, jobs, i, jobs.size());
            return j -> getResult();
            
        }

        std::size_t worker_number() const { return _total_worker;}
        void terminate() {
            while(!_all_jobs.empty()){}
            _stop = true;
            for(Worker* w: _workers ){
                w->_thread.join();
            }

            for(Worker* w: _workers){
                std::cout<<"t"<<w->_id<<": " <<w->timer.count()<<" ms, "<<w->executed<<std::endl;
                std::cout<<"addhit: "<<w->addCacheHit<<", mulhit: "<<w->mulCacheHit<<std::endl;
            }
        }

        Job* steal();
    private:

        std::size_t next_worker();
        
        bool _stop;
        std::atomic_long _current_worker;
        const std::size_t _total_worker;
        std::vector<Worker*> _workers;

        std::mutex _job_queue_lock;
        CLQueue<Job*> _all_jobs;
};



using AddTable = LockFreeMap<Query>;
using MulTable = LockFreeMap<Query>;


