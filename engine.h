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
#include <random>
#include "wsq.hpp"
#include "cache.hpp"



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

            task = std::packaged_task<mEdge(Worker*)>(std::bind(std::forward<F>(f), _1, std::forward<Args>(args)...));
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

        bool available() const {
            using namespace std::chrono_literals;
            if(result.wait_for(0s) == std::future_status::ready) return true;
            else return false;
        }
    private:
        std::packaged_task<mEdge(Worker*)> task;
        std::future<mEdge> result;

        std::size_t _parent;

};



class Engine;

class Worker{
    friend class Engine;
    public:

        Worker(Engine* eng, std::size_t id,std::size_t total ,  bool* stop, QubitCount q): _eng(eng), _id(id),  _stop(stop), worker_dist(0, total-1), _mulCache(q), _addCache(q) {
            rng.seed(_id);
        };

        void run();


        template<typename F, typename... Args>
            Job* submit(F&& f, Args&&... args){
               Job* j = new Job( _id, std::forward<F>(f), std::forward<Args>(args)...); 
                _local_jobs.push(j);
                return j;
            }

        void run_pending();

        Engine* _eng;
        bool* _stop;
        std::thread _thread;
        std::size_t _id;
        WorkStealingQueue<Job*> _local_jobs;


        std::chrono::microseconds timer{std::chrono::microseconds::zero()};
        int executed{0};

        //pick victim
        std::mt19937_64 rng;
        std::uniform_int_distribution<int> worker_dist;
        
        MulCache _mulCache; 
        AddCache _addCache;

};



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
            for(auto i = 0; i < _total_worker; i++) {
                Worker* w = new Worker(this, i+1, _total_worker, &_stop, q);
                _workers.push_back(w);
            }
            for(auto i = 0; i < _total_worker; i++) {
                _workers[i]->run();
            }
        }
        
        template<typename F, typename... Args>
            Job* submit(F&& f, Args&&... args){
               Job* j = new Job( 0, std::forward<F>(f), std::forward<Args>(args)...); 
               _local_jobs.push(j); 
               return j;
            }


        std::size_t worker_number() const { return _total_worker;}

        void terminate() {
            _stop = true;
            for(Worker* w: _workers ){
                w->_thread.join();
            }


        }

        std::optional<Job*> steal(std::size_t target);
    private:

        
        bool _stop;
        std::atomic_long _current_worker;
        const std::size_t _total_worker;
        std::vector<Worker*> _workers;
        WorkStealingQueue<Job*> _local_jobs;

};



using AddTable = LockFreeMap<Query>;
using MulTable = LockFreeMap<Query>;


