#pragma once
#include <vector>
#include <type_traits>
#include "complex.h"
#include "dd.h"
#include <future>
#include <functional>
#include <thread>
#include <unordered_map>
#include <atomic>
#include "queue.hpp"
#include "sem_queue.h"
#include <chrono>



/**
 *  lhs, rhs, op, result
 *
 * */
class Worker;

class Job{
    public:

        template<typename F, typename... Args>
        Job(F&& f, Args&&... args){
            using namespace std::placeholders;
            using return_type = std::invoke_result_t<F, Worker*, Args...>;
            static_assert(std::is_same_v<return_type, mEdge>);

            auto pt = std::packaged_task<mEdge(Worker*)>(std::bind(std::forward<F>(f), _1, std::forward<Args>(args)...));
            task = std::move(pt);
            result = task.get_future();
            
            _worker = nullptr;
        }

        void execute( Worker* w) {
           _worker = w;
           task(_worker);
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
        Worker* _worker;

};

class JobQueue {
    public:

        void push(Job* job);
        Job* pop();
        Job* steal();
        std::size_t size() const;
        bool empty() const;

    private:
        static const unsigned int MAX_JOBS = 1024u;
        static const unsigned int MASK = MAX_JOBS - 1u;

        Job* _jobs[MAX_JOBS];
        long _top{0}, _bottom{0};

};

class Engine;

class Worker{
    friend class Engine;
    public:

        Worker(Engine* eng, std::size_t id,  bool* stop): _eng(eng), _id(id),  _stop(stop), _queue(1024), timer(std::chrono::microseconds::zero()), _mregion(-1), _vregion(-1){};
    //Worker(Engine* eng, std::size_t id,  bool* stop): _eng(eng), _id(id),  _stop(stop), timer(std::chrono::microseconds::zero()), _mregion(-1), _vregion(-1){};

        void run();

        void submit(Job*);

        template<typename T>
        int64_t& get_region();


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
        ComplexCache ccache;

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

struct AddQuery{
    mEdge lhs;
    mEdge rhs;
    int32_t current_var;
    bool available{false};
    mEdge result;

    bool removed{false};

    AddQuery(const mEdge& l, const mEdge& r, int32_t var):lhs(l),rhs(r), current_var(var){}

    AddQuery(const AddQuery& other): lhs(other.lhs), rhs(other.rhs), current_var(other.current_var){}

    ~AddQuery(){}
    inline bool operator==(const AddQuery& other) const noexcept {
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
struct std::hash<AddQuery>{
    std::size_t operator()(const AddQuery& q) const noexcept {
        std::size_t h1 = std::hash<mEdge>()(q.lhs);
        std::size_t h2 = std::hash<mEdge>()(q.rhs);
        std::size_t h3 = std::hash<int32_t>()(q.current_var);
        std::size_t h = h1 + h2 + h3; 
        return h;
    }
};


struct MulQuery{
    Index lhs;
    Index rhs;
    int32_t current_var;
    bool available{false};
    Index result;

    bool removed{false};


    MulQuery(const Index& l, const Index& r, int32_t var):lhs(l), rhs(r), current_var(var){}

    MulQuery(const MulQuery& other):lhs(other.lhs), rhs(other.rhs), current_var(other.current_var){}

    ~MulQuery(){ }

    inline bool operator==(const MulQuery& other) const noexcept {
        return ((lhs == other.lhs && rhs == other.rhs) || (lhs == other.rhs && rhs == other.lhs)) && (current_var && other.current_var);
    }

    void set_result(Worker* w ,const Index& r){
        result = r;
       __atomic_store_n(&available, true, __ATOMIC_RELEASE);
    }

    bool load_result(Worker* w,Index& r){
        bool a =  __atomic_load_n(&available, __ATOMIC_ACQUIRE);
        if(a){
            r = result;
        }
        return a;
    }


    bool is_removed() const {
        bool r = __atomic_load_n(&removed, __ATOMIC_ACQUIRE);
        return r;
    }
};

template<>
struct std::hash<MulQuery>{
    std::size_t operator()(const MulQuery& q) const noexcept {
        std::size_t h1 = std::hash<Index>()(q.lhs);
        std::size_t h2 = std::hash<Index>()(q.rhs);
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
                Worker* w = new Worker(this, i, &_stop);
                _workers.push_back(w);
                w->run();
            }
        }
        
        template<typename F, typename... Args>
            Job* submit(F&& f, Args&&... args){
               Job* j = new Job( std::forward<F>(f), std::forward<Args>(args)...); 

               auto next = next_worker();

               _workers[next]->submit(j);

               return j;
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

};



using AddTable = LockFreeMap<AddQuery>;
using MulTable = LockFreeMap<MulQuery>;


