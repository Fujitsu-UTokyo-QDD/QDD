#pragma once
#include <vector>
#include <type_traits>
#include "complex.h"
#include "dd.h"
#include "table.hpp"
#include <future>
#include <functional>
#include <thread>
#include <unordered_map>
#include <atomic>
#include "queue.hpp"



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
    public:

        Worker(Engine* eng, std::size_t id,  bool* stop): _eng(eng), _id(id),  _stop(stop), _queue(1024){};

        void run();

        void submit(Job*);

        Complex getComplexFromCache(const double_pair&);
        void returnComplexToCache(Complex&);

        Complex getComplexFromTable(const double_pair&); 

        Index uniquefy(const mNode& n);

    private:
        Engine* _eng;
        bool* _stop;
        std::thread _thread;
        std::size_t _id;

        // worker local
        LockFreeQueue<Job*> _queue;
        ComplexCache ccache;
        ComplexTable ctable;
        

};

struct Query{
    mEdge lhs;
    mEdge rhs;
    int32_t current_var;
    bool available{false};
    mEdge result;
    inline bool operator==(const Query& other) const noexcept {
        return ((lhs == other.lhs && rhs == other.rhs) || (lhs == other.rhs && rhs == other.lhs)) && (current_var && other.current_var);
    }

    void set_result(Worker* w ,const mEdge& r){

       result = {w->getComplexFromCache(r.w.getValuePair()), r.n};
       __atomic_store_n(&available, true, __ATOMIC_RELEASE);
    }

    bool load_result(Worker* w,mEdge& r){
        bool a =  __atomic_load_n(&available, __ATOMIC_ACQUIRE);
        if(a)
            r = {w->getComplexFromCache(result.w.getValuePair()), r.n};
        return a;
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


class ComplexReturner {
public:
    ComplexReturner(Worker* w, const Complex& c): _w(w), _c(c){}

    ~ComplexReturner(){
        assert(_c.r != nullptr && _c.i != nullptr);
        _w->returnComplexToCache(_c);
    }
private:
    Complex _c;
    Worker* _w;
};
class Engine {
    public:
    
        Engine(std::size_t workers, QubitCount q, NodeTable* unique = &uniqueTable)
            :_total_worker(workers), _uniqueTable(unique), _current_worker(0), _stop(false){ 
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
        
        mEdge addReduce(const std::vector<Job*>& jobs, std::size_t grain_size){
            grain_size = std::max(jobs.size()/_total_worker, grain_size);
            std::vector<Job*> next_round;
            
            for(auto i = 0; i < jobs.size(); i += grain_size){
                next_round.push_back(this->submit(addSerial, jobs, i, i+grain_size));
            }

            Job* j   = this->submit(addSerial, next_round, 0, next_round.size());
            return j -> getResult();

            
        }
        mEdge mulReduce(const std::vector<Job*>& jobs, std::size_t grain_size){
            grain_size = std::max(jobs.size()/_total_worker, grain_size);
            std::vector<Job*> next_round;
            
            for(auto i = 0; i < jobs.size(); i += grain_size){
                next_round.push_back(this->submit(mulSerial, jobs, i, i+grain_size));
            }

            Job* j   = this->submit(mulSerial, next_round, 0, next_round.size());
            return j -> getResult();

            
        }

        std::size_t worker_number() const { return _total_worker;}
        void terminate() {
            _stop = true;
        }

        Job* steal();
    private:

        std::size_t next_worker();
        
        bool _stop;
        std::atomic_long _current_worker;
        const std::size_t _total_worker;
        NodeTable* _uniqueTable;
        std::vector<Worker*> _workers;

};
using OpTable = CHashTable<Query>;
