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
#include "oneapi/tbb/parallel_reduce.h"
#include "oneapi/tbb/blocked_range.h"

using namespace oneapi::tbb;


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
        std::atomic_llong  _top, _bottom;

};

class Engine;

class Worker{
    public:

        Worker(Engine* eng, std::size_t id,  bool* stop): _eng(eng), _id(id),  _stop(stop){};

        void run();

        void submit(Job*);

        Complex getComplexFromCache(const double_pair&);
        void returnComplexToCache(Complex&&);

        Complex getComplexFromTable(const double_pair&); 

        Index uniquefy(const mNode& n);

    private:
        Engine* _eng;
        bool* _stop;
        std::thread _thread;
        std::size_t _id;

        // worker local
        JobQueue _queue;
        ComplexCache ccache;
        ComplexTable ctable;
        

};

struct AddReduce {
   const std::vector<Job*>& _jobs;
   Engine* _eng;
   mEdge sum;

   AddReduce(const std::vector<Job*>& j, Engine* e):_jobs(j), _eng(e){};
   AddReduce(const AddReduce& other, split): _jobs(other._jobs), _eng(other._eng){};
   void operator()(const blocked_range<size_t>& r);
   void join(const AddReduce& other);
};

struct MulReduct {

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

            for(auto i = 0; i < next_round.size(); i += grain_size){
                assert(next_round[i]->valid());
                next_round[i]->getResult();
            }
            std::cout<<"before"<<std::endl;
            
            mEdge result = add(_workers[0], next_round[0]->getResult(), next_round[1]->getResult());
            std::cout<<"after"<<std::endl;
            for(auto i = 2; i < next_round.size(); i++){
                std::cout<<i<<std::endl;
                result = add(_workers[next_worker()], result, next_round[i]->getResult());
            
            }

            //mEdge result = addSerial(_workers[0], next_round, 0, next_round.size());
            //Job* j   = this->submit(addSerial, next_round, 0, next_round.size());
            //return j -> getResult();
            return result;

            
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
