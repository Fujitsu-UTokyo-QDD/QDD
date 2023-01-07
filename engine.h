#pragma once
#include <vector>
#include <type_traits>
#include "complex.h"
#include "dd.h"
#include "table.hpp"
#include <future>
#include <functional>


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


class Worker{
    public:

        Complex getComplexFromCache(const double_pair&);
        void returnComplexToCache(Complex&&);

        Complex getComplexFromTable(const double_pair&); 

        Index uniquefy(const mNode& n);

    private:
        // worker local
        JobQueue _queue;
        ComplexCache ccache;
        ComplexTable ctable;
        

};


class Engine {
    public:
    private:
        NodeTable* _uniqueTable;
        std::vector<Worker*> _workers;


};
