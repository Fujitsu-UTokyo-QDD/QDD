#pragma once

#include <deque>
#include <semaphore>
#include <chrono>

class Job;

class SemQueue{

    public:
        SemQueue():_sem{0}{}
        bool push(Job*  j){
            _jobs.push_back(j);
            _sem.release();
            return true;
        }
        bool pop( Job* &result){
            using namespace std::chrono_literals; 
            if(!_sem.try_acquire_for(1ms)) return false;
            result = _jobs.front();
            _jobs.pop_front();
            return true;
        
        }
    private:
        std::counting_semaphore<128> _sem;
        std::deque<Job*> _jobs;

};
