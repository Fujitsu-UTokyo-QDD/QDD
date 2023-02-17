#pragma once

#include <boost/fiber/fiber.hpp>
#include <boost/fiber/operations.hpp>
#include <boost/fiber/algo/work_stealing.hpp>
#include <boost/fiber/operations.hpp>
#include <boost/fiber/mutex.hpp>
#include <boost/fiber/condition_variable.hpp>
#include <thread>
#include <atomic>
#include <vector>
#include "dd.h"

struct Scheduler;

struct WorkerThread{
    int _id;
    Scheduler* _sched; 
    std::thread* _thread;
};

class Scheduler{
    friend struct WorkerThread;
    
public:
    Scheduler(int n);
    ~Scheduler();

    void addGate(const mEdge& e);
    vEdge buildCircuit(vEdge v, int gcfreq=1000000);
    mEdge buildUnitary(const std::vector<mEdge>& g);
private:
    void spawn();
    
    const int _nworkers;

    std::vector<WorkerThread> _workers;
    std::vector<mEdge> _gates;

    boost::fibers::condition_variable_any cond_stop;
    boost::fibers::mutex mtx_stop;
};
