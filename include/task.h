#pragma once

#include <boost/fiber/fiber.hpp>
#include <boost/fiber/operations.hpp>
#include <boost/fiber/algo/work_stealing.hpp>
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

struct Scheduler{
    
    Scheduler(int n);
    void spawn();
    ~Scheduler();

    void addGate(const mEdge& e);
    mEdge buildCircuit();
    
    std::atomic_bool _stop;
    std::atomic_bool _ready;
    const int _nworkers;

    std::vector<WorkerThread> _workers;

    std::vector<mEdge> _gates;
};
