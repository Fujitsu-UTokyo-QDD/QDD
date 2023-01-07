#include "common.h"
#include "table.hpp"
#include <iostream>
#include "engine.h"
#include "dd.h"
#include <random>
#include <chrono>


using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main()
{
    Worker w;

    auto t1 = high_resolution_clock::now();
    mEdge e = makeIdent(&w,10);
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms = t2 - t1;
    std::cout<<ms.count()<<" ms"<<std::endl; 
    w.returnComplexToCache(std::move(e.w));
    //e.printMatrix();
    
}

