#include <boost/fiber/operations.hpp>
#include <iostream>
#include "task.h"
#include <functional>
#include <chrono>
#include <thread>
#include <mutex>
#include "algorithms/grover.hpp"
#include <condition_variable>
#include <boost/fiber/future/future.hpp>
#include <boost/fiber/future/packaged_task.hpp>
#include <boost/fiber/buffered_channel.hpp>
using namespace std::chrono_literals;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


static int fiber_counts = 50;
using task = std::packaged_task<int(int)>;
boost::fibers::buffered_channel<task> ch{1024};

int bar(int i){
    for(int j = 0; j < 10000; j++){}
    //boost::this_fiber::sleep_for(std::chrono::seconds(i));
    return i;
}

void foo(int until) {
    std::vector<boost::fibers::future<int>> results;
    for(int i = 0; i < until; i++){
        boost::fibers::packaged_task<int()> pt(std::bind(bar, i)); 
        results.emplace_back(pt.get_future());
        boost::fibers::fiber(std::move(pt)).detach();
    }

    for(int i = 0; i < until; i++){
        std::cout<<i<<": "<<results[i].get()<<std::endl;
    }

}

void work_stealing_test(){
    boost::fibers::condition_variable_any cond;
    boost::fibers::mutex mtx;

    auto t = new std::thread([&](){
        boost::fibers::use_scheduling_algorithm<boost::fibers::algo::work_stealing>(2);
        {
            std::unique_lock<boost::fibers::mutex> lk(mtx);
            cond.wait(lk);
        }

        });



    boost::fibers::use_scheduling_algorithm<boost::fibers::algo::work_stealing>(2);

    boost::fibers::condition_variable_any cond_fin;
    boost::fibers::mutex mtx_fin;
    for(int i = 0; i < fiber_counts; i++){
        boost::fibers::fiber f([&](){
                                   std::cout<<std::this_thread::get_id()<<std::endl;
                                   std::unique_lock<boost::fibers::mutex> lk(mtx_fin);
                                   if(--fiber_counts == 0){
                                        lk.unlock();
                                        cond_fin.notify_one();
                                   }

                               });
        f.detach();
    }

    //boost::this_fiber::yield();
    {
        std::unique_lock<boost::fibers::mutex> lk(mtx_fin);
        cond_fin.wait(lk, [&](){return 0 == fiber_counts;});
        cond.notify_all();
        t->join();
    }


}


void fiber_channel(){
    
    auto worker1 = std::thread(
        []{
            std::cout<<"worker "<<std::this_thread::get_id()<<std::endl;
            // create pool of fibers
            for (int i=0; i<10; ++i) {
                boost::fibers::fiber{
                    []{
                        task tsk;
                        // dequeue and process tasks
                        while (boost::fibers::channel_op_status::closed!=ch.pop(tsk)){
                            tsk(2);
                        }
                    }}.detach();
            }
            task tsk;
            // dequeue and process tasks
            while (boost::fibers::channel_op_status::closed!=ch.pop(tsk)){
                tsk(2);
            }
        });
    auto worker2 = std::thread(
        []{
            std::cout<<"worker "<<std::this_thread::get_id()<<std::endl;
            // create pool of fibers
            for (int i=0; i<10; ++i) {
                boost::fibers::fiber{
                    []{
                        int c = 0;
                        task tsk;
                        // dequeue and process tasks
                        while (boost::fibers::channel_op_status::closed!=ch.pop(tsk)){
                            tsk(4);
                            if(++c < 20){
                                auto pt = std::packaged_task<int(int)>([](int i ){
                                    std::cout<<"new task "<< i <<std::endl; return i--;});    
                                ch.push(std::move(pt));
                            }
                        }
                    }}.detach();
            }
            int c = 0;
            task tsk;
            // dequeue and process tasks
            while (boost::fibers::channel_op_status::closed!=ch.pop(tsk)){
                tsk(4);
            if(++c < 20){
                auto pt = std::packaged_task<int(int)>([](int i ){
                    std::cout<<"new task "<< i <<std::endl; return i--;});    
                ch.push(std::move(pt));
            }
            }
        });
    
    std::vector<std::future<int>> results;
    for(int i = 0; i < 100; i++){
        auto pt = std::packaged_task<int(int)>([](int i){
                    //std::cout<<std::this_thread::get_id()<<" "<<i<<std::endl;
                    return i;});
        results.emplace_back(pt.get_future());
        ch.push(std::move(pt));
        
    }



    std::cout<<"result"<<std::endl;
    int i =0;
    for(std::future<int>& r: results){
        //std::cout<<r.get()<<std::endl;
        i++;
    }
    std::cout<<i<<std::endl;
    ch.close();
    sleep(3);
    worker1.join();
    worker2.join();

}


int main(int argc, char* argv[]){
    const int nworkers = std::stoi(argv[1]);
    const int nqubits = std::stoi(argv[2]);
    std::cout<<"run with "<<nworkers<<" workers, "<<nqubits<<" qubits"<<std::endl;

    Scheduler s(nworkers);
    auto output = groverFiber(s, nqubits);
    std::cout<<"output"<<std::endl;
    output.printVector();
    return 0;


    s.addGate(RZ(15,0,1.7952706710012407));
    s.addGate(RY(15, 0, 1.0056905557557458));
    s.addGate(RZ(15,0,-2.860782987649066));
    s.addGate(RY(15,1,2.9482444835445483));
    s.addGate(RZ(15,0,1.7952706710012407));
    s.addGate(RZ(15,0,1.7952706710012407));
    s.addGate(RY(15, 0, 1.0056905557557458));
    s.addGate(RZ(15,0,-2.860782987649066));
    s.addGate(RY(15,1,2.9482444835445483));
    s.addGate(RY(15, 0, 1.0056905557557458));
    s.addGate(RZ(15,0,-2.860782987649066));
    s.addGate(RY(15,1,2.9482444835445483));

    auto t1 = std::chrono::high_resolution_clock::now();
    vEdge result = s.buildCircuit(makeZeroState(15));
    auto t2 = std::chrono::high_resolution_clock::now();
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" micros"<<std::endl;
    return 0;
}
