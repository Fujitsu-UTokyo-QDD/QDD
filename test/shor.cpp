#include <boost/fiber/operations.hpp>
#include <iostream>
#include <fstream>
#include <regex>
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
#if OPENMP
#include <omp.h>
#endif
using namespace std::chrono_literals;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


static int fiber_counts = 50;
using ttask = std::packaged_task<int(int)>;
boost::fibers::buffered_channel<ttask> ch{1024};

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
                        ttask tsk;
                        // dequeue and process tasks
                        while (boost::fibers::channel_op_status::closed!=ch.pop(tsk)){
                            tsk(2);
                        }
                    }}.detach();
            }
            ttask tsk;
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
                        ttask tsk;
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
            ttask tsk;
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
    #if OPENMP
    omp_set_num_threads(std::atoi(argv[1]));
    std::cout << "OMP threads = " << omp_get_max_threads() << std::endl;
#endif

    const int nworkers = std::stoi(argv[1]);
    const int gcfreq = std::stoi(argv[2]);
    std::cout<<"run with "<<nworkers<<" workers, "<< std::endl;

    Scheduler s(nworkers, gcfreq);
    
    std::ifstream reading_file;
    std::string filename = "data/shor253.qasm";
    reading_file.open(filename, std::ios::in);
    std::string line;
    int count = 0;
    while (std::getline(reading_file, line))
    {
        count++;
        if (line.starts_with("OPENQASM 2.0;" || line.starts_with("include qelib1.inc;")) || line.starts_with("qreg q[34];"))
        {
            continue;
        }
        else if (line.starts_with("x "))
        {
            std::smatch m;
            if(std::regex_match(line,m,std::regex(R"(x q\[(\d+)\];)"))){
                s.addGate(makeGate(34, Xmat, std::stoi(m[1].str())));
                //std::cout << "X " << std::stoi(m[1].str()) << std::endl;
            }
        }
        else if (line.starts_with("h "))
        {
            std::smatch m;
            if(std::regex_match(line,m,std::regex(R"(h q\[(\d+)\];)"))){
                s.addGate(makeGate(34, Hmat, std::stoi(m[1].str())));
                //std::cout << "H " << std::stoi(m[1].str()) << std::endl;
            }
        }
        else if (line.starts_with("cx "))
        {
            std::smatch m;
            if(std::regex_match(line,m,std::regex(R"(cx q\[(\d+)\],q\[(\d+)\];)"))){
                s.addGate(CX(34, std::stoi(m[2].str()), std::stoi(m[1].str())));
                //std::cout << "CX " << std::stoi(m[1].str()) << "," << std::stoi(m[2].str())<< std::endl;
            }
        }
        else if (line.starts_with("DenseMatrix(1,0,"))
        {
            std::smatch m;
            if(std::regex_match(line,m,std::regex(R"(DenseMatrix\(1,0,(.+),(.+),(.+),(.+),(.+),(.+),(.+),(.+)\) q\[(\d+)\];)"))){
                s.addGate(Dense(34,std::stoi(m[9].str()),std::stof(m[1].str()),std::stof(m[2].str()),std::stof(m[3].str()),std::stof(m[4].str()),std::stof(m[5].str()),std::stof(m[6].str()),std::stof(m[7].str()),std::stof(m[8].str())));
                //std::cout << "Dense " << std::stof(m[1].str()) << "," << std::stof(m[2].str()) << "," << std::stof(m[3].str()) << "," << std::stof(m[4].str()) << std::endl << "      " << std::stof(m[5].str()) << "," << std::stof(m[6].str()) << "," << std::stof(m[7].str()) << "," << std::stof(m[8].str()) << std::endl;
                //std::cout << "      " << std::stoi(m[9].str()) << std::endl;
            }
        }
        else if (line.starts_with("DenseMatrix(1,2,"))
        {
            std::smatch m;
            if(std::regex_match(line,m,std::regex(R"(DenseMatrix\(1,2,(.+),(.+),(.+),(.+),(.+),(.+),(.+),(.+),1,1\) q\[(\d+)\],q\[(\d+)\],q\[(\d+)\];)"))){
                s.addGate(Dense1211(34,std::stoi(m[9].str()),std::stoi(m[10].str()),std::stoi(m[11].str()),std::stof(m[1].str()),std::stof(m[2].str()),std::stof(m[3].str()),std::stof(m[4].str()),std::stof(m[5].str()),std::stof(m[6].str()),std::stof(m[7].str()),std::stof(m[8].str())));
//                std::cout << "Dense " << std::stof(m[1].str()) << "," << std::stof(m[2].str()) << "," << std::stof(m[3].str()) << "," << std::stof(m[4].str()) << std::endl << "      " << std::stof(m[5].str()) << "," << std::stof(m[6].str()) << "," << std::stof(m[7].str()) << "," << std::stof(m[8].str()) << std::endl;
  //              std::cout << "      " << std::stoi(m[9].str()) << "," <<std::stoi(m[10].str()) << "," << std::stoi(m[11].str()) << std::endl;
            }
        }else{
            assert(1);
        }
    }
    std::cout << count << " gates" << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    vEdge result = s.buildCircuit(makeZeroState(34));
    auto t2 = std::chrono::high_resolution_clock::now();
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()/1000000<<" seconds"<<std::endl;
    return 0;
}
