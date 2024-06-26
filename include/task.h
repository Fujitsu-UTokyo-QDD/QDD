#pragma once

#include <boost/fiber/fiber.hpp>
#include <boost/fiber/operations.hpp>
#include <boost/fiber/algo/work_stealing.hpp>
#include <boost/fiber/operations.hpp>
#include <boost/fiber/mutex.hpp>
#include <boost/fiber/condition_variable.hpp>
#include <boost/fiber/algo/algorithm.hpp>
#include "wsq.hpp"
#include <thread>
#include <atomic>
#include <vector>
#include "dd.h"
#include "table.hpp"
#include "cache.hpp"

class Scheduler;

namespace boost{
namespace fibers{
namespace algo{

class my_ws : public algorithm {
private:
    static std::atomic< std::uint32_t >                     counter_;
    static std::vector< intrusive_ptr< my_ws > >    schedulers_;

    std::uint32_t                                           id_;
    std::uint32_t                                           thread_count_;
    
#ifdef BOOST_FIBERS_USE_SPMC_QUEUE
    detail::context_spmc_queue                              rqueue_{};
#else
    detail::context_spinlock_queue                          rqueue_{};
#endif

//    WorkStealingQueue<context*>                             rqueue_{2048};
    std::mutex                                              mtx_{};
    std::condition_variable                                 cnd_{};
    bool                                                    flag_{ false };
    bool                                                    suspend_;

    static void init_( std::uint32_t, std::vector< intrusive_ptr< my_ws > > &);

public:

    my_ws( std::uint32_t, bool = false);

    my_ws( my_ws const&) = delete;
    my_ws( my_ws &&) = delete;

    my_ws & operator=( my_ws const&) = delete;
    my_ws & operator=( my_ws &&) = delete;

    void awakened( context *) noexcept override;

    context * pick_next() noexcept override;

    virtual context * steal() noexcept {
        /*
        auto ctx_opt = rqueue_.steal();
        if(ctx_opt){
            return *ctx_opt;
        }else{
            return nullptr;
        }
        */
        return rqueue_.steal();
    }

    bool has_ready_fibers() const noexcept override {
        return ! rqueue_.empty();
    }

    void suspend_until( std::chrono::steady_clock::time_point const&) noexcept override;

    void notify() noexcept override;
};
}}}

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
private:
    void spawn();

    const int _nworkers;

    std::vector<WorkerThread> _workers;

    boost::fibers::condition_variable_any cond_stop;
    boost::fibers::mutex mtx_stop;
};

