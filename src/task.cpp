
#if defined(__linux__)
#include <sched.h>
#endif

#include "task.h"
#include <chrono>
#include <iostream>
#include <table.hpp>

#include <random>

#include <boost/assert.hpp>
#include <boost/context/detail/prefetch.hpp>

#include "boost/fiber/detail/thread_barrier.hpp"
#include "boost/fiber/type.hpp"

using namespace std::chrono_literals;

namespace boost {
namespace fibers {
namespace algo {

std::atomic<std::uint32_t> my_ws::counter_{0};
std::vector<intrusive_ptr<my_ws>> my_ws::schedulers_{};

void my_ws::init_(std::uint32_t thread_count,
                  std::vector<intrusive_ptr<my_ws>> &schedulers) {
    // resize array of schedulers to thread_count, initilized with nullptr
    std::vector<intrusive_ptr<my_ws>>{thread_count, nullptr}.swap(schedulers);
}

my_ws::my_ws(std::uint32_t thread_count, bool suspend)
    : id_{counter_++}, thread_count_{thread_count}, suspend_{suspend} {
    static boost::fibers::detail::thread_barrier b{thread_count};
    // initialize the array of schedulers
    static std::once_flag flag;
    std::call_once(flag, &my_ws::init_, thread_count_, std::ref(schedulers_));
    // register pointer of this scheduler
    schedulers_[id_] = this;
    b.wait();
}

void my_ws::awakened(context *ctx) noexcept {
    if (!ctx->is_context(type::pinned_context)) {
        ctx->detach();
    }
    rqueue_.push(ctx);
}

context *my_ws::pick_next() noexcept {
    context *victim = rqueue_.pop();
    // std::optional<context*> victim_opt = rqueue_.pop();
    // context* victim{nullptr};
    if (nullptr != victim) {
        // if ( victim_opt) {
        //    victim = *victim_opt;
        boost::context::detail::prefetch_range(victim, sizeof(context));
        if (!victim->is_context(type::pinned_context)) {
            context::active()->attach(victim);
        }
    } else {
        std::uint32_t id = 0;
        std::size_t count = 0, size = schedulers_.size();
        static thread_local std::minstd_rand generator{std::random_device{}()};
        std::uniform_int_distribution<std::uint32_t> distribution{
            0, static_cast<std::uint32_t>(thread_count_ - 1)};
        do {
            do {
                ++count;
                // random selection of one logical cpu
                // that belongs to the local NUMA node
                id = distribution(generator);
                // prevent stealing from own scheduler
            } while (id == id_);
            // steal context from other scheduler
            victim = schedulers_[id]->steal();
        } while (nullptr == victim && count < size);
        if (nullptr != victim) {
            boost::context::detail::prefetch_range(victim, sizeof(context));
            BOOST_ASSERT(!victim->is_context(type::pinned_context));
            context::active()->attach(victim);
        }
    }
    return victim;
}

void my_ws::suspend_until(
    std::chrono::steady_clock::time_point const &time_point) noexcept {
    if (suspend_) {
        if ((std::chrono::steady_clock::time_point::max)() == time_point) {
            std::unique_lock<std::mutex> lk{mtx_};
            cnd_.wait(lk, [this]() { return flag_; });
            flag_ = false;
        } else {
            std::unique_lock<std::mutex> lk{mtx_};
            cnd_.wait_until(lk, time_point, [this]() { return flag_; });
            flag_ = false;
        }
    }
}

void my_ws::notify() noexcept {
    if (suspend_) {
        std::unique_lock<std::mutex> lk{mtx_};
        flag_ = true;
        lk.unlock();
        cnd_.notify_all();
    }
}

} // namespace algo
} // namespace fibers
} // namespace boost

void Scheduler::spawn() {

    for (int i = 0; i < _nworkers; i++) {
        _workers[i]._id = i;
        _workers[i]._sched = this;

        _workers[i]._thread = new std::thread(
            [this](int id) {
                boost::fibers::use_scheduling_algorithm<
                    boost::fibers::algo::my_ws>(this->_nworkers + 1);
                {
                    std::unique_lock<boost::fibers::mutex> lk(this->mtx_stop);
                    this->cond_stop.wait(lk);
                }
            },
            i);
    }

#if defined(__linux__)
    for (int i = 0; i < _nworkers; i++) {
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        if (pthread_setaffinity_np(_workers[i]._thread->native_handle(),
                                   sizeof(cpu_set_t), &cpuset)) {
            std::cout << "pthread_setaffinity_np failed" << std::endl;
            exit(1);
        }
    }
#endif

    boost::fibers::use_scheduling_algorithm<boost::fibers::algo::my_ws>(
        this->_nworkers + 1);
}

Scheduler::Scheduler(int n, int gcfreq) : _nworkers(n), _gcfreq(gcfreq) {

    for (int i = 0; i < _nworkers; i++) {
        _workers.emplace_back(WorkerThread());
    }
    this->spawn();
}

Scheduler::~Scheduler() {
    cond_stop.notify_all();

    for (WorkerThread &w : _workers) {
        w._thread->join();
    }

    std::size_t add_lk = 0;
    std::size_t add_h = 0;
    std::size_t mul_lk = 0;
    std::size_t mul_h = 0;
#ifdef CACHE

    for (AddCache &c : _aCache) {
        add_lk += c.lookups;
        add_h += c.hits;
    }
    for (MulCache &c : _mCache) {
        mul_lk += c.lookups;
        mul_h += c.hits;
    }
#endif
#ifdef CACHE_GLOBAL
    add_lk = _aCache_global.lookups;
    add_h = _aCache_global.hits;
    mul_lk = _mCache_global.lookups;
    mul_h = _mCache_global.hits;
#endif
    std::cout << "add cache hit ratio: " << (1.0 * add_h) / add_lk << std::endl;
    std::cout << "mul cache hit ratio: " << (1.0 * mul_h) / mul_lk << std::endl;
}

void Scheduler::addGate(const mEdge &e) { _gates.emplace_back(e); }

void Scheduler::clearCache() {
#ifdef CACHE
    for (auto _alocal : _aCache) {
        // _alocal.hitRatio();
        _alocal.clearAll();
    }
    for (auto _mlocal : _mCache) {
        //_mlocal.hitRatio();
        _mlocal.clearAll();
    }
#endif
#ifdef CACHE_GLOBAL
    _aCache_global.clearAll();
    _mCache_global.clearAll();
#endif
}

vEdge Scheduler::buildCircuit(vEdge input) {

    vEdge v = input;

    for (auto i = 0; i < _gates.size(); i++) {
        //        if(i%100==0)
        //          std::cout << "### " << i << " ###\n";
        v = mv_multiply(_gates[i], v);

        if (i % _gcfreq == 0 && i) {
            std::cout << "gc" << std::endl;
            //            v.incRef();
            //            vUnique.gc();
            //            v.decRef();
            vNodeTable new_table(NQUBITS);
            makeUniqueForV(v, new_table);
            vUnique = std::move(new_table);
            clearCache();
        }
    }

    return v;
}
mEdge Scheduler::buildUnitary(const std::vector<mEdge> &g) {

    if (g.size() == 0) {
        return mEdge();
    }

    mEdge rhs = g[0];
    for (int i = 1; i < g.size(); i++) {
        // std::cout<<"i: "<<i<<std::endl;
        rhs = mm_multiply(g[i], rhs);
    }
    return rhs;
}

void makeUniqueForV(vEdge &root, CHashTable<vNode> &v) {
    if (root.isTerminal())
        return;

    makeUniqueForV(root.n->children[0], v);
    makeUniqueForV(root.n->children[1], v);

    vNode *n = v.getNode();
    n->v = root.n->v;
    n->children = root.n->children;
    n = v.lookup(n);
    root.n = n;
    return;
}
