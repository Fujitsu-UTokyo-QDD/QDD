#pragma once

#include "common.h"
#include "dd.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <mutex>
#include <random>
#include <stdio.h>
#include <mutex>

#ifdef isMT
#include <oneapi/tbb/enumerable_thread_specific.h>
using namespace oneapi::tbb;
#endif

/*
 * CL_MASK and CL_MASK_R are for the probe sequence calculation.
 * With 64 bytes per cacheline, there are 8 64-bit values per cacheline.
 */
static const uint64_t CL_MASK = ~(((LINE_SIZE) / 8) - 1); //  X & CL_MASK = X/8
static const uint64_t CL_MASK_R = ((LINE_SIZE) / 8) - 1;  // X&CL_MASK_R = X%8
                                                          //
                                                          //

template <typename T, typename Hash = std::hash<T>,
          typename ValueEqual = std::equal_to<T>>
class CHashTable {
  public:
    CHashTable(QubitCount n) : _tables{n}, _qn(n){};

    QubitCount getQubitCount() const { return _qn; }

    T *getNode() {
    #ifdef MT_ENABLE
        const std::lock_guard<std::mutex> lock(table_mtx);
    #endif
#ifdef isMT
        Cache& _cache = _caches.local();
#endif
        if (_cache.available != nullptr) {
            T *p = _cache.available;
            _cache.available = p->next;
            p->next = nullptr;
            return p;
        }

        if (_cache.chunkIt == _cache.chunkEndIt) {
            _cache.chunks.emplace_back(_cache.allocationSize);
            _cache.allocations += _cache.allocationSize;
            _cache.allocationSize *= GROWTH_FACTOR;
            _cache.chunkID++;
            _cache.chunkIt = _cache.chunks[_cache.chunkID].begin();
            _cache.chunkEndIt = _cache.chunks[_cache.chunkID].end();
        }

        auto p = &(*_cache.chunkIt);
        p->next = nullptr;
        ++_cache.chunkIt;
        return p;
    }

    void returnNode(T *p) {
    #ifdef MT_ENABLE
        const std::lock_guard<std::mutex> lock(table_mtx);
    #endif
        if constexpr (std::is_same_v<T, mNode>) {
            if (p == mNode::terminal)
                return;
        }

        p->v = -2;
#ifdef isMT
        Cache& _cache = _caches.local();
#endif

        p->next = _cache.available;
        p->previous = nullptr;
        _cache.available = p;
    }

    T *register_wo_lookup(T *node){
    #ifdef MT_ENABLE
        const std::lock_guard<std::mutex> lock(table_mtx);
    #endif
        const auto key = Hash()(*node) % NBUCKETS;
        const Qubit v = node->v;

         T *current = _tables[v]._table[key];
        //T *previous = current;

        if (current == nullptr) {
            _tables[v]._table[key] = node;
            return node;
        }

        node->previous = _tables[v]._table[key];
        _tables[v]._table[key]->next = node;
        _tables[v]._table[key] = node;
        return node;
    }

    T *lookup(T *node) {
    #ifdef MT_ENABLE
        const std::lock_guard<std::mutex> lock(table_mtx);
    #endif
        const auto key = Hash()(*node) % NBUCKETS;
        const Qubit v = node->v;

        T *current = _tables[v]._table[key];
        //T *previous = current;

        if (current == nullptr) {
            _tables[v]._table[key] = node;
            return node;
        }

        while (current != nullptr) {
            if (ValueEqual()(*node, *current)) {
                assert(current->v == node->v);

                returnNode(node);

                return current;
            }
            current = current->previous;
        }

        node->previous = _tables[v]._table[key];
        _tables[v]._table[key]->next = node;
        _tables[v]._table[key] = node;
        return node;
    }

    void dump(){
#ifdef isMT
        Cache& _cache = _caches.local();
#endif
        std::cout << "#chunk = " << _cache.chunkID << std::endl;
    }

    std::size_t get_allocations(){
#ifdef isMT
        Cache& _cache = _caches.local();
#endif
        return _cache.allocations;
    }

  private:
    struct Table {
        T *_table[NBUCKETS] = {nullptr};
    };

    std::vector<Table> _tables;

    struct Cache {
        T *available{};
        std::vector<std::vector<T>> chunks{
            1, std::vector<T>{INITIAL_ALLOCATION_SIZE}};
        std::size_t chunkID{0};
        typename std::vector<T>::iterator chunkIt{chunks[0].begin()};
        typename std::vector<T>::iterator chunkEndIt{chunks[0].end()};
        std::size_t allocationSize{INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR};
        std::size_t allocations = INITIAL_ALLOCATION_SIZE;
    };

    std::size_t collected;

#ifdef isMT
    enumerable_thread_specific<Cache> _caches;
#else 
    Cache _cache;
#endif

    QubitCount _qn;

#ifdef MT_ENABLE
    std::mutex table_mtx;
#endif

};

using mNodeTable = CHashTable<mNode>;
extern mNodeTable mUnique;

using vNodeTable = CHashTable<vNode>;
extern vNodeTable vUnique;
