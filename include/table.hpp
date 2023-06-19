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
        if constexpr (std::is_same_v<T, mNode>) {
            if (p == mNode::terminal)
                return;
        }

        p->v = -2;

        p->next = _cache.available;
        _cache.available = p;
    }

    T *lookup(T *node) {
        const auto key = Hash()(*node) % NBUCKETS;
        const Qubit v = node->v;

        T *current = _tables[v]._table[key];
        T *previous = current;

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

            previous = current;
            current = current->next;
        }

        previous->next = node;
        return node;
    }

    void gc() {

        collected = 0;
        for (Table &t : _tables) {
            for (auto i = 0; i < NBUCKETS; i++) {
                bucket_gc(t, i);
            }
        }

        // std::cout<<"gc collected "<<collected<<std::endl;
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

    Cache _cache;

    QubitCount _qn;
};

using mNodeTable = CHashTable<mNode>;
extern mNodeTable mUnique;

using vNodeTable = CHashTable<vNode>;
extern vNodeTable vUnique;
