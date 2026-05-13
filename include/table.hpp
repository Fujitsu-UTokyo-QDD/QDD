#pragma once

#include "common.h"
#include "dd.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <stdio.h>
#include <type_traits>
#include <utility>

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
/**
 * @brief Unique table and node allocator for decision diagram nodes.
 *
 * The table interns equivalent nodes per qubit variable and reuses returned
 * nodes through an internal free list. This keeps canonical decision diagram
 * nodes shared across operations.
 */
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
            _cache.addChunk(_cache.allocationSize);
            _cache.allocationSize *= GROWTH_FACTOR;
        }

        auto p = _cache.chunkIt;
        std::allocator_traits<typename Cache::Allocator>::construct(
            _cache.allocator, p);
        ++_cache.chunks[_cache.chunkID].constructed;
        p->next = nullptr;
        ++_cache.chunkIt;
        return p;
    }

    void returnNode(T *p) {
        if constexpr (std::is_same_v<T, mNode>) {
            if (p == mNode::terminal)
                return;
        } else if constexpr (std::is_same_v<T, vNode>) {
            if (p == vNode::terminal)
                return;
        }

        p->v = -2;

        p->next = _cache.available;
        _cache.available = p;
    }

    T *register_wo_lookup(T *node){
        // Assuming single threaded
        const auto key = bucket_index(Hash()(*node));
        const Qubit v = node->v;

        node->next = _tables[v]._table[key].node;
        _tables[v]._table[key].node = node;
        return node;
    }

    T *lookup(T *node) {
        const auto key = bucket_index(Hash()(*node));
        const Qubit v = node->v;      
        T *current = _tables[v]._table[key].node;

        if (current == nullptr) {
            _tables[v]._table[key].node = node;
            return node;
        }

        while (current != nullptr) {
            if (ValueEqual()(*node, *current)) {
                assert(current->v == node->v);

                returnNode(node);
                return current;
            }
            current = current->next;
        }

        node->next = _tables[v]._table[key].node;
        _tables[v]._table[key].node = node;
        return node;
    }

    void dump(){
        std::cout << "#chunk = " << _cache.chunkID << std::endl;
    }

    std::size_t get_allocations(){
        return _cache.allocations;
    }

  private:
    static constexpr bool is_power_of_two(std::size_t x) noexcept {
        return x != 0 && (x & (x - 1)) == 0;
    }

    static_assert(is_power_of_two(NBUCKETS),
                  "NBUCKETS must be a power of two");

    static constexpr std::size_t bucket_index(std::size_t h) noexcept {
        return h & (NBUCKETS - 1);
    }

    struct NodeSlot {
        T *node;
    };
    static_assert(sizeof(NodeSlot) == sizeof(T*));
    
    struct Table {
        NodeSlot _table[NBUCKETS] = {nullptr};
    };

    std::vector<Table> _tables;

    struct Cache {
        using Allocator = std::allocator<T>;
        using AllocatorTraits = std::allocator_traits<Allocator>;

        struct Chunk {
            T *data{nullptr};
            std::size_t capacity{0};
            std::size_t constructed{0};
        };

        T *available{};
        std::vector<Chunk> chunks{};
        std::size_t chunkID{0};
        T *chunkIt{nullptr};
        T *chunkEndIt{nullptr};
        std::size_t allocationSize{INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR};
        std::size_t allocations{0};
        Allocator allocator{};

        Cache() {
            addChunk(INITIAL_ALLOCATION_SIZE);
        }

        ~Cache() {
            clear();
        }

        Cache(const Cache &) = delete;
        Cache &operator=(const Cache &) = delete;

        Cache(Cache &&other) noexcept {
            moveFrom(std::move(other));
        }

        Cache &operator=(Cache &&other) noexcept {
            if (this != &other) {
                clear();
                moveFrom(std::move(other));
            }
            return *this;
        }

        void addChunk(std::size_t capacity) {
            T *data = AllocatorTraits::allocate(allocator, capacity);
            Chunk chunk{data, capacity, 0};
            try {
                chunks.push_back(chunk);
            } catch (...) {
                AllocatorTraits::deallocate(allocator, data, capacity);
                throw;
            }

            allocations += capacity;
            chunkID = chunks.size() - 1;
            chunkIt = data;
            chunkEndIt = data + capacity;
        }

      private:
        void clear() noexcept {
            for (auto &chunk : chunks) {
                if (chunk.data == nullptr) {
                    continue;
                }

                if constexpr (!std::is_trivially_destructible_v<T>) {
                    for (std::size_t i = 0; i < chunk.constructed; ++i) {
                        AllocatorTraits::destroy(allocator, chunk.data + i);
                    }
                }
                AllocatorTraits::deallocate(allocator, chunk.data,
                                            chunk.capacity);
            }

            chunks.clear();
            available = nullptr;
            chunkID = 0;
            chunkIt = nullptr;
            chunkEndIt = nullptr;
            allocations = 0;
        }

        void moveFrom(Cache &&other) noexcept {
            available = other.available;
            chunks = std::move(other.chunks);
            chunkID = other.chunkID;
            chunkIt = other.chunkIt;
            chunkEndIt = other.chunkEndIt;
            allocationSize = other.allocationSize;
            allocations = other.allocations;
            allocator = std::move(other.allocator);

            other.available = nullptr;
            for (auto &chunk : other.chunks) {
                chunk = {};
            }
            other.chunks.clear();
            other.chunkID = 0;
            other.chunkIt = nullptr;
            other.chunkEndIt = nullptr;
            other.allocationSize = INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR;
            other.allocations = 0;
        }
    };

    std::size_t collected;
    Cache _cache;

    QubitCount _qn;

};

using mNodeTable = CHashTable<mNode>;
extern mNodeTable mUnique;

using vNodeTable = CHashTable<vNode>;
extern vNodeTable vUnique;
