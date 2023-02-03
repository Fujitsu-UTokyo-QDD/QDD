#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <functional>
#include <stdio.h>
#include "common.h"
#include "dd.h"
#include <oneapi/tbb/enumerable_thread_specific.h>

using namespace oneapi::tbb;


/*
 * CL_MASK and CL_MASK_R are for the probe sequence calculation.
 * With 64 bytes per cacheline, there are 8 64-bit values per cacheline.
 */
static const uint64_t CL_MASK     = ~(((LINE_SIZE) / 8) - 1);  //  X & CL_MASK = X/8
static const uint64_t CL_MASK_R   = ((LINE_SIZE) / 8) - 1;    // X&CL_MASK_R = X%8
                                                              //
                                                              //

template<typename T, typename Hash = std::hash<T>, typename ValueEqual = std::equal_to<T>>
class CHashTable{
public:



    CHashTable(QubitCount n): _tables{n}{};

    ~CHashTable(){
        _tables.clear();

        for(auto it = _caches.begin(); it != _caches.end(); it++){
            it->available = nullptr;

            while(it->chunkID > 0){
                it->chunks.pop_back();
                it->chunkID--;
            }

            it->chunkIt = it->chunks[0].begin();
            it->chunkEndIt = it->chunks[0].end();

            it->allocationSize = INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR;
            it->allocations = INITIAL_ALLOCATION_SIZE;
        
        }

        
    }


   T* getNode() {
        Cache& c = _caches.local();
        if (c.available != nullptr) {
            T* p   = c.available;
            c.available = p->next;
            return p;
        }

        if (c.chunkIt == c.chunkEndIt) {
            c.chunks.emplace_back(c.allocationSize);
            c.allocations += c.allocationSize;
            c.allocationSize *= GROWTH_FACTOR;
            c.chunkID++;
            c.chunkIt    = c.chunks[c.chunkID].begin();
            c.chunkEndIt = c.chunks[c.chunkID].end();
        }

        auto p = &(*c.chunkIt);
        ++c.chunkIt;
        return p;
    }

    void returnNode(T* p) {
        if constexpr(std::is_same_v<T, mNode>){
            if(p == mNode::terminal) return; 
        }
        Cache& c = _caches.local();
        p->next   = c.available;
        c.available = p;
    } 


    T* lookup(T* node){
        const auto key = Hash()(*node) % NBUCKETS;
        const Qubit v = node->v;


        T* current = _tables[v]._table[key];
        T* previous = current;
RELOAD:
        while(current != nullptr){
            if(ValueEqual()(*node, *current)){
                assert(current -> v == node->v);

                returnNode(node);

                return current;
            }



            previous = current;
            current = current ->next;
        }

        if(previous == nullptr){
            if(cas(_tables[v]._table+key, nullptr, node)){
                return node;
            }else{
                current = _tables[v]._table[key];
                goto RELOAD;
            }
        }else{
            if(cas(&(previous->next), nullptr, node )){
                return node;
            }else{
               current = previous->next;
               goto RELOAD;
            }
        }


    }


    
    


private:


    struct Table{
        T* _table[NBUCKETS] = {nullptr};
    };

    std::vector<Table> _tables;
    
    struct Cache {
        T*                                available{};
        std::vector<std::vector<T>>       chunks{1, std::vector<T>{INITIAL_ALLOCATION_SIZE}};
        std::size_t                          chunkID{0};
        typename std::vector<T>::iterator chunkIt{chunks[0].begin()};
        typename std::vector<T>::iterator chunkEndIt{chunks[0].end()};
        std::size_t                          allocationSize{INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR};
        std::size_t                       allocations = INITIAL_ALLOCATION_SIZE;
    };

    enumerable_thread_specific<Cache> _caches;


};





using mNodeTable = CHashTable<mNode>;
extern mNodeTable mUnique;



using vNodeTable = CHashTable<vNode>;
extern vNodeTable vUnique;





