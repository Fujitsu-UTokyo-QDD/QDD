#pragma once
#include <new>
#include <vector>
#include "common.h"
#include <algorithm>
#include <type_traits>
#include "dd.h"
#include <cassert>
#include <random>


#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_constructive_interference_size;
    using std::hardware_destructive_interference_size;
#else
    // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │ ...
    constexpr std::size_t hardware_constructive_interference_size = 64;
    constexpr std::size_t hardware_destructive_interference_size = 64;
#endif


class AddCache{
    public:
        AddCache(QubitCount q):_tables{2}, rng(std::random_device()()), dist(0,1){

            assert(_tables.size() == 2);
            for(auto i = 0; i < _tables.size(); i++){
                _tables[i].resize(q);
            }

        }

        template<typename T>
            T find(T lhs, T rhs){

                if(lhs.getVar() > rhs.getVar()){
                    swap(lhs, rhs);
                }

                //pick the right table
                int idx = 0;
                if constexpr(std::is_same_v<T, mEdge>){ //mm
                    idx = 0;
                }else if constexpr(std::is_same_v<T, vEdge>){ //vv
                    idx = 1; 
                }else{
                    assert(false);
                }

                Qubit lv = lhs.getVar();
                Table& t = _tables[idx][lv];
                

                std::size_t key = this->hash(lhs, rhs) % NBUCKETS; 

                return find_in_bucket<T>(t, key ,lhs, rhs);

            }
        template<typename T>
            void set(T lhs, T rhs, const T& result){

                if(lhs.getVar() > rhs.getVar()){
                    swap(lhs, rhs);
                }

                //pick the right table
                int idx = 0;
                if constexpr(std::is_same_v<T, mEdge>){ //mm
                    idx = 0;
                }else if constexpr(std::is_same_v<T, vEdge>){ //vv
                    idx = 1; 
                }else{
                    assert(false);
                }

                Qubit lv = lhs.getVar();
                Table& t = _tables[idx][lv];
                

                std::size_t key = this->hash(lhs, rhs) % NBUCKETS; 

                return set_in_bucket<T>(t, key ,lhs, rhs, result);
                
            }


    private:


        union Edge{
            mEdge m;
            vEdge v;
            Edge(){}
        };

        struct Entry{
            Edge lhs;
            Edge rhs;
            Edge result;
        };

        struct alignas(hardware_constructive_interference_size) // the same cacheline
        Bucket{
            Entry e;
        };



        static_assert(std::is_default_constructible_v<Bucket>);
        struct Cache {
            std::vector<std::vector<Bucket>>       chunks{1, std::vector<Bucket>{INITIAL_ALLOCATION_SIZE}};
            std::size_t                          chunkID{0};
            typename std::vector<Bucket>::iterator chunkIt{chunks[0].begin()};
            typename std::vector<Bucket>::iterator chunkEndIt{chunks[0].end()};
            std::size_t                          allocationSize{INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR};
            std::size_t                       allocations = INITIAL_ALLOCATION_SIZE;
        };

        Cache c;

        struct Table{
            Bucket* _table[NBUCKETS] = {nullptr};

        };


        std::vector<std::vector<Table>> _tables; //mm,  vv
                                                 //
        std::mt19937_64 rng;
        std::uniform_int_distribution<int> dist;

        Bucket* getBucket() {

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

        template<typename T>
            std::size_t hash(const T& lhs, const T& rhs){
                assert(lhs.getVar() <= rhs.getVar());
                auto h1 = std::hash<T>()(lhs);
                auto h2 = std::hash<T>()(rhs);
                return hash_combine(h1,h2);
            }

        template<typename T>
            T find_in_bucket(Table& t, std::size_t key, const T& l, const T& r){
                if(t._table[key] == nullptr) {t._table[key] = getBucket(); return T{};}
                Bucket* b = t._table[key];

                Entry& e = b->e;

                if constexpr(std::is_same_v<T, mEdge>){
                    if(e.lhs.m == l && e.rhs.m ==r){
                       return e.result.m; 
                    }else{
                        return T{};
                    }
                }else{
                    if(e.lhs.v == l && e.rhs.v ==r){
                       return e.result.v; 
                    }else{
                        return T{}; 
                    }
                
                }
            }

        template<typename T>
            void set(Entry& e, const T& l , const T& r, const T& res){
                    if constexpr(std::is_same_v<T, mEdge>){
                        e.lhs.m = l;
                        e.rhs.m = r;
                        e.result.m = res;
                    }else{
                        e.lhs.v = l;
                        e.rhs.v = r;
                        e.result.v = res;
                    }
            }

        template<typename T>
            void set_in_bucket(Table& t, std::size_t key, const T& l, const T& r, const T& result){
                if(t._table[key] == nullptr) t._table[key] = getBucket();
                Bucket* b = t._table[key];

                Entry& e = b->e;
                return set(e, l, r, result);
            }


};

class MulCache{

    public:
        MulCache(QubitCount q):_tables{3}, rng(std::random_device()()), dist(0,1){

            assert(_tables.size() == 3);
            for(auto i = 0; i < _tables.size(); i++){
                _tables[i].resize(q);
            }

        }

        template<typename LT, typename RT, typename RetT = std::conditional_t<std::is_same_v<RT, vNode>, vEdge, mEdge>>
            RetT find(const LT* lhs, const RT* rhs){
                //pick the right table
                int idx = 0;
                if constexpr(std::is_same_v<RT, mNode>){ //mm
                    idx = 0;
                }else if constexpr(std::is_same_v<LT, mNode>){ //mv
                    idx = 1; 
                }else if constexpr(std::is_same_v<LT, vNode>){ //vv
                    idx = 2; 
                }else{
                    assert(false);
                }

                Qubit lv = lhs->v;
                Table& t = _tables[idx][lv];
                
                uintptr_t l = reinterpret_cast<uintptr_t>(lhs);
                uintptr_t r = reinterpret_cast<uintptr_t>(rhs);

                std::size_t key = murmur_hash(l) % NBUCKETS; 

                return find_in_bucket<RetT>(t, key ,l, r);

            }

        template<typename LT, typename RT, typename ResT = std::conditional_t<std::is_same_v<RT, vNode>, vEdge, mEdge>>
            void set(const LT* lhs, const RT* rhs, const ResT& result){
                //pick the right table
                int idx = 0;
                if constexpr(std::is_same_v<RT, mNode>){ //mm
                    idx = 0;
                }else if constexpr(std::is_same_v<LT, mNode>){ //mv
                    idx = 1; 
                }else if constexpr(std::is_same_v<LT, vNode>){ //vv
                    idx = 2; 
                }else{
                    assert(false);
                }

                Qubit lv = lhs->v;
                Table& t = _tables[idx][lv];
                
                uintptr_t l = reinterpret_cast<uintptr_t>(lhs);
                uintptr_t r = reinterpret_cast<uintptr_t>(rhs);

                std::size_t key = murmur_hash(l) % NBUCKETS; 

                return set_in_bucket<ResT>(t, key ,l, r, result);
                
            }


    private:


        union Edge{
            mEdge m;
            vEdge v;
            Edge(){}
        };

        struct Entry{
            uintptr_t lhs;
            uintptr_t rhs;
            Edge result;
        };

        struct alignas(hardware_constructive_interference_size) // the same cacheline
        Bucket{
            Entry e1;
            Entry e2;
        };



        static_assert(std::is_default_constructible_v<Bucket>);
        struct Cache {
            std::vector<std::vector<Bucket>>       chunks{1, std::vector<Bucket>{INITIAL_ALLOCATION_SIZE}};
            std::size_t                          chunkID{0};
            typename std::vector<Bucket>::iterator chunkIt{chunks[0].begin()};
            typename std::vector<Bucket>::iterator chunkEndIt{chunks[0].end()};
            std::size_t                          allocationSize{INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR};
            std::size_t                       allocations = INITIAL_ALLOCATION_SIZE;
        };

        Cache c;

        struct Table{
            Bucket* _table[NBUCKETS] = {nullptr};

        };


        std::vector<std::vector<Table>> _tables; //mm, mv, vv
                                                 //
        std::mt19937_64 rng;
        std::uniform_int_distribution<int> dist;

        Bucket* getBucket() {

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

        template<typename RET>
            RET find_in_bucket(Table& t, std::size_t key, uintptr_t l, uintptr_t r){
                if(t._table[key] == nullptr) {t._table[key] = getBucket(); return RET{};}
                Bucket* b = t._table[key];

                Entry& e1 = b->e1;
                Entry& e2 = b->e2;

                if constexpr(std::is_same_v<RET, mEdge>){
                    if(e1.result.m.n != nullptr && e1.lhs == l && e1.rhs == r){
                        return e1.result.m;
                    }else if(e2.result.m.n != nullptr && e2.lhs == l && e2.rhs == r){
                        return e2.result.m;
                    }else{
                        return RET{};
                    }

                }else{
                    if(e1.result.v.n != nullptr && e1.lhs == l && e1.rhs == r){
                        return e1.result.v;
                    }else if(e2.result.m.n != nullptr && e2.lhs == l && e2.rhs == r){
                        return e2.result.v;
                    }else{
                        return RET{};
                    }
                
                }
            }

        template<typename T>
            void set(Entry& e, uintptr_t& l , uintptr_t& r, const T& res){
                    if constexpr(std::is_same_v<T, mEdge>){
                        e.lhs = l;
                        e.rhs = r;
                        e.result.m = res;
                    }else{
                        e.lhs = l;
                        e.rhs = r;
                        e.result.v = res;
                    }
            }

        template<typename ResT>
            void set_in_bucket(Table& t, std::size_t key, uintptr_t l, uintptr_t r, const ResT& result){
                if(t._table[key] == nullptr) t._table[key] = getBucket();
                Bucket* b = t._table[key];

                Entry& e1 = b->e1;
                Entry& e2 = b->e2;


                if constexpr(std::is_same_v<ResT, mEdge>){
                    if(e1.result.m.n == nullptr){
                        set(e1, l, r, result);
                        return;
                    }else if(e2.result.m.n == nullptr ){
                        set(e2, l, r, result);
                        return;
                    }else{
                        //currently just evict randomly
                       if(dist(rng) == 0) return set(e1, l, r, result);
                       else return set(e2, l, r, result);
                    }

                }else{
                    if(e1.result.v.n == nullptr){
                        set(e1, l, r, result);
                        return;
                    }else if(e2.result.m.n == nullptr){
                        set(e2, l, r, result);
                        return;
                    }else{
                       if(dist(rng) == 0) return set(e1, l, r, result);
                       else return set(e2, l, r, result);
                    }
                
                }
            }


};



