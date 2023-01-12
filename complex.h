#pragma once 
#include <cstdint>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include "common.h"
#include "lockfree_hashmap.hpp"

struct Value {
    double v{0.0};
    Value* next{nullptr};
    
    static Value zero;
    static Value one;

    static constexpr bool isApproximatelyEqual(double lhs, double rhs, double tol = std::numeric_limits<double>::epsilon()){
    double diff = std::fabs(lhs - rhs);
    if(diff <= tol)
        return true;

    if(diff < std::fmax(std::fabs(lhs), std::fabs(rhs)) * tol)
        return true;

    return false;

    }

    static constexpr bool isApproximatelyEqual(const Value* lhs, const Value* rhs){
        return lhs == rhs || isApproximatelyEqual(lhs->v, rhs->v);
    }

    static constexpr bool isApproximatelyZero(double v, double tol = std::numeric_limits<double>::epsilon()){
    if(std::fabs(v) <= tol)
        return true;
    return false;
    }

    static constexpr bool isApproximatelyZero(const Value* v){
    return isApproximatelyZero(v->v);
    }


    static constexpr bool isApproximatelyOne(double v, double tol = std::numeric_limits<double>::epsilon()){
    return isApproximatelyEqual(v, 1.0);
    }

    static constexpr bool isApproximatelyOne(const Value* v){
    return isApproximatelyOne(v->v);
    }
    
    inline bool operator==(const Value& other) const noexcept {
        return v == other.v;
    }

};

static_assert(std::is_aggregate_v<Value>);

template<>
struct std::hash<Value>{
    std::size_t operator()(const Value& v) const noexcept {
        return std::hash<double>()(v.v);   
    }
};

struct  Complex {
    enum Source:uint8_t {
        Cache,
        Table,
    };

    Value* r{nullptr};
    Value* i{nullptr};
    
    Source src{Cache};

    static Complex one;
    static Complex zero;
    
    Complex() = default;
    Complex(Value* rr, Value* ii, Source ss):r(rr),i(ii), src(ss){};
    Complex(Value* rr, Value* ii):r(rr), i(ii){};

    

    double_pair getValuePair() const {
        assert(r != nullptr && i != nullptr);
        return {r->v, i->v};
    }



    inline bool valueEqual(const Complex& other) const{
        return r->v == other.r->v && i->v == other.i->v;
    }
    inline bool objectEqual(const Complex& other) const{
        return r == other.r && i == other.i;
    }

    inline bool operator==(const Complex& other) const {
        return this->valueEqual(other) && src == other.src;
    }


    inline bool isApproximatelyEqual(const Complex& other) const {
        return Value::isApproximatelyEqual(r, other.r) && Value::isApproximatelyEqual(i, other.i);
    }

    inline bool isApproximatelyZero() const {
        return Value::isApproximatelyZero(r) && Value::isApproximatelyZero(i);
    }

    inline bool isApproximatelyOne() const {
        return Value::isApproximatelyOne(r) && Value::isApproximatelyOne(i);
    }


    static double_pair add(const Complex& lhs, const Complex& rhs);
    static double_pair sub(const Complex& lhs, const Complex& rhs);
    static double_pair mul(const Complex& lhs, const Complex& rhs);
    static double_pair div(const Complex& lhs, const Complex& rhs);
    static double_pair div(const Complex& lhs, const double_pair& rhs);
    static double_pair div(const double_pair& lhs, const double_pair& rhs);


    
    inline double mag() const {
        double rr = r->v;
        double ii = i->v;
        return std::sqrt(rr * rr + ii * ii);
    }
    inline double mag2() const {
        double mag = this->mag();
        return mag * mag;

    }

    inline bool operator<(const Complex& other) const {
        return this->mag2() < other.mag2();
    }
     
};



template<>
struct std::hash<Complex>{
    std::size_t operator()(const Complex& v) const noexcept {
        auto h1 =std::hash<Value>()(*(v.r));
        auto h2 =std::hash<Value>()(*(v.i));
        return hash_combine(h1,h2);
    }
};
inline std::ostream& operator<<(std::ostream& os, const Complex& c){
    const auto rv = c.r->v;
    const auto iv = c.i->v;
    if(iv != 0)
        os << std::fixed<<std::setprecision(4)<<rv<<" + "<<iv<<"i";    
    else
        os << std::fixed<<std::setprecision(4)<<rv;    
    return os;

}





class ComplexCache {
    public:
        ComplexCache(): allocationSize(ALLOCATION_SIZE) {
            chunks.emplace_back(allocationSize);
            allocations += allocationSize;
            allocationSize *= GROWTH;
            chunkIt = chunks[0].begin();
            chunkEndIt = chunks[0].end();
        }

        ~ComplexCache() = default;

        Complex getCachedComplex() {
           if (available != nullptr) {
                assert(available->next != nullptr);
                auto entry = Complex{available, available->next, Complex::Cache};
                available  = entry.i->next;
                count += 2;
                return entry;
            }

            // new chunk has to be allocated
            if (chunkIt == chunkEndIt) {
                chunks.emplace_back(allocationSize);
                allocations += allocationSize;
                allocationSize *= GROWTH;
                chunkID++;
                chunkIt    = chunks[chunkID].begin();
                chunkEndIt = chunks[chunkID].end();
            }

            Complex c{};
            c.r = &(*chunkIt);
            ++chunkIt;
            c.i = &(*chunkIt);
            ++chunkIt;
            count += 2;
            return c; 
        }

        void returnCached(Complex& c) {
            if(c.r == Complex::one.r || c.r == Complex::zero.r)
                return;
            assert(available != c.r && available != c.i);

            c.i->next = available;
            c.r->next = c.i;
            available = c.r;
            count -= 2;

            c.r = c.i = nullptr;
        }

        void clear() {
            available = nullptr;

            while (chunkID > 0) {
                chunks.pop_back();
                chunkID--;
            }
            chunkIt        = chunks[0].begin();
            chunkEndIt     = chunks[0].end();
            allocationSize = ALLOCATION_SIZE * GROWTH;
            allocations    = ALLOCATION_SIZE;

            count     = 0;
            peakCount = 0;
        };

        Complex getCached();
        Complex getCached_0();
        Complex getCached_1();
        Complex getCached_v(const double_pair& p);

    private:
        static std::size_t ALLOCATION_SIZE;
        static std::size_t GROWTH; 

        Value* available;
        std::vector<std::vector<Value>>  chunks;
        std::size_t        chunkID;
        typename std::vector<Value>::iterator  chunkIt;
        typename std::vector<Value>::iterator  chunkEndIt;
        std::size_t             allocationSize;


        std::size_t allocations = 0;
        std::size_t count = 0;
        std::size_t peakCount = 0;

};

// thread local
/*
class ComplexTable {
    public:
        Complex find_or_insert(const double_pair&);
        ~ComplexTable();
    private:
        Value* find_or_insert(double);
        std::unordered_map<double, Value*> map;
};
*/
//global
class ComplexTable {
    public:
        Complex find_or_insert(const double_pair&);
    private:
       Value* find_or_insert(double);
       LockFreeMap<Value> _map; 
};


extern ComplexTable ctable;



