#pragma once 
#include <cstdint>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include "common.h"

struct Value {
    double v{0.0};
    Value* next{nullptr};


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

struct Complex {
    enum Source:uint8_t {
        Cache,
        Table,
    };

    Value* r{nullptr};
    Value* i{nullptr};
    
    Source src{Cache};
    
    Complex() = default;
    Complex(Value* rr, Value* ii, Source ss):r(rr),i(ii), src(ss){};

    

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


// thread local
class ComplexCache {
public:
    
    ~ComplexCache(){
        for(auto it = available.begin(); it != available.end(); it++){
            delete *it;
        }
    }

    Complex getCached();
    Complex getCached_1();
    Complex getCached_0();
    Complex getCached_v(const double_pair&);
    void returnCached(Complex&&); 


private:
    
    static const std::size_t ALLOCATION_SIZE = 512;
    static_assert((ALLOCATION_SIZE&0x1) == 0);

    std::vector<Value*> available;
    
    std::pair<Value*, Value*> request();

};

// thread local
class ComplexTable {
    public:
        Complex find_or_insert(const double_pair&);
        ~ComplexTable();
    private:
        Value* find_or_insert(double);
        std::unordered_map<double, Value*> map;
};
