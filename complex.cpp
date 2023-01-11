#include "complex.h"
#include "common.h"
#include <algorithm>
#include <cstdio>

Value Value::zero{0.0};
Value Value::one{1,0};

std::size_t ComplexCache::ALLOCATION_SIZE = 1024;
std::size_t ComplexCache::GROWTH = 2;

Complex Complex::one = Complex(&Value::one, &Value::zero);
Complex Complex::zero = Complex(&Value::zero, &Value::zero);

Complex ComplexCache::getCached() {

    return this->getCachedComplex();

}


Complex ComplexCache::getCached_0() {
    Complex c = this->getCached();
    c.r->v = 0.0;
    c.i->v = 0.0;
    return c;

}

Complex ComplexCache::getCached_1() {
    Complex c = this->getCached();
    c.r->v = 1.0;
    c.i->v = 0.0;
    return c;



}

Complex ComplexCache::getCached_v(const double_pair& p) {
    Complex c = this->getCached();
    c.r->v = p.first;
    c.i->v = p.second;
    return c;

}


/*
Value* ComplexTable::find_or_insert(double val){
    Value* v = new Value{ .v = val, .next = nullptr};

    auto t = map.try_emplace(val, v);

    if(t.second){
        return v;
    }else{
        delete v;
        return t.first->second;
    }

}

Complex ComplexTable::find_or_insert(const double_pair& p){
    Value* r = find_or_insert(p.first);
    Value* i = find_or_insert(p.second);
    
    return {r, i, Complex::Source::Table};
}
*/

Value* ComplexTable::find_or_insert(double val){

    Value* v = _map.find_or_insert(double_to_bits(val), {.v = val});

    return v;
}


Complex ComplexTable::find_or_insert(const double_pair& p){
    Value* r = find_or_insert(p.first);
    Value* i = find_or_insert(p.second);

    return {r, i, Complex::Source::Table};
}

double_pair Complex::add(const Complex& lhs, const Complex& rhs){
    
    return {lhs.r->v + rhs.r->v, lhs.i->v + rhs.i->v};
}
double_pair Complex::sub(const Complex& lhs, const Complex& rhs){
    return {lhs.r->v - rhs.r->v, lhs.i->v - rhs.i->v};
}
double_pair Complex::mul(const Complex& lhs, const Complex& rhs){
    if(lhs.isApproximatelyZero() || rhs.isApproximatelyZero()){
        return {0.0, 0.0};
    }

    if(lhs.isApproximatelyOne()){
        return {rhs.r->v, rhs.i->v};
    }

    if(rhs.isApproximatelyOne()){
        return {lhs.r->v, lhs.i->v};
    }

    const auto lr = lhs.r->v;
    const auto li = lhs.i->v;
    const auto rr = rhs.r->v;
    const auto ri = rhs.i->v;

    return { lr * rr - li * ri, lr * ri + li * rr};
}


double_pair Complex::div( const double_pair& lhs, const double_pair& rhs){


    const auto lr = lhs.first;
    const auto li = lhs.second;
    const auto rr = rhs.first;
    const auto ri = rhs.second;

    const auto rmag2 = rr * rr + ri * ri;
    
    const auto r = (lr*rr+li*ri)/rmag2;
    const auto i = (li*rr-lr*ri)/rmag2;

    return {r,i};
}

double_pair Complex::div( const Complex& lhs, const Complex& rhs){
    assert(!rhs.isApproximatelyZero() && rhs.mag2() != 0.0);

    if(lhs.isApproximatelyEqual(rhs)){
        return {1.0, 0.0};
    }

    if(lhs.isApproximatelyZero()){
       return {0.0, 0.0};
    }

    const auto lr = lhs.r->v;
    const auto li = lhs.i->v;
    const auto rr = rhs.r->v;
    const auto ri = rhs.i->v;

    const auto rmag2 = rhs.mag2();

    return {(lr*rr+li*ri)/rmag2,  (li*rr-lr*ri)/rmag2};
}
double_pair Complex::div( const Complex& lhs, const double_pair& rhs){
    double rmag2 = rhs.first * rhs.first + rhs.second * rhs.second;
    assert( rmag2 != 0.0);



    const auto lr = lhs.r->v;
    const auto li = lhs.i->v;
    const auto rr = rhs.first;
    const auto ri = rhs.second;


    return {(lr*rr+li*ri)/rmag2,  (li*rr-lr*ri)/rmag2};
}


