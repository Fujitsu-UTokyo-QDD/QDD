#include "complex.h"
#include <algorithm>


std::pair<Value*, Value*> ComplexCache::request() {
    
    if(available.size() == 0){
        available.resize(ALLOCATION_SIZE);
        for(auto i = 0; i < ALLOCATION_SIZE; i++) available[i] = new Value();
        //std::generate_n(available.begin(), ALLOCATION_SIZE, [](){ return new Value();}); 
    }

    assert(available.size() >= 2);

    Value* r = available.back();
    available.pop_back();
    Value* i = available.back();
    available.pop_back();
    assert(r != nullptr && i !=nullptr);
    return std::make_pair(r, i);


}

Complex ComplexCache::getCached() {

    auto values = request();    

    return {values.first, values.second, Complex::Source::Cache};

}


Complex ComplexCache::getCached_0() {

    auto values = request();
    values.first->v = 0.0;
    values.second->v = 0.0;

    return {values.first, values.second, Complex::Source::Cache};

}

Complex ComplexCache::getCached_1() {

    auto values = request();
    values.first->v = 1.0;
    values.second->v = 0.0;

    return {values.first, values.second, Complex::Source::Cache};

}

Complex ComplexCache::getCached_v(const double_pair& p) {

    auto values = request();
    values.first->v = p.first;
    values.second->v = p.second;

    return {values.first, values.second, Complex::Source::Cache};

}


void ComplexCache::returnCached(Complex && c){
    Value* r = c.r;
    Value* i = c.i;

    r->v = i->v = 0.0;
    r->next = i->next = nullptr;

    available.push_back(r);
    available.push_back(i);

    c.r = c.i = nullptr;
}



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

ComplexTable::~ComplexTable(){
    for(auto it = map.begin(); it != map.end(); it++){
        delete it->second;
    }
}
