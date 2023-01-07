#pragma once

#include "complex.h"
#include "common.h"

struct mNode;
struct Worker;

struct __attribute__ ((packed)) mEdge {

    Qubit getVar() const;
    bool isTerminal() const;
    mNode* getNode() const;
    
    void printMatrix() const;

    inline bool operator==(const mEdge& e) const noexcept {
        return w == e.w && n == e.n;
    }

    Complex w;
    Index n;
};

static_assert(std::is_aggregate_v<mEdge>);
template<>
struct std::hash<mEdge>{
    std::size_t operator()(const mEdge& v) const noexcept {
        auto h1 = std::hash<Complex>()(v.w);
        auto h2 = std::hash<Index>()(v.n);
        return hash_combine(h1,h2);
    }
};

inline std::ostream& operator<<(std::ostream& os, const mEdge& c){
    os << "weight: "<<c.w << " index: "<< c.n;
    return os;

}

struct __attribute__ ((packed)) mNode {
    

    mEdge operator[](std::size_t i){
        return children[i]; 
    }

    mEdge getEdge(std::size_t i){
        return this->operator[](i);
    }

    inline bool operator==(const mNode& n) const noexcept{
        return v== n.v && children == n.children;
    }

    Qubit v;
    std::array<mEdge, 4> children;

};

static_assert(std::is_aggregate_v<mNode>);


template<>
struct std::hash<mNode>{
    std::size_t operator()(const mNode& n) const noexcept {
        std::size_t h = std::hash<Qubit>()(n.v);
        for(const mEdge& e: n.children){
            h = hash_combine(h, std::hash<mEdge>()(e));
        }
        return h;
    }
};

mEdge makeEdge(Worker* w, Qubit q, const std::array<mEdge, 4> c);
mEdge makeIdent(Worker* w, Qubit q);
