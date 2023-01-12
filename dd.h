#pragma once

#include "common.h"
#include <complex>

using std_complex = std::complex<float>;

struct mNode;
struct Worker;

struct mEdge {

    Qubit getVar() const;
    bool isTerminal() const;
    mNode* getNode() const;
    
    void printMatrix() const;

    inline bool operator==(const mEdge& e) const noexcept {
        return w == e.w && n == e.n;
    }

    std_complex w;
    Index n{TERMINAL};
};

static_assert(std::is_aggregate_v<mEdge>);


struct vEdge {


};


template<>
struct std::hash<std_complex>{
    std::size_t operator()(const std_complex& v) const noexcept {
        auto h1 = std::hash<float>()(v.real());
        auto h2 = std::hash<float>()(v.imag());
        return hash_combine(h1,h2);
    }
};



template<>
struct std::hash<mEdge>{
    std::size_t operator()(const mEdge& v) const noexcept {
        auto h1 = std::hash<std_complex>()(v.w);
        auto h2 = std::hash<Index>()(v.n);
        return hash_combine(h1,h2);
    }
};

inline std::ostream& operator<<(std::ostream& os, const mEdge& c){
    os << "weight: "<<c.w << " index: "<< c.n;
    return os;

}

struct mNode {
    

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

struct compare_node_ut{
    bool operator()(const mNode& lhs, const mNode& rhs)const {
        return (lhs.v == rhs.v) && (lhs.children == rhs.children);
    }

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


extern std::vector<mEdge> identityTable;



struct Job;

mEdge makeEdge(Worker* w, Qubit q, const std::array<mEdge, 4> c);
mEdge makeIdent(Worker* w, QubitCount q);
mEdge makeGate(Worker* w, GateMatrix g, QubitCount q, Qubit target, const Controls& c );


mEdge add(Worker* w, const mEdge& lhs, const mEdge& rhs);
mEdge multiply(Worker* w, const mEdge& lhs, const mEdge& rhs);


// serial addition of jobs in the vector of the rang [start,end)
mEdge addSerial(Worker* w, const std::vector<Job*> jobs, std::size_t start, std::size_t end);

mEdge mulSerial(Worker* w, const std::vector<Job*> jobs, std::size_t start, std::size_t end);


