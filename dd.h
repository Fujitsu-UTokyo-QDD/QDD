#pragma once

#include "common.h"
#include <complex>
#include <vector>



struct Complex{
    float r; 
    float i;

    Complex& operator+=(const Complex& rhs) noexcept{
        r += rhs.r;
        i += rhs.i;
        return *this;
    }
    Complex& operator-=(const Complex& rhs) noexcept{
        r -= rhs.r;
        i -= rhs.i;
        return *this;
    }
    Complex& operator*=(const Complex& rhs) noexcept{
        float a = r, b = i;
        float p = rhs.r, q = rhs.i;

        r = (a*p - b*q);
        i = (a*q + b*p);
        return *this;
    }
    Complex& operator/=(const Complex& rhs) noexcept{

        float mag = rhs.r*rhs.r + rhs.i*rhs.i;
        float a = r, b = i;
        float p = rhs.r, q = rhs.i;
        r = (a*p+b*q)/mag;
        i = (b*p-a*q)/mag;
        return *this;
    }

    bool operator==(const Complex& rhs) const noexcept {
        return r == rhs.r && i == rhs.i;
    }

    float real() const noexcept { return r;}
    float imag() const noexcept { return i;}
};

inline Complex operator+(Complex lhs, const Complex& rhs) noexcept{

    lhs += rhs;
    return lhs;
}
inline Complex operator-(Complex lhs, const Complex& rhs) noexcept{
    lhs -= rhs;
    return lhs;

}
inline Complex operator*(Complex lhs, const Complex& rhs) noexcept{
    lhs *= rhs;
    return lhs;

}
inline Complex operator/(Complex lhs, const Complex& rhs) noexcept{
    lhs /= rhs;
    return lhs;
}

inline std::ostream& operator<<(std::ostream& os, const Complex& c) noexcept {
    return os<<"("<<c.r<<" , "<<c.i<<")";
}

inline float norm(const Complex& c){
    return std::sqrt(c.r*c.r + c.i*c.i);
}


//using std_complex = std::complex<float>;
using std_complex = Complex;



struct mNode;
struct Worker;
struct vNode;

struct vEdge {

    Qubit getVar() const;
    bool isTerminal() const;
    vNode* getNode() const;

    void printVector() const;

    inline bool operator==(const vEdge& e) const noexcept {
        return w == e.w && n == e.n;
    }

    std_complex w;
    vNode* n;


};

struct vNode {

    vEdge operator[](std::size_t i){
        return children[i]; 
    }

    vEdge getEdge(std::size_t i){
        return this->operator[](i);
    }

    inline bool operator==(const vNode& n) const noexcept{
        return v== n.v && children == n.children;
    }

    static vNode       terminalNode;
    constexpr static vNode*   terminal{&terminalNode};

    Qubit v;
    std::array<vEdge, 2> children;
};

struct mEdge {

    static mEdge one;
    static mEdge zero;

    Qubit getVar() const;
    bool isTerminal() const;
    mNode* getNode() const;
    
    vEdge get_column(std::size_t col) const;
    void printMatrix() const;



    inline bool operator==(const mEdge& e) const noexcept {
        return w == e.w && n == e.n;
    }

    std_complex w;
    mNode*    n;
};

static_assert(std::is_aggregate_v<mEdge>);



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
    std::size_t operator()(const mEdge& e) const noexcept {
        auto h1 = std::hash<std_complex>()(e.w);
        auto h2 = std::hash<std::size_t>()(reinterpret_cast<std::size_t>(e.n));
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

    static mNode       terminalNode;
    constexpr static mNode*   terminal{&terminalNode};




    Qubit v;
    std::array<mEdge, 4> children;
    mNode*   next{nullptr};

};





static_assert(std::is_aggregate_v<mNode>);


template<>
struct std::hash<mNode>{
    std::size_t operator()(const mNode& n) const noexcept {
        std::size_t h = 0; // don't hash the Qubit since each qubit has its own uniqueTable
        for(const mEdge& e: n.children){
            h = hash_combine(h, std::hash<mEdge>()(e));
        }
        return h;
    }
};


extern std::vector<mEdge> identityTable;



struct Job;

mEdge makeEdge(Worker* w, Qubit q, const std::array<mEdge, 4>& c);
mEdge makeIdent(Worker* w, QubitCount q);
mEdge makeGate(Worker* w, GateMatrix g, QubitCount q, Qubit target, const Controls& c );


mEdge add(Worker* w, const mEdge& lhs, const mEdge& rhs);
mEdge multiply(Worker* w, const mEdge& lhs, const mEdge& rhs);
mEdge kronecker(Worker* w, const mEdge& lhs, const mEdge& rhs);


// serial addition of jobs in the vector of the rang [start,end)
mEdge addSerial(Worker* w,  const std::vector<Job*>& jobs, std::size_t start, std::size_t end);

mEdge mulSerial(Worker* w,  const std::vector<Job*>& jobs, std::size_t start, std::size_t end);


