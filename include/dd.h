#pragma once

#include "Eigen/Dense"
#include "common.h"
#include <atomic>
#include <complex>
#include <random>
#include <vector>
using Eigen::MatrixXcf;
using Eigen::VectorXcf;


#ifdef isMPI
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>
namespace bmpi = boost::mpi;
#endif

struct Complex {
#ifdef isMPI
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar &r;
        ar &i;
    }
#endif
    double r{0.0};
    double i{0.0};

    static inline const double TOLERANCE =
        std::numeric_limits<double>::epsilon() * 1024;

    Complex() = default;
    Complex(const std::complex<float> &c) : r(c.real()), i(c.imag()) {}
    Complex(double rr, double ii) : r(rr), i(ii) {}

    Complex &operator+=(const Complex &rhs) noexcept {
        r += rhs.r;
        i += rhs.i;
        return *this;
    }
    Complex &operator-=(const Complex &rhs) noexcept {
        r -= rhs.r;
        i -= rhs.i;
        return *this;
    }
    Complex &operator*=(const Complex &rhs) noexcept {
        double a = r, b = i;
        double p = rhs.r, q = rhs.i;

        r = (a * p - b * q);
        i = (a * q + b * p);
        return *this;
    }
    Complex &operator/=(const Complex &rhs) noexcept {

        double mag = rhs.r * rhs.r + rhs.i * rhs.i;
        double a = r, b = i;
        double p = rhs.r, q = rhs.i;
        r = (a * p + b * q) / mag;
        i = (b * p - a * q) / mag;
        return *this;
    }

    double mag2() const { return r * r + i * i; }

    bool isZero() const noexcept { return r == 0.0 && i == 0.0; }

    bool isApproximatelyZero() const noexcept {
        return std::abs(r) <= TOLERANCE && std::abs(i) <= TOLERANCE;
    }

    bool isApproximatelyEqual(const Complex &rhs) const noexcept {
        if (r == rhs.r && i == rhs.i)
            return true;
        if (std::abs(r - rhs.r) <= TOLERANCE &&
            std::abs(i - rhs.i) <= TOLERANCE)
            return true;

        return false;
    }

    bool isApproximatelyOne() const noexcept {
        return isApproximatelyEqual({1.0, 0.0});
    }

    bool isOne() const noexcept { return r == 1.0 && i == 0.0; }

    bool operator==(const Complex &rhs) const noexcept {
        return r == rhs.r && i == rhs.i;
    }

    bool operator==(const std::complex<double> rhs) const noexcept {
        return r == rhs.real() && i == rhs.imag();
    }

    double real() const noexcept { return r; }
    double imag() const noexcept { return i; }
};

inline Complex operator+(Complex lhs, const Complex &rhs) noexcept {

    lhs += rhs;
    return lhs;
}
inline Complex operator-(Complex lhs, const Complex &rhs) noexcept {
    lhs -= rhs;
    return lhs;
}
inline Complex operator*(Complex lhs, const Complex &rhs) noexcept {
    lhs *= rhs;
    return lhs;
}
inline Complex operator/(Complex lhs, const Complex &rhs) noexcept {
    lhs /= rhs;
    return lhs;
}

inline std::ostream &operator<<(std::ostream &os, const Complex &c) noexcept {
    return os << "(" << c.r << " , " << c.i << ")";
}

inline double norm(const Complex &c) {
    return std::sqrt(c.r * c.r + c.i * c.i);
}

// using std_complex = std::complex<double>;
using std_complex = Complex;

struct mNode;
struct vNode;

struct vEdge {
#ifdef isMPI
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar &w;
        ar &n;
    }
#endif

    static vEdge one;
    static vEdge zero;

    Qubit getVar() const;
    bool isTerminal() const;
    vNode *getNode() const { return n; };

    void printVector() const;
#ifdef isMPI
    void printVectorMPI(bmpi::communicator &world) const;
#endif
    void printVector_sparse() const;
    std_complex *getVector(std::size_t *dim) const;

    inline bool operator==(const vEdge &e) const noexcept {
        return w.isApproximatelyEqual(e.w) && n == e.n;
    }

    std_complex w;
    vNode *n{nullptr};
};
inline void swap(vEdge &lhs, vEdge &rhs) {
    using std::swap;
    swap(lhs.w, rhs.w);
    swap(lhs.n, rhs.n);
}

struct vNode {
#ifdef isMPI
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar &v;
        ar &children;
        ar &next;
    }
#endif

    vNode() = default;
    vNode(const vNode &vv) : v(vv.v), children(vv.children) {}
    vNode(Qubit q, const std::array<vEdge, 2> &c, vNode *n)
        : v(q), children(c), next(n) {}

    vEdge &operator[](std::size_t i) { return children[i]; }

    vEdge &getEdge(std::size_t i) { return this->operator[](i); }

    inline bool operator==(const vNode &n) const noexcept {
        return v == n.v && children == n.children;
    }

    static vNode terminalNode;
    constexpr static vNode *terminal{&terminalNode};

    Qubit v;
    std::array<vEdge, 2> children;
    vNode *next{nullptr};

};

struct mEdge {
#ifdef isMPI
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar &w;
        ar &n;
    }
#endif

    static mEdge one;
    static mEdge zero;

    Qubit getVar() const;
    bool isTerminal() const;
    mNode *getNode() const { return n; };

    void printMatrix() const;

    std_complex **getMatrix(std::size_t *dim) const;
    MatrixXcf getEigenMatrix();

    inline bool operator==(const mEdge &e) const noexcept {

        // because nullptr == nullptr
        return w.isApproximatelyEqual(e.w) && n == e.n;
    }

    inline bool operator!=(const mEdge &e) const noexcept {
        return !(this->operator==(e));
    }

    std_complex w;

    mNode *n{nullptr};
};

inline void swap(mEdge &lhs, mEdge &rhs) {
    using std::swap;
    swap(lhs.w, rhs.w);
    swap(lhs.n, rhs.n);
}

template <> struct std::hash<std_complex> {
    std::size_t operator()(const std_complex &v) const noexcept {
        auto h1 = std::hash<double>()(v.real());
        auto h2 = std::hash<double>()(v.imag());
        return hash_combine(h1, h2);
    }
};

template <> struct std::hash<mEdge> {
    std::size_t operator()(const mEdge &e) const noexcept {
        auto h1 = std::hash<std_complex>()(e.w);
        std::size_t h2;
        h2 = std::hash<std::size_t>()(reinterpret_cast<std::size_t>(e.n));

        return hash_combine(h1, h2);
    }
};

template <> struct std::hash<vEdge> {
    std::size_t operator()(const vEdge &e) const noexcept {
        auto h1 = std::hash<std_complex>()(e.w);
        auto h2 = std::hash<std::size_t>()(reinterpret_cast<std::size_t>(e.n));
        return hash_combine(h1, h2);
    }
};

inline std::ostream &operator<<(std::ostream &os, const mEdge &c) {
    os << "weight: " << c.w << " index: " << c.n;
    return os;
}

struct mNode {
#ifdef isMPI
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar &v;
        ar &children;
        ar &next;
    }
#endif

    mNode() = default;
    mNode(Qubit q, const std::array<mEdge, 4> &c, mNode *n)
        : v(q), children(c), next(n){}
    mNode(const mNode &vv) : v(vv.v), children(vv.children) {}
    mEdge &operator[](std::size_t i) { return children[i]; }

    mEdge &getEdge(std::size_t i) { return this->operator[](i); }

    inline bool operator==(const mNode &n) const noexcept {
        return v == n.v && children == n.children;
    }

    static mNode terminalNode;
    constexpr static mNode *terminal{&terminalNode};

    Qubit v;
    std::array<mEdge, 4> children;
    mNode *next{nullptr};

};

template <> struct std::hash<mNode> {
    std::size_t operator()(const mNode &n) const noexcept {
        std::size_t h =
            0; // don't hash the Qubit since each qubit has its own uniqueTable
        for (const mEdge &e : n.children) {
            h = hash_combine(h, std::hash<mEdge>()(e));
        }
        return h;
    }
};

template <> struct std::hash<vNode> {
    std::size_t operator()(const vNode &n) const noexcept {
        std::size_t h =
            0; // don't hash the Qubit since each qubit has its own uniqueTable
        for (const vEdge &e : n.children) {
            h = hash_combine(h, std::hash<vEdge>()(e));
        }
        return h;
    }
};

struct Control {
    enum class Type : bool { pos = true, neg = false };
    Qubit qubit;
    Type type = Type::pos;
};
inline bool operator<(const Control &lhs, const Control &rhs) {
    return lhs.qubit < rhs.qubit ||
           (lhs.qubit == rhs.qubit && lhs.type < rhs.type);
}
inline bool operator==(const Control &lhs, const Control &rhs) {
    return lhs.qubit == rhs.qubit && lhs.type == rhs.type;
}
inline bool operator!=(const Control &lhs, const Control &rhs) {
    return !(lhs == rhs);
}

struct ControlComparator {
    inline bool operator()(const Control &lhs, const Control &rhs) const {
        return lhs < rhs;
    }

    inline bool operator()(Qubit lhs, const Control &rhs) const {
        return lhs < rhs.qubit;
    }

    inline bool operator()(const Control &lhs, Qubit rhs) const {
        return lhs.qubit < rhs;
    }
};
struct vContent {
#ifdef isMPI
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar &v;
        ar &w;
        ar &index;
    }
#endif

    Qubit v;
    std::array<std_complex, 2> w;
    std::array<int, 2> index;

    vContent(Qubit qpos, std_complex w1, std_complex w2, int i1, int i2)
        : v(qpos), w({w1, w2}), index({i1, i2}){};
    vContent()
        : v(-1), w({cf_zero, cf_zero}), index({-1,-1}){};
};

using Controls = std::set<Control, ControlComparator>;

extern std::vector<mEdge> identityTable;

mEdge makeMEdge(Qubit q, const std::array<mEdge, 4> &c);
mEdge makeIdent(Qubit q);
mEdge makeGate(QubitCount q, GateMatrix g, Qubit target, const Controls &c);
mEdge makeGate(QubitCount q, GateMatrix g, Qubit target);
mEdge makeSwap(QubitCount q, Qubit target0, Qubit target1);

mEdge getMPIGate(mEdge root, int row, int col, int world_size);

mEdge mm_add(const mEdge &lhs, const mEdge &rhs);
mEdge mm_multiply(const mEdge &lhs, const mEdge &rhs);

mEdge mm_kronecker(const mEdge &lhs, const mEdge &rhs);

vEdge vv_add(const vEdge &lhs, const vEdge &rhs);
vEdge vv_kronecker(const vEdge &lhs, const vEdge &rhs);

vEdge mv_multiply(mEdge lhs, vEdge rhs);
#ifdef isMPI
vEdge mv_multiply_MPI(mEdge lhs, vEdge rhs, bmpi::communicator &world);
#endif

vEdge makeVEdge(Qubit q, const std::array<vEdge, 2> &c);
vEdge makeZeroState(QubitCount q);
vEdge makeOneState(QubitCount q);
#ifdef isMPI
vEdge makeZeroStateMPI(QubitCount q, bmpi::communicator &world);
vEdge makeOneStateMPI(QubitCount q, bmpi::communicator &world);
#endif

std::string measureAll(vEdge &rootEdge, const bool collapse,
                       std::mt19937_64 &mt, double epsilon = 0.001);
char measureOneCollapsing(vEdge &rootEdge, const Qubit index,
                          std::mt19937_64 &mt, double epsilon = 0.001);

mEdge RX(QubitCount qnum, int target, float angle);

mEdge RY(QubitCount qnum, int target, float angle);

mEdge RZ(QubitCount qnum, int target, float angle);

mEdge CX(QubitCount qnum, int target, int control);

GateMatrix rx(float angle);
GateMatrix ry(float angle);
GateMatrix rz(float angle);
GateMatrix u3(float theta, float phi, float lambda);
GateMatrix u1(float lambda);
GateMatrix u2(float phi, float lambda);
GateMatrix u(float theta, float phi, float lambda);
GateMatrix p(float angle);
GateMatrix r(float theta, float phi);

std::string genDot(vEdge &rootEdge);
std::string genDot(mEdge &rootEdge);

#ifdef isMPI
vEdge receive_dd(boost::mpi::communicator &world, int source_node_id, bool isBlocking = true);
void send_dd(boost::mpi::communicator &world, vEdge e, int dest_node_id, bool isBlocking = true);
#endif

void dump(boost::mpi::communicator &world, vEdge e, int cycle);