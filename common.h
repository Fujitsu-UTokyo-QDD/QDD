#ifndef COMMON_H
#define COMMON_H

#include <utility>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <set>
#ifndef LINE_SIZE
#define LINE_SIZE 64
#endif

#ifndef WORKERS
#define WORKERS 8
#endif

#ifndef cas
#define cas(ptr, old, new) (__sync_bool_compare_and_swap((ptr),(old),(new)))
#endif


using Qubit = int32_t;
using QubitCount = uint32_t;
using Index = std::size_t;
using Controls = std::set<Qubit>;

inline Index operator""_idx(unsigned long long int idx){
    return static_cast<Index>(idx);
}
inline Qubit operator""_q(unsigned long long int q){
    return static_cast<Qubit>(q);
}

const Index TERMINAL = 0_idx;


using double_pair = std::pair<double, double>;

const uint64_t INDEX_MASK  = 0x00000000FFFFFFFFLL; // the last 4 bytes
const uint64_t HASH_MASK   = 0xFFFFFFFF00000000LL; // the first 4 bytes
const uint64_t HASH_SHIFT  = 32;



constexpr std::size_t hash_combine(std::size_t lhs, std::size_t rhs) {
    lhs ^= rhs + 0x9e3779b97f4a7c15ULL + (lhs << 6) + (lhs >> 2);
    return lhs;
}

constexpr std::size_t murmur_hash(std::size_t k){
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

const std::size_t REGION_SIZE = 512;
const uint64_t UT_INIT_SZ =   8 * REGION_SIZE;
const uint64_t UT_MAX_SZ  =   16 * REGION_SIZE;

constexpr double SQRT2 = 0.707106781186547524400844362104849039284835937688474036588L;
constexpr double mSQRT2 = -0.707106781186547524400844362104849039284835937688474036588L;

static_assert(sizeof(double) == sizeof(unsigned long long));


struct c{
    float r;
    float i;
};



constexpr c complex_one       = {1., 0.};
constexpr c complex_mone      = {-1., 0.};
constexpr c complex_zero      = {0., 0.};
constexpr c complex_i         = {0., 1.};
constexpr c complex_mi        = {0., -1.};
constexpr c complex_SQRT2_2   = {SQRT2, 0.};
constexpr c complex_mSQRT2_2  = {-SQRT2, 0.};
constexpr c complex_iSQRT2_2  = {0., SQRT2};
constexpr c complex_miSQRT2_2 = {0., -SQRT2};
constexpr c complex_1plusi    = {SQRT2, SQRT2};
constexpr c complex_1minusi   = {SQRT2, -SQRT2};
constexpr c complex_1plusi_2  = {0.5, 0.5};
constexpr c complex_1minusi_2 = {0.5, -0.5};

    // Gate matrices
using GateMatrix = std::array<c, 4>;
constexpr GateMatrix Imat{complex_one, complex_zero, complex_zero, complex_one};
constexpr GateMatrix Hmat{complex_SQRT2_2, complex_SQRT2_2, complex_SQRT2_2, complex_mSQRT2_2};
constexpr GateMatrix Xmat{complex_zero, complex_one, complex_one, complex_zero};
constexpr GateMatrix Ymat{complex_zero, complex_mi, complex_i, complex_zero};
constexpr GateMatrix Zmat{complex_one, complex_zero, complex_zero, complex_mone};
constexpr GateMatrix Smat{complex_one, complex_zero, complex_zero, complex_i};
constexpr GateMatrix Sdagmat{complex_one, complex_zero, complex_zero, complex_mi};
constexpr GateMatrix Tmat{complex_one, complex_zero, complex_zero, complex_1plusi};
constexpr GateMatrix Tdagmat{complex_one, complex_zero, complex_zero, complex_1minusi};
constexpr GateMatrix SXmat{complex_1plusi_2, complex_1minusi_2, complex_1minusi_2, complex_1plusi_2};
constexpr GateMatrix SXdagmat{complex_1minusi_2, complex_1plusi_2, complex_1plusi_2, complex_1minusi_2};
constexpr GateMatrix Vmat{complex_SQRT2_2, complex_miSQRT2_2, complex_miSQRT2_2, complex_SQRT2_2};
constexpr GateMatrix Vdagmat{complex_SQRT2_2, complex_iSQRT2_2, complex_iSQRT2_2, complex_SQRT2_2};

#endif
