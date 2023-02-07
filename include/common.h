#ifndef COMMON_H
#define COMMON_H

#include <utility>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <set>
#include <variant>
#include <complex>
#include <chrono>
#ifndef LINE_SIZE
#define LINE_SIZE 64
#endif

#ifndef WORKERS
#define WORKERS 8
#endif

#ifndef cas
#define cas(ptr, old, new) (__sync_bool_compare_and_swap((ptr),(old),(new)))
#endif

using duration_micro = std::chrono::duration<double, std::micro>;

const std::size_t NBUCKETS = 32768;
//const std::size_t NBUCKETS = 524288;
//const std::size_t NBUCKETS = 1024288;
const std::size_t INITIAL_ALLOCATION_SIZE = 2048;
const std::size_t GROWTH_FACTOR = 2;

const unsigned int MAX_REF = std::numeric_limits<unsigned int>::max();

using Qubit = int32_t;
using QubitCount = uint32_t;
using Index = std::size_t;

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
    c(const std::complex<float>& cf): r(cf.real()), i(cf.imag()){}
};

constexpr std::complex<float> cf_one       = {1., 0.};
constexpr std::complex<float> cf_mone      = {-1., 0.};
constexpr std::complex<float> cf_zero      = {0., 0.};
constexpr std::complex<float> cf_i         = {0., 1.};
constexpr std::complex<float> cf_mi        = {0., -1.};
constexpr std::complex<float> cf_SQRT2_2   = {SQRT2, 0.};
constexpr std::complex<float> cf_mSQRT2_2  = {-SQRT2, 0.};
constexpr std::complex<float> cf_iSQRT2_2  = {0., SQRT2};
constexpr std::complex<float> cf_miSQRT2_2 = {0., -SQRT2};
constexpr std::complex<float> cf_1plusi    = {SQRT2, SQRT2};
constexpr std::complex<float> cf_1minusi   = {SQRT2, -SQRT2};
constexpr std::complex<float> cf_1plusi_2  = {0.5, 0.5};
constexpr std::complex<float> cf_1minusi_2 = {0.5, -0.5};





    // Gate matrices
using GateMatrix = std::array<std::complex<float>, 4>;
constexpr GateMatrix Imat{cf_one, cf_zero, cf_zero, cf_one};
constexpr GateMatrix Hmat{cf_SQRT2_2, cf_SQRT2_2, cf_SQRT2_2, cf_mSQRT2_2};
constexpr GateMatrix Xmat{cf_zero, cf_one, cf_one, cf_zero};
constexpr GateMatrix Ymat{cf_zero, cf_mi, cf_i, cf_zero};
constexpr GateMatrix Zmat{cf_one, cf_zero, cf_zero, cf_mone};
constexpr GateMatrix Smat{cf_one, cf_zero, cf_zero, cf_i};
constexpr GateMatrix Sdagmat{cf_one, cf_zero, cf_zero, cf_mi};
constexpr GateMatrix Tmat{cf_one, cf_zero, cf_zero, cf_1plusi};
constexpr GateMatrix Tdagmat{cf_one, cf_zero, cf_zero, cf_1minusi};
constexpr GateMatrix SXmat{cf_1plusi_2, cf_1minusi_2, cf_1minusi_2, cf_1plusi_2};
constexpr GateMatrix SXdagmat{cf_1minusi_2, cf_1plusi_2, cf_1plusi_2, cf_1minusi_2};
constexpr GateMatrix Vmat{cf_SQRT2_2, cf_miSQRT2_2, cf_miSQRT2_2, cf_SQRT2_2};
constexpr GateMatrix Vdagmat{cf_SQRT2_2, cf_iSQRT2_2, cf_iSQRT2_2, cf_SQRT2_2};



// ----------------------------------------------------------------------------
// std::variant
// ----------------------------------------------------------------------------
template <typename T, typename>
struct get_index;

template <size_t I, typename... Ts>
struct get_index_impl {};

template <size_t I, typename T, typename... Ts>
struct get_index_impl<I, T, T, Ts...> : std::integral_constant<size_t, I>{};

template <size_t I, typename T, typename U, typename... Ts>
struct get_index_impl<I, T, U, Ts...> : get_index_impl<I+1, T, Ts...>{};

template <typename T, typename... Ts>
struct get_index<T, std::variant<Ts...>> : get_index_impl<0, T, Ts...>{};

template <typename T, typename... Ts>
constexpr auto get_index_v = get_index<T, Ts...>::value;

#endif
