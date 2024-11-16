#include "dd.h"

#include <math.h>

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>

#include "cache.hpp"
#include "common.h"
#include "table.hpp"

#ifdef isMPI
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/utility.hpp>
#endif

#ifdef isMT
#include <oneapi/tbb/enumerable_thread_specific.h>

#include <boost/fiber/future/future.hpp>
#include <boost/fiber/future/packaged_task.hpp>
#endif

#define SUBTASK_THRESHOLD 5

mNodeTable mUnique(NQUBITS);
vNodeTable vUnique(NQUBITS);
std::unordered_map<vNode *, double> probs;

std::vector<mEdge> identityTable(NQUBITS);

mNode mNode::terminalNode = mNode(-1, {}, nullptr, nullptr);
vNode vNode::terminalNode = vNode(-1, {}, nullptr, nullptr);

mEdge mEdge::one{.w = {1.0, 0.0}, .n = mNode::terminal};
mEdge mEdge::zero{.w = {0.0, 0.0}, .n = mNode::terminal};
vEdge vEdge::one{.w = {1.0, 0.0}, .n = vNode::terminal};
vEdge vEdge::zero{.w = {0.0, 0.0}, .n = vNode::terminal};

#if defined(isMT) && !defined(CACHE_GLOBAL)
oneapi::tbb::enumerable_thread_specific<AddCache> _aCaches(NQUBITS);
oneapi::tbb::enumerable_thread_specific<MulCache> _mCaches(NQUBITS);
#else
AddCache _aCache(NQUBITS);
MulCache _mCache(NQUBITS);
#endif

static int LIMIT = 10000;
const int MINUS = 6;

static mEdge normalizeM(const mEdge &e) {
    // check for all zero weights
    if (std::all_of(e.n->children.begin(), e.n->children.end(),
                    [](const mEdge &e) { return norm(e.w) == 0.0; })) {
        mUnique.returnNode(e.n);
        return mEdge::zero;
    }

    auto result = std::max_element(e.n->children.begin(), e.n->children.end(),
                                   [](const mEdge &lhs, const mEdge &rhs) {
                                       return norm2(lhs.w) < norm2(rhs.w);
                                   });

    std_complex max_weight = result->w;
    const std::size_t idx = std::distance(e.n->children.begin(), result);

    // parents weight
    std_complex new_weight = max_weight * e.w;
    if (new_weight.isApproximatelyZero()) {
        mUnique.returnNode(e.n);
        return mEdge::zero;
    } else if (new_weight.isApproximatelyOne()) {
        new_weight = {1.0, 0.0};
    }

    for (int i = 0; i < 4; i++) {
        std_complex r = e.n->children[i].w / max_weight;
        if (r.isApproximatelyZero()) {
            e.n->children[i] = mEdge::zero;
        } else {
            if (r.isApproximatelyOne()) {
                r = {1.0, 0.0};
            }
            e.n->children[i].w = r;
        }
    }

    mNode *n = mUnique.lookup(e.n);
    assert(n->v >= -1);

    return {.w = new_weight, .n = n};
}

static vEdge normalizeV(const vEdge &e) {
    // check for all zero weights
    if (std::all_of(e.n->children.begin(), e.n->children.end(),
                    [](const vEdge &e) { return e.w.isApproximatelyZero(); })) {
        vUnique.returnNode(e.n);
        return vEdge::zero;
    }

    auto result = std::max_element(e.n->children.begin(), e.n->children.end(),
                                   [](const vEdge &lhs, const vEdge &rhs) {
                                       return norm2(lhs.w) < norm2(rhs.w);
                                   });
    std_complex max_weight = result->w;

    // parents weight
    std_complex new_weight = max_weight * e.w;
    if (new_weight.isApproximatelyZero()) {
        vUnique.returnNode(e.n);
        return vEdge::zero;
    } else if (new_weight.isApproximatelyOne()) {
        new_weight = {1.0, 0.0};
    }

    // child weight (larger one)
    size_t max_idx = (result == &(e.n->children[1])) ? 1 : 0;
    e.n->children[max_idx].w = {1.0, 0.0};

    // child weight (smaller one)
    size_t min_idx = (max_idx == 1) ? 0 : 1;
    e.n->children[min_idx].w = e.n->children[min_idx].w / max_weight;
    if (e.n->children[min_idx].w.isApproximatelyOne()) {
        e.n->children[min_idx].w = {1.0, 0.0};
    } else if (e.n->children[min_idx].w.isApproximatelyZero()) {
        e.n->children[min_idx] = vEdge::zero;
    }

    // making new node
    vNode *n = vUnique.lookup(e.n);
    return {new_weight, n};
}

mEdge makeMEdge(Qubit q, const std::array<mEdge, 4> &c) {
    // Identity Stripping
    if( c[0].n==c[3].n && c[1].w.isApproximatelyZero() && c[2].w.isApproximatelyZero() && c[0].w.isApproximatelyEqual(c[3].w)){
        return {c[0].w, c[0].n};
    }

    mNode *node = mUnique.getNode();
    node->v = q;
    node->children = c;

    mEdge e = normalizeM({.w = {1.0, 0.0}, .n = node});

    return e;
}

vEdge makeVEdge(Qubit q, const std::array<vEdge, 2> &c) {
    vNode *node = vUnique.getNode();
    node->v = q;
    node->children = c;

    for (int i = 0; i < 2; i++) {
        assert(&node->children[i] != &vEdge::one &&
               (&node->children[i] != &vEdge::zero));
    }

    vEdge e = normalizeV({{1.0, 0.0}, node});

    assert(e.getVar() == q || e.isTerminal());

    return e;
}

Qubit mEdge::getVar() const { return n->v; }

Qubit vEdge::getVar() const { return n->v; }

bool mEdge::isTerminal() const { return n == mNode::terminal; }

bool vEdge::isTerminal() const { return n == vNode::terminal; }

static void fillMatrix(const mEdge &edge, size_t row, size_t col,
                       const std_complex &w, uint64_t dim, std_complex **m) {
    if ( (1 << (edge.getVar()+1)) == dim) {
        std_complex wp = edge.w * w;

        if (edge.isTerminal()) {
            for (auto i = row; i < row + dim; i++) {
                for (auto j = col; j < col + dim; j++) {
                    m[i][j] = wp;
                }
            }
            return;
        }

        mNode *node = edge.getNode();
        fillMatrix(node->getEdge(0), row, col, wp, dim / 2, m);
        fillMatrix(node->getEdge(1), row, col + dim / 2, wp, dim / 2, m);
        fillMatrix(node->getEdge(2), row + dim / 2, col, wp, dim / 2, m);
        fillMatrix(node->getEdge(3), row + dim / 2, col + dim / 2, wp, dim / 2, m);
    }else{
        std_complex wp = w;

        fillMatrix(edge, row, col, wp, dim / 2, m);
        fillMatrix(edge, row, col + dim / 2, {0.0,0.0}, dim / 2, m);
        fillMatrix(edge, row + dim / 2, col, {0.0,0.0}, dim / 2, m);
        fillMatrix(edge, row + dim / 2, col + dim / 2, wp, dim / 2, m);
    }
}

void mEdge::printMatrix(Qubit nQubits) const {
    if(nQubits==-1){
        nQubits = this->getVar() + 1;
    }
    std::size_t dim = 1 << nQubits;

    std_complex **matrix = new std_complex *[dim];
    for (std::size_t i = 0; i < dim; i++) matrix[i] = new std_complex[dim];

    fillMatrix(*this, 0, 0, {1.0, 0.0}, dim, matrix);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < dim; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

std_complex **mEdge::getMatrix(std::size_t *dim, Qubit nQubits) const {
    if(nQubits==-1){
        nQubits = this->getVar() + 1;
    }
    std::size_t d = 1 << nQubits;

    std_complex **matrix = new std_complex *[d];
    for (std::size_t i = 0; i < d; i++) matrix[i] = new std_complex[d];

    fillMatrix(*this, 0, 0, {1.0, 0.0}, d, matrix);
    if (dim != nullptr) *dim = d;
    return matrix;
}

struct MatrixGuard {
    MatrixGuard(std_complex **m, std::size_t dim) : _m(m), _dim(dim) {}
    ~MatrixGuard() {
        for (size_t i = 0; i < _dim; i++) {
            delete[] _m[i];
        }
        delete[] _m;
    }

    std_complex **_m;
    const std::size_t _dim;
};

MatrixXcf mEdge::getEigenMatrix(Qubit nQubit) {
    std::size_t dim;
    auto m = getMatrix(&dim, nQubit);
    MatrixGuard g(m, dim);
    MatrixXcf M(dim, dim);

    for (auto i = 0; i < dim; i++) {
        for (auto j = 0; j < dim; j++) {
            M(i, j) = std::complex<double>{m[i][j].r, m[i][j].i};
        }
    }
    return M;
}

mEdge makeIdent(Qubit q) {
    return mEdge::one;
}

vEdge makeZeroState(QubitCount q) {
    vEdge e = makeVEdge(0, {vEdge::one, vEdge::zero});
    for (Qubit i = 1; i < q; i++) {
        e = makeVEdge(i, {{e, vEdge::zero}});
    }
    return e;
}

#ifdef isMPI
vEdge makeZeroStateMPI(QubitCount q, bmpi::communicator &world) {
    if (world.rank() == 0) {
        int shift = std::log2(world.size());
        assert(1 << shift == world.size());
        return makeZeroState(q - shift);
    } else {
        return vEdge::zero;
    }
}
#endif

vEdge makeOneState(QubitCount q) {
    vEdge e = makeVEdge(0, {vEdge::zero, vEdge::one});
    for (Qubit i = 1; i < q; i++) {
        e = makeVEdge(i, {{vEdge::zero, e}});
    }
    return e;
}

#ifdef isMPI
vEdge makeOneStateMPI(QubitCount q, bmpi::communicator &world) {
    if (world.rank() == world.size() - 1) {
        int shift = std::log2(world.size());
        assert(1 << shift == world.size());
        return makeOneState(q - shift);
    } else {
        return vEdge::zero;
    }
}
#endif

mEdge makeGate(QubitCount q, GateMatrix g, Qubit target) {
    return makeGate(q, g, target, {});
}

mEdge makeGate(QubitCount q, GateMatrix g, Qubit target, const Controls &c) {
    std::array<mEdge, 4> edges;

    for (auto i = 0; i < 4; i++)
        edges[i] = mEdge{{g[i].real(), g[i].imag()}, mNode::terminal};

    auto it = c.begin();

    Qubit z = 0;
    for (; z < target; z++) {
        for (int b1 = 0; b1 < 2; b1++) {
            for (int b0 = 0; b0 < 2; b0++) {
                std::size_t i = (b1 << 1) | b0;
                if (it != c.end() && it->qubit == z) {
                    if (it->type == Control::Type::neg)
                        edges[i] = makeMEdge(
                            z, {edges[i], mEdge::zero, mEdge::zero,
                                (b1 == b0) ? makeIdent(z - 1) : mEdge::zero});
                    else
                        edges[i] = makeMEdge(
                            z, {(b1 == b0) ? makeIdent(z - 1) : mEdge::zero,
                                mEdge::zero, mEdge::zero, edges[i]});

                } else {
                    edges[i] = makeMEdge(
                        z, {edges[i], mEdge::zero, mEdge::zero, edges[i]});
                }
            }
        }

        if (it != c.end() && it->qubit == z) ++it;
    }

    auto e = makeMEdge(z, edges);

    for (z = z + 1; z < q; z++) {
        if (it != c.end() && it->qubit == z) {
            if (it->type == Control::Type::neg)
                e = makeMEdge(z,
                              {e, mEdge::zero, mEdge::zero, makeIdent(z - 1)});
            else
                e = makeMEdge(z,
                              {makeIdent(z - 1), mEdge::zero, mEdge::zero, e});
            ++it;
        } else {
            e = makeMEdge(z, {e, mEdge::zero, mEdge::zero, e});
        }
    }

    return e;
}

mEdge makeTwoQubitGate(QubitCount q, TwoQubitGateMatrix g, Qubit target0,
                       Qubit target1) {
    return makeTwoQubitGate(q, g, target0, target1, {});
};

mEdge makeTwoQubitGate(QubitCount q, TwoQubitGateMatrix g, Qubit target0,
                       Qubit target1, const Controls &controls) {
    std::array<std::array<mEdge, 4>, 4> edges_4x4;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            edges_4x4[i][j] =
                mEdge{{g[i][j].real(), g[i][j].imag()}, mNode::terminal};
        }
    }

    auto it = controls.begin();
    const auto smaller_target = std::min(target0, target1);
    Qubit z = 0;

    for (; z < smaller_target; z++) {
        for (int b1 = 0; b1 < 4; b1++) {
            for (int b0 = 0; b0 < 4; b0++) {
                std::size_t i = (b1 << 2) | b0;
                if (it != controls.end() && it->qubit == z) {
                    if (it->type == Control::Type::neg) {
                        edges_4x4[b1][b0] = makeMEdge(
                            z, {edges_4x4[b1][b0], mEdge::zero, mEdge::zero,
                                (b1 == b0) ? makeIdent(z - 1) : mEdge::zero});
                    } else {
                        edges_4x4[b1][b0] = makeMEdge(
                            z, {(b1 == b0) ? makeIdent(z - 1) : mEdge::zero,
                                mEdge::zero, mEdge::zero, edges_4x4[b1][b0]});
                    }
                } else {
                    edges_4x4[b1][b0] =
                        makeMEdge(z, {edges_4x4[b1][b0], mEdge::zero,
                                      mEdge::zero, edges_4x4[b1][b0]});
                }
            }
        }
        if (it != controls.end() && it->qubit == z) {
            ++it;
        }
    }

    std::array<mEdge, 4> edges;
    for (int b1 = 0; b1 < 2; b1++) {
        for (int b0 = 0; b0 < 2; b0++) {
            std::array<mEdge, 4> edges_4x4_sub;
            if (target0 > target1) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        edges_4x4_sub[i * 2 + j] =
                            edges_4x4[b1 * 2 + i][b0 * 2 + j];
                    }
                }
            } else {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        edges_4x4_sub[i * 2 + j] =
                            edges_4x4[i * 2 + b1][j * 2 + b0];
                    }
                }
            }
            edges[b1 * 2 + b0] = makeMEdge(smaller_target, edges_4x4_sub);
        }
    }

    const auto larger_target = std::max(target0, target1);
    for (z = z + 1; z < larger_target; z++) {
        for (int b1 = 0; b1 < 2; b1++) {
            for (int b0 = 0; b0 < 2; b0++) {
                std::size_t i = (b1 << 1) | b0;
                if (it != controls.end() && it->qubit == z) {
                    if (it->type == Control::Type::neg)
                        edges[i] = makeMEdge(
                            z, {edges[i], mEdge::zero, mEdge::zero,
                                (b1 == b0) ? makeIdent(z - 1) : mEdge::zero});
                    else
                        edges[i] = makeMEdge(
                            z, {(b1 == b0) ? makeIdent(z - 1) : mEdge::zero,
                                mEdge::zero, mEdge::zero, edges[i]});

                } else {
                    edges[i] = makeMEdge(
                        z, {edges[i], mEdge::zero, mEdge::zero, edges[i]});
                }
            }
        }
        if (it != controls.end() && it->qubit == z) ++it;
    }

    auto e = makeMEdge(z, edges);

    for (z = z + 1; z < q; z++) {
        if (it != controls.end() && it->qubit == z) {
            if (it->type == Control::Type::neg) {
                e = makeMEdge(z,
                              {e, mEdge::zero, mEdge::zero, makeIdent(z - 1)});
            } else {
                e = makeMEdge(z,
                              {makeIdent(z - 1), mEdge::zero, mEdge::zero, e});
            }
            ++it;
        } else {
            e = makeMEdge(z, {e, mEdge::zero, mEdge::zero, e});
        }
    }

    return e;
}

mEdge makeLargeGate(QubitCount q, ComplexMatrix &g) {
    std::vector<Qubit> targets;
    for (int i = 0; i < int(log2(g.size())); i++) {
        targets.push_back(i);
    }
    return makeLargeGate(q, g, targets);
}

mEdge makeLargeGate(QubitCount q, ComplexMatrix &g,
                    const std::vector<Qubit> &targets) {
    if (g.empty()) {
        return mEdge::one;
    }

    const auto &rows = g.size();
    if ((rows & rows - 1) != 0) {
        throw std::invalid_argument("The number of rows must be a power of 2.");
    }

    const auto &cols = g[0].size();
    if (rows != cols) {
        throw std::invalid_argument("The matrix must be square.");
    }

    if (rows == 1) {
        return mEdge{{g[0][0].real(), g[0][0].imag()}, mNode::terminal};
    }

    if (int(log2(rows)) != targets.size()) {
        throw std::invalid_argument(
            "The number of targets is incompatible with the matrix size.");
    }

    // make conversion map of indices
    std::vector<size_t> target_compresed;
    for (int i = 0; i < targets.size(); i++) {
        target_compresed.push_back(i);
    }
    std::sort(
        target_compresed.begin(), target_compresed.end(),
        [&targets](size_t a, size_t b) { return targets[a] < targets[b]; });
    std::vector<std::size_t> conversionMap;
    conversionMap.push_back(0);
    for (int i = 0; i < rows; i++) {
        size_t addition = 1 << target_compresed[i];
        int rep = conversionMap.size();
        for (int j = 0; j < rep; j++) {
            conversionMap.push_back(conversionMap[j] + addition);
        }
    }
    // swap rows and columns
    std::vector<std::vector<std::complex<double>>> g_tmp;
    for (int i = 0; i < rows; i++) {
        g_tmp.push_back(g[conversionMap[i]]);
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            g[i][j] = g_tmp[i][conversionMap[j]];
        }
    }

    Qubit target_min = *std::min_element(targets.begin(), targets.end());
    std::vector<bool> skip(q, true);
    for (auto target : targets) {
        skip[target] = false;
    }

    return makeLargeGate(g, q - 1, 0, rows, 0, cols, skip, target_min);
}
mEdge makeLargeGate(ComplexMatrix &g, const Qubit level,
                    const std::size_t rowStart, const std::size_t rowEnd,
                    const std::size_t colStart, const std::size_t colEnd,
                    const std::vector<bool> skip, const Qubit target_min) {
    if (level == target_min) {
        assert(rowEnd - rowStart == 2 && colEnd - colStart == 2);
        return makeMEdge(
            target_min,
            {mEdge{{g[rowStart][colStart].real(), g[rowStart][colStart].imag()},
                   mNode::terminal},
             mEdge{{g[rowStart][colStart + 1].real(),
                    g[rowStart][colStart + 1].imag()},
                   mNode::terminal},
             mEdge{{g[rowStart + 1][colStart].real(),
                    g[rowStart + 1][colStart].imag()},
                   mNode::terminal},
             mEdge{{g[rowStart + 1][colStart + 1].real(),
                    g[rowStart + 1][colStart + 1].imag()},
                   mNode::terminal}});
    }

    const auto l = level - 1;

    if (skip[level]) {
        auto e = makeLargeGate(g, l, rowStart, rowEnd, colStart, colEnd, skip,
                               target_min);
        return makeMEdge(level, {e, mEdge::zero, mEdge::zero, e});
    }

    const auto rowMid = (rowStart + rowEnd) / 2;
    const auto colMid = (colStart + colEnd) / 2;

    return makeMEdge(level, {makeLargeGate(g, l, rowStart, rowMid, colStart,
                                           colMid, skip, target_min),
                             makeLargeGate(g, l, rowStart, rowMid, colMid,
                                           colEnd, skip, target_min),
                             makeLargeGate(g, l, rowMid, rowEnd, colStart,
                                           colMid, skip, target_min),
                             makeLargeGate(g, l, rowMid, rowEnd, colMid, colEnd,
                                           skip, target_min)});
}

#ifdef isMPI
mEdge getMPIGate(mEdge root, int row, int col, int world_size) {
    if (root.isTerminal() || world_size <= 1) {
        return root;
    }

    // Check world_size is 2^n
    assert((world_size & world_size - 1) == 0);
    assert(row < world_size && col < world_size);

    /*
    sub-matrix index
        0|1
        ---
        2|3
    */
    int index = 0;
    int border = world_size / 2;
    if (row >= border) {
        index += 2;
    }
    if (col >= border) {
        index += 1;
    }
    mEdge tmp = root.getNode()->children[index];
    tmp.w *= root.w;
    return getMPIGate(tmp, row % border, col % border, border);
}
#endif

int vNode_to_vec(vNode *node, std::vector<vContent> &table,
                 std::unordered_map<vNode *, int> &map);
vNode *vec_to_vNode(std::vector<vContent> &table, vNodeTable &uniqTable,
                    bool no_lookup = false);
int mNode_to_vec(mNode *node, std::vector<mContent> &table,
                 std::unordered_map<mNode *, int> &map);
mNode *vec_to_mNode(std::vector<mContent> &table, mNodeTable &uniqTable,
                    bool no_lookup = false);

#ifdef isMPI
vEdge mv_multiply_MPI_org(mEdge lhs, vEdge rhs, bmpi::communicator &world) {
    int row = world.rank();
    int world_size = world.size();
    int left_neighbor = (world.rank() - 1) % world_size;
    int right_neighbor = (world.rank() + 1) % world_size;

    std_complex send_w, recv_w;
    std::vector<vContent> send_buffer, recv_buffer;
    std::unordered_map<vNode *, int> rhs_map;

    send_w = rhs.w;
    vNode_to_vec(rhs.n, send_buffer, rhs_map);
    mEdge gate = getMPIGate(lhs, row, row, world_size);
    vEdge result = mv_multiply(gate, rhs);

    for (int i = 1; i < world_size; i++) {
        std::vector<bmpi::request> recv_reqs;
        std::vector<bmpi::request> send_reqs;
        recv_buffer.clear();
        send_reqs.push_back(world.isend(right_neighbor, 2 * i - 2, send_w));
        send_reqs.push_back(
            world.isend(right_neighbor, 2 * i - 1, send_buffer));
        recv_reqs.push_back(world.irecv(left_neighbor, 2 * i - 2, recv_w));
        recv_reqs.push_back(world.irecv(left_neighbor, 2 * i - 1, recv_buffer));
        int col = (row - i + world.size()) % world_size;
        gate = getMPIGate(lhs, row, col, world_size);
        bmpi::wait_all(std::begin(recv_reqs), std::end(recv_reqs));
        vEdge received = {recv_w, vec_to_vNode(recv_buffer, vUnique)};
        result = vv_add(result, mv_multiply(gate, received));
        bmpi::wait_all(std::begin(send_reqs), std::end(send_reqs));
        send_buffer = recv_buffer;
        send_w = recv_w;
    }
    return result;
}

vEdge mv_multiply_MPI(mEdge lhs, vEdge rhs, bmpi::communicator &world,
                      std::size_t total_qubits, std::size_t largest_qubit) {
    int row = world.rank();
    int world_size = world.size();

    std_complex send_w, recv_w;
    // std::vector<vContent> send_buffer, recv_buffer;
    std::pair<std_complex, std::vector<vContent>> send_data, recv_data;
    std::unordered_map<vNode *, int> rhs_map;

    send_data.first = rhs.w;
    mEdge gate = getMPIGate(lhs, row, row, world_size);
    vEdge result = mv_multiply(gate, rhs);
    if (largest_qubit < total_qubits - int(log2(world_size))) {
        return result;
    }
    vNode_to_vec(rhs.n, send_data.second, rhs_map);
    int small_world_size =
        std::pow(2, 1 + largest_qubit - total_qubits + int(log2(world_size)));
    int offset = row - (row + world_size) % small_world_size;
    int left_neighbor = (row - 1 + world_size) % small_world_size + offset;
    int right_neighbor = (row + 1 + world_size) % small_world_size + offset;

    for (int i = 1; i < small_world_size; i++) {
        // std::vector<bmpi::request> recv_reqs;
        std::vector<bmpi::request> send_reqs;
        send_reqs.push_back(world.isend(right_neighbor, i, send_data));
        world.recv(left_neighbor, i, recv_data);
        int col = (row - i + world.size()) % small_world_size + offset;
        gate = getMPIGate(lhs, row, col, world_size);
        // bmpi::wait_all(std::begin(recv_reqs), std::end(recv_reqs));
        vEdge received = {recv_data.first,
                          vec_to_vNode(recv_data.second, vUnique)};
        result = vv_add(result, mv_multiply(gate, received));
        bmpi::wait_all(std::begin(send_reqs), std::end(send_reqs));
        send_data = recv_data;
    }
    return result;
}

vEdge mv_multiply_MPI_new(mEdge lhs, vEdge rhs, bmpi::communicator &world) {
    int row = world.rank();
    int world_size = world.size();
    int left_neighbor = (world.rank() - 1) % world_size;
    int right_neighbor = (world.rank() + 1) % world_size;

    std::vector<std_complex> buffer_w(world_size);
    buffer_w[0] = rhs.w;

    std::vector<std::vector<vContent>> buffer_table(world_size);
    std::unordered_map<vNode *, int> rhs_map;
    vNode_to_vec(rhs.n, buffer_table[0], rhs_map);

    mEdge gate = getMPIGate(lhs, row, row, world_size);
    vEdge result = mv_multiply(gate, rhs);

    for (int i = 1; i < world_size; i++) {
        std::vector<bmpi::request> recv_reqs;
        world.isend(right_neighbor, 2 * i - 2, buffer_w[i - 1]);
        world.isend(right_neighbor, 2 * i - 1, buffer_table[i - 1]);
        recv_reqs.push_back(world.irecv(left_neighbor, 2 * i - 2, buffer_w[i]));
        recv_reqs.push_back(
            world.irecv(left_neighbor, 2 * i - 1, buffer_table[i]));
        bmpi::wait_all(std::begin(recv_reqs), std::end(recv_reqs));
    }

    for (int i = 1; i < world_size; i++) {
        int col = (row - i + world.size()) % world_size;
        gate = getMPIGate(lhs, row, col, world_size);
        vEdge received = {buffer_w[i], vec_to_vNode(buffer_table[i], vUnique)};
        result = vv_add(result, mv_multiply(gate, received));
    }

    return result;
}

vEdge mv_multiply_MPI_bcast(mEdge lhs, vEdge rhs, bmpi::communicator &world) {
    int row = world.rank();
    int world_size = world.size();
    int left_neighbor = (world.rank() - 1) % world_size;
    int right_neighbor = (world.rank() + 1) % world_size;

    std::vector<std_complex> buffer_w(world_size);
    buffer_w[row] = rhs.w;

    std::vector<std::vector<vContent>> buffer_table(world_size);
    std::unordered_map<vNode *, int> rhs_map;
    vNode_to_vec(rhs.n, buffer_table[row], rhs_map);

    mEdge gate = getMPIGate(lhs, row, row, world_size);
    vEdge result = mv_multiply(gate, rhs);

    for (int i = 0; i < world_size; i++) {
        bmpi::broadcast(world, buffer_w[i], i);
        bmpi::broadcast(world, buffer_table[i], i);
    }

    for (int i = 0; i < world_size; i++) {
        int col = i;
        if (col == row) continue;
        gate = getMPIGate(lhs, row, col, world_size);
        vEdge received = {buffer_w[i], vec_to_vNode(buffer_table[i], vUnique)};
        result = vv_add(result, mv_multiply(gate, received));
    }

    return result;
}

vEdge mv_multiply_MPI_bcast2(mEdge lhs, vEdge rhs, bmpi::communicator &world) {
    int row = world.rank();
    int world_size = world.size();
    int left_neighbor = (world.rank() - 1) % world_size;
    int right_neighbor = (world.rank() + 1) % world_size;

    // prepare data to be sent
    std_complex send_w = rhs.w;
    std::vector<vContent> send_table;
    std::unordered_map<vNode *, int> rhs_map;
    vNode_to_vec(rhs.n, send_table, rhs_map);

    // calculate initial result
    mEdge gate = getMPIGate(lhs, row, row, world_size);
    vEdge result = mv_multiply(gate, rhs);

    for (int i = 0; i < world_size; i++) {
        if (row == i) {
            bmpi::broadcast(world, send_w, i);
            bmpi::broadcast(world, send_table, i);
        } else {
            std_complex buffer_w;
            std::vector<vContent> buffer_table;
            bmpi::broadcast(world, buffer_w, i);
            bmpi::broadcast(world, buffer_table, i);
            gate = getMPIGate(lhs, row, i, world_size);
            vEdge received = {buffer_w, vec_to_vNode(buffer_table, vUnique)};
            result = vv_add(result, mv_multiply(gate, received));
        }
    }

    return result;
}

vEdge mv_multiply_MPI_bcast3(mEdge lhs, vEdge rhs, bmpi::communicator &world,
                             std::size_t total_qubits,
                             std::size_t largest_qubit) {
    int row = world.rank();
    int world_size = world.size();
    int left_neighbor = (world.rank() - 1) % world_size;
    int right_neighbor = (world.rank() + 1) % world_size;

    // calculate initial result
    mEdge gate = getMPIGate(lhs, row, row, world_size);
    vEdge result = mv_multiply(gate, rhs);
    if (largest_qubit < total_qubits - int(log2(world_size))) {
        return result;
    }

    // prepare data to be sent
    std::pair<std_complex, std::vector<vContent>> send_data;
    send_data.first = rhs.w;
    std::unordered_map<vNode *, int> rhs_map;
    vNode_to_vec(rhs.n, send_data.second, rhs_map);

    for (int i = 0; i < world_size; i++) {
        if (row == i) {
            bmpi::broadcast(world, send_data, i);
        } else {
            std::pair<std_complex, std::vector<vContent>> recv_data;
            bmpi::broadcast(world, recv_data, i);
            gate = getMPIGate(lhs, row, i, world_size);
            vEdge received = {recv_data.first,
                              vec_to_vNode(recv_data.second, vUnique)};
            result = vv_add(result, mv_multiply(gate, received));
        }
    }

    return result;
}

#endif

static Qubit rootVar(const mEdge &lhs, const mEdge &rhs) {
    assert(!(lhs.isTerminal() && rhs.isTerminal()));

    return (lhs.isTerminal() ||
            (!rhs.isTerminal() && (rhs.getVar() > lhs.getVar())))
               ? rhs.getVar()
               : lhs.getVar();
}

mEdge mm_add2(const mEdge &lhs, const mEdge &rhs, int32_t current_var) {
    if (lhs.w.isApproximatelyZero()) {
        return rhs;
    } else if (rhs.w.isApproximatelyZero()) {
        return lhs;
    }

    if (current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, mNode::terminal};
    }
    if (lhs.n == rhs.n) {
        return {lhs.w + rhs.w, lhs.n};
    }

    mEdge result;

#if defined(isMT) && !defined(CACHE_GLOBAL)
    AddCache &_aCache = _aCaches.local();
#endif
    result = _aCache.find(lhs, rhs);
    if (result.n != nullptr) {
        if (result.w.isApproximatelyZero()) {
            return mEdge::zero;
        } else {
            return result;
        }
    }

    mEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    mNode *lnode = lhs.getNode();
    mNode *rnode = rhs.getNode();

    std::array<mEdge, 4> edges;

    for (auto i = 0; i < 4; i++) {
        if (lv == current_var && !lhs.isTerminal()) {
            x = lnode->getEdge(i);
            x.w = lhs.w * x.w;
        } else {
            if(i==1 || i==2){
                x = mEdge::zero;
            }else{
                x = lhs;
            }
        }
        if (rv == current_var && !rhs.isTerminal()) {
            y = rnode->getEdge(i);
            y.w = rhs.w * y.w;
        } else {
            if(i==1 || i==2){
                y = mEdge::zero;
            }else{
                y = rhs;
            }
        }

        edges[i] = mm_add2(x, y, current_var - 1);
    }

    result = makeMEdge(current_var, edges);
    _aCache.set(lhs, rhs, result);

    return result;
}

mEdge mm_add(const mEdge &lhs, const mEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w + rhs.w, mNode::terminal};
    }

    Qubit root = rootVar(lhs, rhs);
    return mm_add2(lhs, rhs, root);
}

mEdge mm_multiply2(const mEdge &lhs, const mEdge &rhs, int32_t current_var) {
    if (lhs.w.isApproximatelyZero() || rhs.w.isApproximatelyZero()) {
        return mEdge::zero;
    }

    if (current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, mNode::terminal};
    }

#if defined(isMT) && !defined(CACHE_GLOBAL)
    MulCache &_mCache = _mCaches.local();
#endif
    mEdge result;
    result = _mCache.find(lhs.n, rhs.n);
    if (result.n != nullptr) {
        if (result.w.isApproximatelyZero()) {
            return mEdge::zero;
        } else {
            result.w = result.w * lhs.w * rhs.w;
            if (result.w.isApproximatelyZero())
                return mEdge::zero;
            else
                return result;
        }
    }

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    assert(lv <= current_var && rv <= current_var);
    mNode *lnode = lhs.getNode();
    mNode *rnode = rhs.getNode();
    mEdge x, y;
    mEdge lcopy = lhs;
    mEdge rcopy = rhs;
    lcopy.w = {1.0, 0.0};
    rcopy.w = {1.0, 0.0};

    std::array<mEdge, 4> edges;

    for (auto i = 0; i < 4; i++) {
        std::size_t row = i >> 1;
        std::size_t col = i & 0x1;

        std::array<mEdge, 2> product;
        for (auto k = 0; k < 2; k++) {
            if (lv == current_var && !lhs.isTerminal()) {
                x = lnode->getEdge((row << 1) | k);
            } else {
                auto index = row << 1 | k;
                if(index == 1 || index == 2){
                    x = mEdge::zero;
                }else{
                    x = lcopy;
                }
            }

            if (rv == current_var && !rhs.isTerminal()) {
                y = rnode->getEdge((k << 1) | col);
            } else {
                auto index = k << 1 | col;
                if(index == 1 || index == 2){
                    y = mEdge::zero;
                }else{
                    y = rcopy;
                }
            }

            product[k] = mm_multiply2(x, y, current_var - 1);
        }
        edges[i] = mm_add2(product[0], product[1], current_var - 1);
    }

    result = makeMEdge(current_var, edges);
    _mCache.set(lhs.n, rhs.n, result);

    result.w = result.w * lhs.w * rhs.w;
    if (result.w.isApproximatelyZero()) return mEdge::zero;
    if (result.w.isApproximatelyOne()) result.w = {1.0, 0.0};
    return result;
}

mEdge mm_multiply(const mEdge &lhs, const mEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w * rhs.w, mNode::terminal};
    }

    Qubit root = rootVar(lhs, rhs);
    mEdge result = mm_multiply2(lhs, rhs, root);
    return result;
}

static void printVector2(const vEdge &edge, std::size_t row,
                         const std_complex &w, uint64_t left, std_complex *m) {
    std_complex wp = edge.w * w;

    if (edge.isTerminal() && left == 0) {
        m[row] = wp;
        return;
    } else if (edge.isTerminal()) {
        row = row << left;

        for (std::size_t i = 0; i < (1 << left); i++) {
            m[row | i] = wp;
        }
        return;
    }

    vNode *node = edge.getNode();
    printVector2(node->getEdge(0), (row << 1) | 0, wp, left - 1, m);
    printVector2(node->getEdge(1), (row << 1) | 1, wp, left - 1, m);
}

mEdge mm_kronecker2(const mEdge &lhs, const mEdge &rhs) {
    if (lhs.isTerminal()) {
        return {lhs.w * rhs.w, rhs.n};
    }

    std::array<mEdge, 4> edges;
    mEdge x;
    mNode *lnode = lhs.getNode();

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    for (auto i = 0; i < 4; i++) {
        x = lnode->getEdge(i);
        x.w = lhs.w * x.w;
        edges[i] = mm_kronecker2(x, rhs);
    }

    mEdge ret = makeMEdge(lv + rv + 1, edges);
    return ret;
}

mEdge mm_kronecker(const mEdge &lhs, const mEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w * rhs.w, mNode::terminal};
    }

    return mm_kronecker2(lhs, rhs);
}

vEdge vv_add2(const vEdge &lhs, const vEdge &rhs, int32_t current_var) {
    if (lhs.w.isApproximatelyZero()) {
        if (rhs.w.isApproximatelyZero()) {
            return vEdge::zero;
        }
        return rhs;
    } else if (rhs.w.isApproximatelyZero()) {
        return lhs;
    }

    if (current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, vNode::terminal};
    }
    if (lhs.n == rhs.n) {
        return {lhs.w + rhs.w, lhs.n};
    }

    vEdge result;

#if defined(isMT) && !defined(CACHE_GLOBAL)
    AddCache &_aCache = _aCaches.local();
#endif
    result = _aCache.find(lhs, rhs);
    if (result.n != nullptr) {
        if (result.w.isApproximatelyZero()) {
            return vEdge::zero;
        }
        return result;
    }

    vEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    vNode *lnode = lhs.getNode();
    vNode *rnode = rhs.getNode();
    std::array<vEdge, 2> edges;
#ifdef isMT
    std::vector<boost::fibers::future<vEdge>> sums;
#endif

    for (auto i = 0; i < 2; i++) {
        if (lv == current_var && !lhs.isTerminal()) {
            x = lnode->getEdge(i);
            x.w = lhs.w * x.w;
        } else {
            x = lhs;
        }
        if (rv == current_var && !rhs.isTerminal()) {
            y = rnode->getEdge(i);
            y.w = rhs.w * y.w;
        } else {
            y = rhs;
        }
#ifdef isMT
        if (current_var > LIMIT) {
            boost::fibers::packaged_task<vEdge()> pt(
                std::bind(vv_add2, x, y, current_var - 1));
            sums.emplace_back(pt.get_future());
            boost::fibers::fiber f(std::move(pt));
            f.detach();
        } else {
#endif
            edges[i] = vv_add2(x, y, current_var - 1);
#ifdef isMT
        }
#endif
    }

#ifdef isMT
    if (current_var > LIMIT) {
        assert(sums.size() == 2);
        edges[0] = sums[0].get();
        edges[1] = sums[1].get();
    }
#endif
    result = makeVEdge(current_var, edges);
    _aCache.set(lhs, rhs, result);

    return result;
}

vEdge vv_add(const vEdge &lhs, const vEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w + rhs.w, vNode::terminal};
    }

    return vv_add2(lhs, rhs, rhs.getVar());
}

vEdge vv_kronecker2(const vEdge &lhs, const vEdge &rhs) {
    if (lhs.isTerminal()) {
        return {lhs.w * rhs.w, rhs.n};
    }
    if (rhs.isTerminal()) {
        return {lhs.w * rhs.w, lhs.n};
    }

    std::array<vEdge, 2> edges;
    vEdge x;
    vNode *lnode = lhs.getNode();

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    for (auto i = 0; i < 2; i++) {
        x = lnode->getEdge(i);
        x.w = lhs.w * x.w;
        edges[i] = vv_kronecker2(x, rhs);
    }

    vEdge ret = makeVEdge(lv + rv + 1, edges);
    return ret;
}

vEdge vv_kronecker(const vEdge &lhs, const vEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w * rhs.w, vNode::terminal};
    }

    return vv_kronecker2(lhs, rhs);
}

void vEdge::printVector() const {
    if (this->isTerminal()) {
        std::cout << this->w << std::endl;
        return;
    }
    Qubit q = this->getVar();
    std::size_t dim = 1 << (q + 1);

    std_complex *vector = new std_complex[dim];

    printVector2(*this, 0, {1.0, 0.0}, q + 1, vector);

    for (size_t i = 0; i < dim; i++) {
        std::cout << vector[i] << " ";
        std::cout << "\n";
    }
    std::cout << std::endl;

    delete[] vector;
}

#ifdef isMPI
void vEdge::printVectorMPI(bmpi::communicator &world) const {
    std::stringstream ss;
    if (this->isTerminal()) {
        ss << "(all 0)" << std::endl;
    } else {
        Qubit q = this->getVar();
        std::size_t dim = 1 << (q + 1);
        std_complex *vector = new std_complex[dim];
        printVector2(*this, 0, {1.0, 0.0}, q + 1, vector);
        for (size_t i = 0; i < dim; i++) {
            ss << vector[i] << std::endl;
        }
        delete[] vector;
    }

    if (world.rank() == 0) {
        std::cout << ss.str();
        for (int i = 1; i < world.size(); i++) {
            std::string msg;
            world.recv(i, i, msg);
            std::cout << msg;
        }
    } else {
        world.send(0, world.rank(), ss.str());
    }
    return;
}
#endif

static void printVector_sparse2(const vEdge &edge, std::size_t row,
                                const std_complex &w, uint64_t left,
                                std::map<int, std_complex> &m) {
    std_complex wp = edge.w * w;

    const double thr = 0.001;

    if (edge.isTerminal() && left == 0) {
        if (norm(wp) > thr) m[row] = wp;
        return;
    } else if (edge.isTerminal()) {
        if (norm(wp) > thr) {
            row = row << left;
            for (std::size_t i = 0; i < (1 << left); i++) {
                m[row | i] = wp;
            }
        }
        return;
    }

    vNode *node = edge.getNode();
    printVector_sparse2(node->getEdge(0), (row << 1) | 0, wp, left - 1, m);
    printVector_sparse2(node->getEdge(1), (row << 1) | 1, wp, left - 1, m);
}

void vEdge::printVector_sparse() const {
    if (this->isTerminal()) {
        std::cout << this->w << std::endl;
        return;
    }
    Qubit q = this->getVar();
    std::size_t dim = 1 << (q + 1);

    std::map<int, std_complex> map;

    printVector_sparse2(*this, 0, {1.0, 0.0}, q + 1, map);

    for (auto itr = map.begin(); itr != map.end(); itr++) {
        std::cout << std::bitset<30>(itr->first) << ": " << itr->second
                  << std::endl;
    }
}

static void fillVector(const vEdge &edge, std::size_t row, const std_complex &w,
                       uint64_t left, std_complex *m) {
    std_complex wp = edge.w * w;

    if (edge.isTerminal() && left == 0) {
        m[row] = wp;
        return;
    } else if (edge.isTerminal()) {
        row = row << left;

        for (std::size_t i = 0; i < (1 << left); i++) {
            m[row | i] = wp;
        }
        return;
    }

    vNode *node = edge.getNode();
    fillVector(node->getEdge(0), (row << 1) | 0, wp, left - 1, m);
    fillVector(node->getEdge(1), (row << 1) | 1, wp, left - 1, m);
}

std_complex *vEdge::getVector(std::size_t *dim) const {
    assert(!this->isTerminal());

    Qubit q = this->getVar();
    std::size_t d = 1 << (q + 1);

    std_complex *vector = new std_complex[d];
    fillVector(*this, 0, {1.0, 0.0}, q + 1, vector);
    if (dim != nullptr) *dim = d;
    return vector;
}

VectorXcf vEdge::getEigenVector() {
    std::size_t dim;
    auto v = getVector(&dim);
    VectorXcf V(dim);

    for (auto i = 0; i < dim; i++) {
        V(i) = std::complex<double>{v[i].r, v[i].i};
    }
    return V;
}

vEdge mv_multiply2(const mEdge &lhs, const vEdge &rhs, int32_t current_var) {
    if (lhs.w.isApproximatelyZero() || rhs.w.isApproximatelyZero()) {
        return vEdge::zero;
    }

    if (current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, vNode::terminal};
    }

#if defined(isMT) && !defined(CACHE_GLOBAL)
    MulCache &_mCache = _mCaches.local();
#endif
    vEdge result;

    result = _mCache.find(lhs.n, rhs.n);
    if (result.n != nullptr) {
        if (result.w.isApproximatelyZero()) {
            return vEdge::zero;
        } else {
            result.w = result.w * lhs.w * rhs.w;
            if (result.w.isApproximatelyZero())
                return vEdge::zero;
            else
                return result;
        }
    }

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    mNode *lnode = lhs.getNode();
    vNode *rnode = rhs.getNode();
    mEdge x;
    vEdge y;
    mEdge lcopy = lhs;
    vEdge rcopy = rhs;
    lcopy.w = {1.0, 0.0};
    rcopy.w = {1.0, 0.0};

    std::array<vEdge, 2> edges;

    for (auto i = 0; i < 2; i++) {
#ifdef isMT
        std::vector<boost::fibers::future<vEdge>> products;
#endif
        std::array<vEdge, 2> product;
        for (auto k = 0; k < 2; k++) {
            if (lv == current_var && !lhs.isTerminal()) {
                x = lnode->getEdge((i << 1) | k);
            } else {
                if( ((i << 1) | k)==0 || ((i << 1) | k)==3 )
                    x = lcopy;
                else
                    x = mEdge::zero;
            }

            if (rv == current_var && !rhs.isTerminal()) {
                y = rnode->getEdge(k);
            } else {
                y = rcopy;
            }
#ifdef isMT
            if (current_var > LIMIT) {
                boost::fibers::packaged_task<vEdge()> pt(
                    std::bind(mv_multiply2, x, y, current_var - 1));
                products.emplace_back(pt.get_future());
                boost::fibers::fiber f(std::move(pt));
                f.detach();
            } else {
#endif
                product[k] = mv_multiply2(x, y, current_var - 1);
#ifdef isMT
            }
#endif
        }
#ifdef isMT
        if (current_var > LIMIT) {
            product[0] = products[0].get();
            product[1] = products[1].get();
        }
#endif
        edges[i] = vv_add2(product[0], product[1], current_var - 1);
    }

    result = makeVEdge(current_var, edges);
    _mCache.set(lhs.n, rhs.n, result);
    result.w = result.w * lhs.w * rhs.w;
    if (result.w.isApproximatelyZero()) {
        return vEdge::zero;
    }
    if (result.w.isApproximatelyOne()) {
        result.w = {1.0, 0.0};
    }
    return result;
}

vEdge mv_multiply(mEdge lhs, vEdge rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w * rhs.w, vNode::terminal};
    }
    LIMIT = rhs.getVar() - MINUS;
    vEdge v = mv_multiply2(lhs, rhs, rhs.getVar());
    //_mCache.hitRatio();
    //_aCache.hitRatio();
    // vUnique.dump();
    // mUnique.dump();
    // genDot(v);
    return v;
}

double assignProbabilities(const vEdge &edge) {
    auto it = probs.find(edge.n);
    if (it != probs.end()) {
        return edge.w.mag2() * it->second;
    }
    double sum{1};
    if (!edge.isTerminal()) {
        sum = assignProbabilities(edge.n->children.at(0)) +
              assignProbabilities(edge.n->children.at(1));
    }

    probs.insert({edge.n, sum});

    return edge.w.mag2() * sum;
}

std::vector<double> probabilities(const vEdge &rootEdge) {
    const auto nqubits = static_cast<QubitCount>(rootEdge.getVar() + 1);
    assignProbabilities(rootEdge);
    std::vector<double> result;
    for (auto i = 0; i < pow(2, nqubits); i++) {
        vEdge cur = rootEdge;
        double prob = 1.0;
        for (Qubit j = rootEdge.getVar(); j >= 0; j--) {
            double p0 = probs[cur.n->getEdge(0).n] * cur.n->getEdge(0).w.mag2();
            double p1 = probs[cur.n->getEdge(1).n] * cur.n->getEdge(1).w.mag2();
            double tmp = p0 + p1;
            if (tmp > 0) {
                p0 /= tmp;
            }
            if ((i >> j) % 2 == 0) {
                prob *= p0;
            } else {
                prob *= (1 - p0);
            }
            if (prob == 0) {
                continue;
            }
            cur = cur.n->children.at((i >> j) % 2);
        }
        result.push_back(prob);
    }
    return result;
}

template <typename T>
std::vector<T> flatten(std::vector<std::vector<T>> const &vec) {
    std::vector<T> flattened;
    for (auto const &v : vec) {
        flattened.insert(flattened.end(), v.begin(), v.end());
    }
    return flattened;
}

#ifdef isMPI
std::vector<double> probabilitiesMPI(bmpi::communicator &world,
                                     const vEdge &rootEdge) {
    const double each_prob = assignProbabilities(rootEdge);
    std::vector<double> each_result = probabilities(rootEdge);
    const double all_prob =
        bmpi::all_reduce(world, each_prob, std::plus<double>());
    std::for_each(each_result.begin(), each_result.end(),
                  [each_prob, all_prob](double &num) {
                      num = num * each_prob / all_prob;
                  });

    std::vector<std::vector<double>> all_result;
    bmpi::all_gather(world, each_result, all_result);
    auto final_result = flatten(all_result);

    return final_result;
}

std::string measureAllMPI(bmpi::communicator &world, vEdge &rootEdge,
                          const bool collapse, std::mt19937_64 &mt,
                          double epsilon) {
    assert(collapse == false);
    if (probs.find(rootEdge.n) == probs.end()) {
        assignProbabilities(rootEdge);
    }
    double each_prob = probs[rootEdge.n] * rootEdge.w.mag2();
    std::string each_result = measureAll(rootEdge, collapse, mt, epsilon);
    std::string result = "";
    std::vector<double> vec_prob;
    std::vector<std::string> vec_strings;

    if (world.rank() == 0) {
        bmpi::gather(world, each_prob, vec_prob, 0);
        bmpi::gather(world, each_result, vec_strings, 0);
        double sum = 0;
        for (int i = 0; i < world.size(); i++) {
            sum += vec_prob[i];
        }
        // std::cout << "sum=" << sum << std::endl;
        std::uniform_real_distribution<double> dist(0.0, sum);
        const double threshold = dist(mt);
        double current = 0;
        int target_rank = 0;
        for (int i = 0; i < world.size(); i++) {
            current += vec_prob[i];
            target_rank = i;
            if (threshold <= current) {
                break;
            }
        }
        const int pre_length = log2(world.size());
        // std::cout << "thre="<< threshold << " target_rank=" << target_rank <<
        // " " << std::endl;
        int tmp = target_rank;
        for (int i = pre_length - 1; i >= 0; i--) {
            if (tmp / (int)pow(2, i) >= 1) {
                tmp -= (int)pow(2, i);
                result += '1';
            } else {
                result += '0';
            }
        }
        result += vec_strings[target_rank];
    } else {
        bmpi::gather(world, each_prob, vec_prob, 0);
        bmpi::gather(world, each_result, vec_strings, 0);
    }
    bmpi::broadcast(world, result, 0);
    return result;
}
#endif

std::string measureAll(vEdge &rootEdge, const bool collapse,
                       std::mt19937_64 &mt, double epsilon) {
    vEdge cur = rootEdge;
    const auto nqubits = static_cast<QubitCount>(rootEdge.getVar() + 1);
    std::string result(nqubits, '0');
    std::uniform_real_distribution<double> dist(0.0, 1.0L);

    if (probs.find(rootEdge.n) == probs.end()) {
        assignProbabilities(rootEdge);
    }

    for (Qubit i = rootEdge.getVar(); i >= 0; --i) {
        double p0 = probs[cur.n->getEdge(0).n] * cur.n->getEdge(0).w.mag2();
        double p1 = probs[cur.n->getEdge(1).n] * cur.n->getEdge(1).w.mag2();
        double tmp = p0 + p1;

        if (tmp != 0) {
            p0 /= tmp;
        } else {
            p0 = 0.5;
        }

        const double threshold = dist(mt);

        if (threshold < p0) {
            cur = cur.n->getEdge(0);
        } else {
            result[cur.n->v] = '1';
            cur = cur.n->getEdge(1);
        }
    }

    if (collapse) {
        vEdge e = vEdge::one;
        std::array<vEdge, 2> edges{};
        for (Qubit p = 0; p < nqubits; p++) {
            if (result[p] == '0') {
                edges[0] = e;
                edges[1] = vEdge::zero;
            } else {
                edges[0] = vEdge::zero;
                edges[1] = e;
            }
            e = makeVEdge(p, edges);
        }
        rootEdge = e;
    }
    return std::string{result.rbegin(), result.rend()};
}

std::pair<double, double> determineMeasurementProbabilities(
    const vEdge &rootEdge, const Qubit index) {
    std::unordered_map<vNode *, double> probsMone;
    std::set<vNode *> visited;
    std::queue<vNode *> q;

    probsMone[rootEdge.n] = rootEdge.w.mag2();
    visited.insert(rootEdge.n);
    q.push(rootEdge.n);

    while (q.front()->v != index) {
        vNode *ptr = q.front();
        q.pop();
        const double prob = probsMone[ptr];

        if (!ptr->children.at(0).w.isApproximatelyZero()) {
            const double tmp1 = prob * ptr->children.at(0).w.mag2();

            if (visited.find(ptr->children.at(0).n) != visited.end()) {
                probsMone[ptr->children.at(0).n] =
                    probsMone[ptr->children.at(0).n] + tmp1;
            } else {
                probsMone[ptr->children.at(0).n] = tmp1;
                visited.insert(ptr->children.at(0).n);
                q.push(ptr->children.at(0).n);
            }
        }

        if (!ptr->children.at(1).w.isApproximatelyZero()) {
            const double tmp1 = prob * ptr->children.at(1).w.mag2();

            if (visited.find(ptr->children.at(1).n) != visited.end()) {
                probsMone[ptr->children.at(1).n] =
                    probsMone[ptr->children.at(1).n] + tmp1;
            } else {
                probsMone[ptr->children.at(1).n] = tmp1;
                visited.insert(ptr->children.at(1).n);
                q.push(ptr->children.at(1).n);
            }
        }
    }

    double pzero{0};
    double pone{0};

    if (probs.find(rootEdge.n) == probs.end()) {
        assignProbabilities(rootEdge);
    }

    while (!q.empty()) {
        vNode *ptr = q.front();
        q.pop();

        if (!ptr->children.at(0).w.isApproximatelyZero()) {
            pzero += probsMone[ptr] * probs[ptr->children.at(0).n] *
                     ptr->children.at(0).w.mag2();
        }

        if (!ptr->children.at(1).w.isApproximatelyZero()) {
            pone += probsMone[ptr] * probs[ptr->children.at(1).n] *
                    ptr->children.at(1).w.mag2();
        }
    }
    return {pzero, pone};
}

#ifdef isMPI
char measureOneCollapsingMPI(bmpi::communicator &world, vEdge &rootEdge,
                             const Qubit index, const Qubit n_qubits,
                             std::mt19937_64 &mt, bool collapse,
                             double epsilon) {
    const int global_size = log2(world.size());     // e.g. 1
    const int local_size = n_qubits - global_size;  // e.g. 3-1=2
    const int rank = world.rank();
    char result = '1';
    // std::cout << "global=" << global_size << " local=" << local_size << "
    // index=" << index << " n_qubits=" << n_qubits << std::endl;
    double pzero = 0.0, pone = 0.0;

    double p_each = assignProbabilities(rootEdge);

    if (index >= local_size) {
        int id = (rank << index) & (1 << index);
        if (id == 0) {
            pzero = p_each;
        } else {
            pone = p_each;
        }
        // std::cout << "Global: rank=" << rank << " (" << pzero << ", " << pone
        // << ")" << std::endl;
    } else {
        if (p_each > 0) {
            auto tmp = determineMeasurementProbabilities(rootEdge, index);
            pzero = tmp.first;
            pone = tmp.second;
        }
        // std::cout << "Local: rank=" << rank << " (" << pzero << ", " << pone
        // << ")" << std::endl;
    }

    double pzero_all = bmpi::all_reduce(world, pzero, std::plus<double>());
    double pone_all = bmpi::all_reduce(world, pone, std::plus<double>());
    assert(pzero_all != 0.0 || pone_all != 0.0);

    if (rank == 0) {
        std::uniform_real_distribution<double> dist(0.0, pzero_all + pone_all);
        double threshold = dist(mt);
        if (threshold < pzero_all) {
            result = '0';
        }
    }

    bmpi::broadcast(world, result, 0);
    if (collapse == false) {
        return result;
    }

    double normalizationFactor;
    GateMatrix measurementMatrix{cf_zero, cf_zero, cf_zero, cf_zero};
    if (result == '0') {
        measurementMatrix[0] = cf_one;
        normalizationFactor = pzero_all;
    } else {
        measurementMatrix[3] = cf_one;
        normalizationFactor = pone_all;
    }
    assert(normalizationFactor > 0);
    mEdge measurementGate = makeGate(n_qubits, measurementMatrix, index);
    vEdge e =
        mv_multiply_MPI(measurementGate, rootEdge, world, n_qubits, index);
    std_complex c = {std::sqrt((pzero_all + pone_all) / normalizationFactor),
                     0};
    c = e.w * c;
    e.w = c;
    rootEdge = e;

    return result;
}

double measureOneMPI(bmpi::communicator &world, vEdge &rootEdge,
                     const Qubit index, const Qubit n_qubits,
                     std::mt19937_64 &mt, double epsilon) {
    const int global_size = log2(world.size());     // e.g. 1
    const int local_size = n_qubits - global_size;  // e.g. 3-1=2
    const int rank = world.rank();
    char result = '1';
    // std::cout << "global=" << global_size << " local=" << local_size << "
    // index=" << index << " n_qubits=" << n_qubits << std::endl;
    double pzero = 0.0, pone = 0.0;

    double p_each = assignProbabilities(rootEdge);

    if (index >= local_size) {
        int id = (rank << index) & (1 << index);
        if (id == 0) {
            pzero = p_each;
        } else {
            pone = p_each;
        }
        // std::cout << "Global: rank=" << rank << " (" << pzero << ", " << pone
        // << ")" << std::endl;
    } else {
        if (p_each > 0) {
            auto tmp = determineMeasurementProbabilities(rootEdge, index);
            pzero = tmp.first;
            pone = tmp.second;
        }
        // std::cout << "Local: rank=" << rank << " (" << pzero << ", " << pone
        // << ")" << std::endl;
    }

    double pzero_all = bmpi::all_reduce(world, pzero, std::plus<double>());
    double pone_all = bmpi::all_reduce(world, pone, std::plus<double>());
    assert(pzero_all != 0.0 || pone_all != 0.0);

    return pzero_all / (pzero_all + pone_all);
}
#endif

char measureOneCollapsing(vEdge &rootEdge, const Qubit index,
                          std::mt19937_64 &mt, double epsilon) {
    const auto &[pzero, pone] =
        determineMeasurementProbabilities(rootEdge, index);
    const double sum = pzero + pone;
    if (std::abs(sum - 1) > epsilon) {
        throw std::runtime_error(
            "Numerical instability occurred during measurement: |alpha|^2 + "
            "|beta|^2 = " +
            std::to_string(pzero) + " + " + std::to_string(pone) + " = " +
            std::to_string(pzero + pone) + ", but should be 1!");
    }
    GateMatrix measurementMatrix{cf_zero, cf_zero, cf_zero, cf_zero};

    std::uniform_real_distribution<double> dist(0.0, 1.0L);

    double threshold = dist(mt);
    double
        normalizationFactor;  // NOLINT(cppcoreguidelines-init-variables) always
                              // assigned a value in the following block
    char result;  // NOLINT(cppcoreguidelines-init-variables) always assigned a
                  // value in the following block

    if (threshold < pzero / sum) {
        measurementMatrix[0] = cf_one;
        normalizationFactor = pzero;
        result = '0';
    } else {
        measurementMatrix[3] = cf_one;
        normalizationFactor = pone;
        result = '1';
    }

    mEdge measurementGate =
        makeGate(rootEdge.getVar() + 1, measurementMatrix, index);

    vEdge e = mv_multiply(measurementGate, rootEdge);

    std_complex c = {std::sqrt(1.0 / normalizationFactor), 0};
    c = e.w * c;
    e.w = c;
    rootEdge = e;

    return result;
}

double measureOne(vEdge &rootEdge, const Qubit index, std::mt19937_64 &mt,
                  double epsilon) {
    const auto &[pzero, pone] =
        determineMeasurementProbabilities(rootEdge, index);
    const double sum = pzero + pone;
    if (std::abs(sum - 1) > epsilon) {
        throw std::runtime_error(
            "Numerical instability occurred during measurement: |alpha|^2 + "
            "|beta|^2 = " +
            std::to_string(pzero) + " + " + std::to_string(pone) + " = " +
            std::to_string(pzero + pone) + ", but should be 1!");
    }

    return pzero / sum;
}

mEdge RX(QubitCount qnum, int target, double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), 0};
    std::complex<double> i2 = {0, -std::sin(angle / 2)};
    return makeGate(qnum, GateMatrix{i1, i2, i2, i1}, target);
}
mEdge RY(QubitCount qnum, int target, double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), 0};
    std::complex<double> i2 = {-std::sin(angle / 2), 0};
    std::complex<double> i3 = {std::sin(angle / 2), 0};
    return makeGate(qnum, GateMatrix{i1, i2, i3, i1}, target);
}
mEdge RZ(QubitCount qnum, int target, double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), -std::sin(angle / 2)};
    std::complex<double> i2 = {std::cos(angle / 2), std::sin(angle / 2)};
    return makeGate(qnum, GateMatrix{i1, cf_zero, cf_zero, i2}, target);
}

mEdge CX(QubitCount qnum, int target, int control) {
    std::complex<double> zero = {0, 0};
    std::complex<double> one = {1, 0};
    Controls controls;
    controls.emplace(Control{control, Control::Type::pos});
    return makeGate(qnum, GateMatrix{zero, one, one, zero}, target, controls);
}

mEdge makeSwap(QubitCount qnum, int target0, int target1) {
    Controls c1{{Control{target0, Control::Type::pos}}};
    mEdge e1 = makeGate(qnum, Xmat, target1, c1);
    Controls c2{{Control{target1, Control::Type::pos}}};
    mEdge e2 = makeGate(qnum, Xmat, target0, c2);
    mEdge e3 = mm_multiply(e2, e1);
    e3 = mm_multiply(e1, e3);
    return e3;
}

mEdge RXX(QubitCount qnum, int target0, int target1, double angle) {
    TwoQubitGateMatrix matrix = rxx_matrix(angle);
    return makeTwoQubitGate(qnum, matrix, target0, target1);
}

mEdge RYY(QubitCount qnum, int target0, int target1, double angle) {
    TwoQubitGateMatrix matrix = ryy_matrix(angle);
    return makeTwoQubitGate(qnum, matrix, target0, target1);
}

mEdge RZZ(QubitCount qnum, int target0, int target1, double angle) {
    TwoQubitGateMatrix matrix = rzz_matrix(angle);
    return makeTwoQubitGate(qnum, matrix, target0, target1);
}

mEdge RZX(QubitCount qnum, int target0, int target1, double angle) {
    TwoQubitGateMatrix matrix = rzx_matrix(angle);
    return makeTwoQubitGate(qnum, matrix, target0, target1);
}

mEdge SWAP(QubitCount qnum, int target0, int target1) {
    return makeSwap(qnum, target0, target1);
}

mEdge ISWAP(QubitCount qnum, int target0, int target1) {
    TwoQubitGateMatrix matrix = iswap_matrix();
    return makeTwoQubitGate(qnum, matrix, target0, target1);
}

mEdge CSWAP(QubitCount qnum, int target0, int target1, int control) {
    Controls controls;
    controls.emplace(Control{control, Control::Type::pos});
    TwoQubitGateMatrix matrix = swap_matrix();
    return makeTwoQubitGate(qnum, matrix, target0, target1, controls);
}

GateMatrix rx(double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), 0};
    std::complex<double> i2 = {0, -std::sin(angle / 2)};
    return GateMatrix{i1, i2, i2, i1};
}

GateMatrix ry(double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), 0};
    std::complex<double> i2 = {-std::sin(angle / 2), 0};
    std::complex<double> i3 = {std::sin(angle / 2), 0};
    return GateMatrix{i1, i2, i3, i1};
}

GateMatrix rz(double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), -std::sin(angle / 2)};
    std::complex<double> i2 = {std::cos(angle / 2), std::sin(angle / 2)};
    return GateMatrix{i1, cf_zero, cf_zero, i2};
}

GateMatrix u3(double theta, double phi, double lambda) {
    return u(theta, phi, lambda, 0);
}

GateMatrix u1(double lambda) { return u(0, 0, lambda, 0); }

GateMatrix u2(double phi, double lambda) { return u(PI / 2, phi, lambda, 0); }

GateMatrix u(double theta, double phi, double lambda, double gamma) {
    std::complex<double> i1 =
        std::exp(std::complex<double>(0, gamma)) * std::cos(theta / 2);
    std::complex<double> i2 =
        -std::exp(std::complex<double>(0, lambda + gamma)) *
        std::sin(theta / 2);
    std::complex<double> i3 =
        std::exp(std::complex<double>(0, phi + gamma)) * std::sin(theta / 2);
    std::complex<double> i4 =
        std::exp(std::complex<double>(0, lambda + phi + gamma)) *
        std::cos(theta / 2);
    return GateMatrix{i1, i2, i3, i4};
}

GateMatrix p(double angle) {
    std::complex<double> i1 = {std::cos(angle), std::sin(angle)};
    return GateMatrix{1, 0, 0, i1};
}

GateMatrix r(double theta, double phi) {
    std::complex<double> i1 = {std::cos(theta / 2), 0};
    std::complex<double> i2 = std::complex<double>(0, -1) *
                              std::exp(std::complex<double>(0, -phi)) *
                              std::sin(theta / 2);
    std::complex<double> i3 = std::complex<double>(0, -1) *
                              std::exp(std::complex<double>(0, phi)) *
                              std::sin(theta / 2);
    return GateMatrix{i1, i2, i3, i1};
}

TwoQubitGateMatrix rxx_matrix(double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), 0};
    std::complex<double> i2 = {0, -std::sin(angle / 2)};
    return TwoQubitGateMatrix{{{i1, cf_zero, cf_zero, i2},
                               {cf_zero, i1, i2, cf_zero},
                               {cf_zero, i2, i1, cf_zero},
                               {i2, cf_zero, cf_zero, i1}}};
}

TwoQubitGateMatrix ryy_matrix(double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), 0};
    std::complex<double> i2 = {0, std::sin(angle / 2)};
    return TwoQubitGateMatrix{{{i1, cf_zero, cf_zero, i2},
                               {cf_zero, i1, -i2, cf_zero},
                               {cf_zero, -i2, i1, cf_zero},
                               {i2, cf_zero, cf_zero, i1}}};
}

TwoQubitGateMatrix rzz_matrix(double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), -std::sin(angle / 2)};
    std::complex<double> i2 = {std::cos(angle / 2), std::sin(angle / 2)};
    return TwoQubitGateMatrix{{{i1, cf_zero, cf_zero, cf_zero},
                               {cf_zero, i2, cf_zero, cf_zero},
                               {cf_zero, cf_zero, i2, cf_zero},
                               {cf_zero, cf_zero, cf_zero, i1}}};
}

TwoQubitGateMatrix rzx_matrix(double angle) {
    std::complex<double> i1 = {std::cos(angle / 2), 0};
    std::complex<double> i2 = {0, std::sin(angle / 2)};
    return TwoQubitGateMatrix{{{i1, cf_zero, -i2, cf_zero},
                               {cf_zero, i1, cf_zero, i2},
                               {-i2, cf_zero, i1, cf_zero},
                               {cf_zero, i2, cf_zero, i1}}};
}

TwoQubitGateMatrix swap_matrix() {
    return TwoQubitGateMatrix{{{cf_one, cf_zero, cf_zero, cf_zero},
                               {cf_zero, cf_zero, cf_one, cf_zero},
                               {cf_zero, cf_one, cf_zero, cf_zero},
                               {cf_zero, cf_zero, cf_zero, cf_one}}};
}

TwoQubitGateMatrix iswap_matrix() {
    return TwoQubitGateMatrix{{{cf_one, cf_zero, cf_zero, cf_zero},
                               {cf_zero, cf_zero, cf_i, cf_zero},
                               {cf_zero, cf_i, cf_zero, cf_zero},
                               {cf_zero, cf_zero, cf_zero, cf_one}}};
}

void genDot2(vNode *node, std::vector<std::string> &result, int depth,
             std::unordered_set<vNode *> &done) {
    if (done.find(node) != done.end()) {
        return;
    }
    done.insert(node);
    std::stringstream node_ss;
    node_ss << (uint64_t)node << " [label=\"q" << depth << "\"]";
    result.push_back(node_ss.str());
    for (int i = 0; i < node->children.size(); i++) {
        std::stringstream ss;
        ss << (uint64_t)node << " -> " << (uint64_t)node->children[i].n
           << " [label=\"" << i << node->children[i].w << "\"]";
        result.push_back(ss.str());
        if (!node->children[i].isTerminal())
            genDot2(node->children[i].n, result, depth + 1, done);
    }
}

std::string genDot(vEdge &rootEdge) {
    if (rootEdge.isTerminal()) {
        return "";
    }
    std::vector<std::string> result;
    std::unordered_set<vNode *> done;
    genDot2(rootEdge.n, result, 0, done);

    // vNode::terminal
    std::stringstream node_ss;
    node_ss << (uint64_t)vNode::terminal << " [label=\"Term w=" << rootEdge.w
            << "\"]";
    result.push_back(node_ss.str());

    std::stringstream finalresult;
    finalresult << "digraph qdd {" << std::endl;
    for (std::string line : result) {
        finalresult << "  " << line << std::endl;
    }
    finalresult << "}" << std::endl;
    std::cout << done.size() << " nodes" << std::endl;
    return finalresult.str();
}

void genDot2(mNode *node, std::vector<std::string> &result, int depth,
             std::unordered_set<mNode *> &done) {
    if (done.find(node) != done.end()) {
        return;
    }
    done.insert(node);
    std::stringstream node_ss;
    node_ss << (uint64_t)node << " [label=\"q" << node->v << "\"]";
    result.push_back(node_ss.str());
    for (int i = 0; i < node->children.size(); i++) {
        if (node->children[i].w.isApproximatelyZero()) {
            continue;
        }
        std::stringstream ss;
        ss << (uint64_t)node << " -> " << (uint64_t)node->children[i].n
           << " [label=\"" << i << node->children[i].w << "\"]";
        result.push_back(ss.str());
        if (!node->children[i].isTerminal())
            genDot2(node->children[i].n, result, depth + 1, done);
    }
}

std::string genDot(mEdge &rootEdge) {
    assert(!rootEdge.isTerminal());
    std::vector<std::string> result;
    std::unordered_set<mNode *> done;
    genDot2(rootEdge.n, result, 0, done);

    // vNode::terminal
    std::stringstream node_ss;
    node_ss << (uint64_t)mNode::terminal << " [label=\"Term w=" << rootEdge.w
            << "\"]";
    result.push_back(node_ss.str());

    std::stringstream finalresult;
    finalresult << "digraph qdd {" << std::endl;
    for (std::string line : result) {
        finalresult << "  " << line << std::endl;
    }
    finalresult << "}" << std::endl;
    std::cout << done.size() << " nodes" << std::endl;
    return finalresult.str();
}

#ifdef isMPI
void save_binary(vEdge edge, std::string file_name) {
    std::vector<vContent> v;
    std::unordered_map<vNode *, int> map;
    int nNodes = vNode_to_vec(edge.n, v, map);
    std::ofstream binary_ofs(file_name);
    boost::archive::binary_oarchive binary_oa(binary_ofs);
    binary_oa << v;
    binary_ofs.close();
    std::cout << "nNodes=" << nNodes << std::endl;
}

vEdge load_binary(std::string file_name) {
    std::vector<vContent> v;
    std::ifstream binary_ifs(file_name);
    boost::archive::binary_iarchive binary_ia(binary_ifs);
    binary_ia >> v;
    binary_ifs.close();
    vEdge result = {{1.0, 0.0}, vec_to_vNode(v, vUnique)};
    return result;
}
#endif

int vNode_to_vec(vNode *node, std::vector<vContent> &table,
                 std::unordered_map<vNode *, int> &map) {
    /*
    This function is to serialize vNode* recursively.
    'table' is the outcome for serialization.
    */

    // If table is empty, always add a terminal node as id=0.
    // 'map' remember the processed vNode*, so you can avoid adding the same
    // vNode* for multiple times.
    if (table.size() == 0) {
        vContent terminal(-1, {0.0, 0.0}, {0.0, 0.0}, 0, 0);
        table.push_back(terminal);
        map[node->terminal] = 0;
    }
    if (map.find(node) != map.end()) {
        return map[node];
    }

    // If the given vNode* is not included in 'table', new data is pushed.
    int i0 = vNode_to_vec(node->children[0].getNode(), table, map);
    int i1 = vNode_to_vec(node->children[1].getNode(), table, map);
    vContent nodeData(node->v, node->children[0].w, node->children[1].w, i0,
                      i1);
    table.push_back(nodeData);
    map[node] = table.size() - 1;
    return table.size() - 1;
}

int mNode_to_vec(mNode *node, std::vector<mContent> &table,
                 std::unordered_map<mNode *, int> &map) {
    if (table.size() == 0) {
        mContent terminal(-1, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 0,
                          0, 0, 0);
        table.push_back(terminal);
        map[node->terminal] = 0;
    }
    if (map.find(node) != map.end()) {
        return map[node];
    }

    // If the given vNode* is not included in 'table', new data is pushed.
    int i0 = mNode_to_vec(node->children[0].getNode(), table, map);
    int i1 = mNode_to_vec(node->children[1].getNode(), table, map);
    int i2 = mNode_to_vec(node->children[2].getNode(), table, map);
    int i3 = mNode_to_vec(node->children[3].getNode(), table, map);
    mContent nodeData(node->v, node->children[0].w, node->children[1].w,
                      node->children[2].w, node->children[3].w, i0, i1, i2, i3);
    table.push_back(nodeData);
    map[node] = table.size() - 1;
    return table.size() - 1;
}

vNode *vec_to_vNode(std::vector<vContent> &table, vNodeTable &uniqTable,
                    bool no_lookup) {
    /*
    This function is to de-serialize table into vNode*.
    You can specify which uniqueTable to be used. (usually vUnique ?)
    */

    // The node(0) must be terminal node.
    std::unordered_map<int, vNode *> map;
    map[0] = &vNode::terminalNode;

    for (int i = 1; i < table.size(); i++) {
        vNode *node = uniqTable.getNode();
        node->v = table[i].v;
        vNode *i0 = map[table[i].index[0]];
        vNode *i1 = map[table[i].index[1]];
        vEdge e0 = {table[i].w[0], i0};
        vEdge e1 = {table[i].w[1], i1};
        node->children = {e0, e1};
        node = no_lookup ? uniqTable.register_wo_lookup(node)
                         : uniqTable.lookup(node);
        map[i] = node;
    }
    return map[table.size() - 1];
}

mNode *vec_to_mNode(std::vector<mContent> &table, mNodeTable &uniqTable,
                    bool no_lookup) {
    // The node(0) must be terminal node.
    std::unordered_map<int, mNode *> map;
    map[0] = &mNode::terminalNode;

    for (int i = 1; i < table.size(); i++) {
        mNode *node = uniqTable.getNode();
        node->v = table[i].v;
        mNode *i0 = map[table[i].index[0]];
        mNode *i1 = map[table[i].index[1]];
        mNode *i2 = map[table[i].index[2]];
        mNode *i3 = map[table[i].index[3]];
        mEdge e0 = {table[i].w[0], i0};
        mEdge e1 = {table[i].w[1], i1};
        mEdge e2 = {table[i].w[2], i2};
        mEdge e3 = {table[i].w[3], i3};
        node->children = {e0, e1, e2, e3};
        node = no_lookup ? uniqTable.register_wo_lookup(node)
                         : uniqTable.lookup(node);
        map[i] = node;
    }
    return map[table.size() - 1];
}

#ifdef isMPI
vEdge receive_dd(boost::mpi::communicator &world, int source_node_id,
                 bool isBlocking) {
    std::vector<vContent> v;
    std_complex w;
    if (isBlocking) {
        world.recv(source_node_id, 0, w);
        world.recv(source_node_id, 1, v);
    } else {
        world.irecv(source_node_id, 0, w);
        world.irecv(source_node_id, 1, v);
    }
    return {w, vec_to_vNode(v, vUnique)};
}

void send_dd(boost::mpi::communicator &world, vEdge e, int dest_node_id,
             bool isBlocking) {
    std::vector<vContent> v;
    std::unordered_map<vNode *, int> mpi_map;
    vNode_to_vec(e.n, v, mpi_map);
    if (isBlocking) {
        world.send(dest_node_id, 0, e.w);
        world.send(dest_node_id, 1, v);
    } else {
        world.isend(dest_node_id, 0, e.w);
        world.isend(dest_node_id, 1, v);
    }
    return;
};

void dump(boost::mpi::communicator &world, vEdge e, int cycle) {
    int rank = world.rank();
    std::size_t nNode;
    std::size_t send_Byte;
    std::size_t uniq_nNode;
    std::size_t uniq_Byte;

    std::vector<vContent> v;
    std::unordered_map<vNode *, int> map;
    vNode_to_vec(e.n, v, map);
    nNode = v.size();
    send_Byte = nNode * sizeof(vContent);

    uniq_nNode = vUnique.get_allocations();
    uniq_Byte = sizeof(vNode) * uniq_nNode;

    std::cout << rank << " " << cycle << " " << nNode << " " << send_Byte << " "
              << uniq_nNode << " " << uniq_Byte << std::endl;
}

double adjust_weight(bmpi::communicator &world, vEdge rootEdge) {
    const auto &[pzero, pone] = determineMeasurementProbabilities(rootEdge, 0);
    double amp = pzero + pone;
    double amp_sum = bmpi::all_reduce(world, amp, std::plus<double>());
    assert(amp_sum > 0);
    rootEdge.w = rootEdge.w / std_complex(std::sqrt(amp_sum), 0);
    return 1 / std::sqrt(amp_sum);
}
#endif

int get_nNodes(vEdge e) {
    std::vector<vContent> v;
    std::unordered_map<vNode *, int> map;
    int num = vNode_to_vec(e.n, v, map);
    return num;
}

int GC_SIZE = 1024 * 64;
vEdge gc(vEdge state, bool force) {
    if (vUnique.get_allocations() < GC_SIZE && force == false) {
        return state;
    }

    probs.clear();
    std::vector<vContent> v;
    std::unordered_map<vNode *, int> map;
    int nNodes = vNode_to_vec(state.n, v, map);
    if (nNodes > GC_SIZE) {
        GC_SIZE *= 2;
    }
    std::cout << "GC_SIZE = " << GC_SIZE << " Current nNodes = " << nNodes
              << std::endl;

    vNodeTable new_table(NQUBITS);
    vUnique = std::move(new_table);
    state.n = vec_to_vNode(v, vUnique, true);

#if defined(isMT) && !defined(CACHE_GLOBAL)
    _aCaches.clear();
    _mCaches.clear();
#else
    // AddCache newA(NQUBITS);
    // MulCache newM(NQUBITS);
    // _aCache = std::move(newA);
    // _mCache = std::move(newM);
    _mCache.clearAll();
    _aCache.clearAll();
#endif
    std::cout << "gc_done " << std::endl;
    return state;
}

int GC_SIZE_M = 1024 * 64;
mEdge gc_mat(mEdge mat, bool force) {
    if (mUnique.get_allocations() < GC_SIZE_M && force == false) {
        return mat;
    }

    std::vector<mContent> m;
    std::unordered_map<mNode *, int> map;
    int nNodes = mNode_to_vec(mat.n, m, map);
    if (nNodes > GC_SIZE_M) {
        GC_SIZE_M += nNodes;
    }
    std::cout << " Current nNodes = " << nNodes << std::endl;

    mNodeTable new_table(NQUBITS);
    mUnique = std::move(new_table);
    mat.n = vec_to_mNode(m, mUnique, true);

    // Clear identityTable
    std::vector<mEdge> new_identityTable(NQUBITS);
    identityTable = std::move(new_identityTable);

    // Clear cache
#if defined(isMT) && !defined(CACHE_GLOBAL)
    _aCaches.clear();
    _mCaches.clear();
#else
    // AddCache newA(NQUBITS);
    // MulCache newM(NQUBITS);
    // _mCache = std::move(newM);
    // _aCache = std::move(newA);
    _mCache.clearAll();
    _aCache.clearAll();
#endif
    std::cout << "gc_mat_done " << std::endl;
    return mat;
}

void set_gc_thr(int gc_v, int gc_m) {
    GC_SIZE = gc_v;
    GC_SIZE_M = gc_m;
    return;
}

int prune(vEdge &v, double thr, std_complex num) {
    std_complex current = num * v.w;
    double mag = std::sqrt(current.mag2());
    int result = 0;
    if (mag < thr) {
        v.w = {0.0, 0.0};
        v.n = vNode::terminal;
        result = 1;
    } else {
        result += prune(v.n->children[0], thr, current);
        result += prune(v.n->children[1], thr, current);
    }
    return result;
}

int prune(mEdge &m, double thr, std_complex num) {
    std_complex current = num * m.w;
    double mag = std::sqrt(current.mag2());
    int result = 0;
    if (mag < thr) {
        m.w = {0.0, 0.0};
        m.n = mNode::terminal;
        result = 1;
    } else {
        result += prune(m.n->children[0], thr, current);
        result += prune(m.n->children[1], thr, current);
        result += prune(m.n->children[2], thr, current);
        result += prune(m.n->children[3], thr, current);
    }
    return result;
}
