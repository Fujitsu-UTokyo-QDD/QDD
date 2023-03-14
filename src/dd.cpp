#include "dd.h"
#include "cache.hpp"
#include "common.h"
#include "table.hpp"
#include <algorithm>
#include <bitset>
#include <map>

#define SUBTASK_THRESHOLD 5

mNodeTable mUnique(40);
vNodeTable vUnique(40);

std::vector<mEdge> identityTable(40);

mNode mNode::terminalNode = mNode(-1, {}, nullptr, MAX_REF);
vNode vNode::terminalNode = vNode(-1, {}, nullptr, MAX_REF);

mEdge mEdge::one{.w = {1.0, 0.0}, .n = mNode::terminal, .q = -1};
mEdge mEdge::zero{.w = {0.0, 0.0}, .n = mNode::terminal, .q = -1};
vEdge vEdge::one{.w = {1.0, 0.0}, .n = vNode::terminal};
vEdge vEdge::zero{.w = {0.0, 0.0}, .n = vNode::terminal};

AddCache _aCache(40);
MulCache _mCache(40);

static int LIMIT = 10000;
const int MINUS = 3;

static mEdge normalizeM(const mEdge &e) {

    // check for all zero weights
    if (std::all_of(e.n->children.begin(), e.n->children.end(),
                    [](const mEdge &e) { return norm(e.w) == 0.0; })) {
        mUnique.returnNode(e.n);
        return mEdge::zero;
    }

    auto result = std::max_element(e.n->children.begin(), e.n->children.end(),
                                   [](const mEdge &lhs, const mEdge &rhs) {
                                       return norm(lhs.w) < norm(rhs.w);
                                   });

    std_complex max_weight = result->w;
    const std::size_t idx = std::distance(e.n->children.begin(), result);

    for (int i = 0; i < 4; i++) {
        std_complex r = e.n->children[i].w / max_weight;
        e.n->children[i].w = r;
    }

    mNode *n = mUnique.lookup(e.n);
    assert(n->v >= -1);

    return {.w = max_weight * e.w, .n = n, .q = e.q};
}

static vEdge normalizeV(const vEdge &e) {

    // check for all zero weights
    if (std::all_of(e.n->children.begin(), e.n->children.end(),
                    [](const vEdge &e) { return norm(e.w) == 0.0; })) {
        vUnique.returnNode(e.n);
        return vEdge::zero;
    }

    auto result = std::max_element(e.n->children.begin(), e.n->children.end(),
                                   [](const vEdge &lhs, const vEdge &rhs) {
                                       return norm(lhs.w) < norm(rhs.w);
                                   });

    std_complex max_weight = result->w;
    const std::size_t idx = std::distance(e.n->children.begin(), result);

    for (int i = 0; i < 2; i++) {
        std_complex r = e.n->children[i].w / max_weight;
        e.n->children[i].w = r;
    }
    vNode *n = vUnique.lookup(e.n);

    return {max_weight * e.w, n};
}

mEdge makeMEdge(Qubit q, const std::array<mEdge, 4> &c) {

    mNode *node = mUnique.getNode();
    node->v = q;
    node->children = c;

    mEdge e = normalizeM({.w = {1.0, 0.0}, .n = node, .q = q});

    assert(e.getVar() == q || e.isTerminal());

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

Qubit vEdge::getVar() const { return n->v; }

bool mEdge::isTerminal() const { return n == mNode::terminal; }

bool vEdge::isTerminal() const { return n == vNode::terminal; }

static void fillMatrix(const mEdge &edge, size_t row, size_t col,
                       const std_complex &w, uint64_t dim, std_complex **m) {

    std_complex wp = edge.w * w;

    if (edge.isTerminal()) {
        for (auto i = row; i < row + dim; i++) {
            for (auto j = col; j < col + dim; j++) {
                m[i][j] = wp;
            }
        }
        return;
    } else if (edge.isStateVector()) {
        assert(dim == edge.dim);
        for (auto i = 0; i < dim; i++) {
            for (auto j = 0; j < dim; j++) {
                m[row + i][col + j] = w * edge.mat[i][j];
            }
        }

        return;
    }

    mNode *node = edge.getNode();
    fillMatrix(node->getEdge(0), row, col, wp, dim / 2, m);
    fillMatrix(node->getEdge(1), row, col + dim / 2, wp, dim / 2, m);
    fillMatrix(node->getEdge(2), row + dim / 2, col, wp, dim / 2, m);
    fillMatrix(node->getEdge(3), row + dim / 2, col + dim / 2, wp, dim / 2, m);
}

void mEdge::printMatrix() const {
    if (this->isTerminal()) {
        std::cout << this->w << std::endl;
        return;
    }
    Qubit q = this->getVar();
    std::size_t dim = 1 << (q + 1);

    std_complex **matrix = new std_complex *[dim];
    for (std::size_t i = 0; i < dim; i++)
        matrix[i] = new std_complex[dim];

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

std_complex **mEdge::getMatrix(std::size_t *dim) const {
    assert(!this->isTerminal());

    Qubit q = this->getVar();
    std::size_t d = 1 << (q + 1);

    std_complex **matrix = new std_complex *[d];
    for (std::size_t i = 0; i < d; i++)
        matrix[i] = new std_complex[d];

    fillMatrix(*this, 0, 0, {1.0, 0.0}, d, matrix);
    if (dim != nullptr)
        *dim = d;
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

MatrixXcf mEdge::getEigenMatrix() {
    std::size_t dim;
    auto m = getMatrix(&dim);
    MatrixGuard g(m, dim);
    MatrixXcf M(dim, dim);

    for (auto i = 0; i < dim; i++) {
        for (auto j = 0; j < dim; j++) {
            M(i, j) = std::complex<double>{m[i][j].r, m[i][j].i};
        }
    }
    return M;
}

bool mEdge::compareNumerically(const mEdge &other) const noexcept {
    if (this->getVar() != other.getVar())
        return false;

    auto m1 = this->getMatrix(nullptr);
    auto m2 = other.getMatrix(nullptr);
    Qubit q = this->getVar();
    std::size_t dim = 1 << (q + 1);
    MatrixGuard mg1(m1, dim);
    MatrixGuard mg2(m2, dim);

    for (auto i = 0; i < dim; i++) {
        for (auto j = 0; j < dim; j++) {
            if (m1[i][j] != m2[i][j])
                return false;
        }
    }
    return true;
}

mEdge makeIdent(Qubit q) {

    if (q < 0)
        return mEdge::one;

    if (identityTable[q].n != nullptr) {
        assert(identityTable[q].n->v > -1);
        return identityTable[q];
    }

    mEdge e = makeMEdge(0, {mEdge::one, mEdge::zero, mEdge::zero, mEdge::one});
    for (Qubit i = 1; i <= q; i++) {
        e = makeMEdge(i, {{e, mEdge::zero, mEdge::zero, e}});
    }

    identityTable[q] = e;
    // e.incRef();
    return e;
}

vEdge makeZeroState(QubitCount q) {
    vEdge e = makeVEdge(0, {vEdge::one, vEdge::zero});
    for (Qubit i = 1; i < q; i++) {
        e = makeVEdge(i, {{e, vEdge::zero}});
    }
    return e;
}

vEdge makeOneState(QubitCount q) {
    vEdge e = makeVEdge(0, {vEdge::zero, vEdge::one});
    for (Qubit i = 1; i < q; i++) {
        e = makeVEdge(i, {{vEdge::zero, e}});
    }
    return e;
}
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
        if (it != c.end() && it->qubit == z)
            ++it;
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

static void set_4_diagonal_submatrice_with(Complex **mat, size_t row,
                                           size_t col, size_t dim,
                                           const GateMatrix &gate) {

    assert(row == col);

    const size_t half = dim / 2;
    // top half
    for (size_t i = row; i < row + half; i++) {
        mat[i][i] = gate[0];
        mat[i][i + half] = gate[1];
    }
    // bottom half
    for (size_t i = row + half; i < row + dim; i++) {
        mat[i][i - half] = gate[2];
        mat[i][i] = gate[3];
    }
}

//     q's < threshold are represented using matrix
//     e.g., threshold = 4, then q0, q1,  q2 and q3,are represented by a 16x16
//     matrix
mEdge make_hybrid_with_small_target(QubitCount q, GateMatrix g, Qubit target,
                                    Qubit threshold) {

    // We currently only support  target below threshold
    assert(threshold < q && target < threshold);
    size_t dim = 1 << threshold;
    Complex **mat = new Complex *[dim];
    for (auto i = 0; i < dim; i++)
        mat[i] = new Complex[dim];

    // suppose target = q2, then this is a 8x8 matrix
    size_t target_dim = 1 << (target + 1);

    // set the diagonal submatrices of size target_dim x target_dim
    for (auto i = 0; i < dim; i += target_dim) {
        set_4_diagonal_submatrice_with(mat, i, i, target_dim, g);
    }

    // the matrix is ready! begin to build the tree..
    mEdge e;
    e.w = {1.0, 0.0};
    e.mat = mat;
    e.dim = dim;
    e.q = threshold - 1;

    for (Qubit z = threshold; z < q; z++) {
        e = makeMEdge(z, {e, mEdge::zero, mEdge::zero, e});
    }

    return e;
}

static void set_diagonal(Complex **mat, size_t dim, const Complex &c) {
    for (auto i = 0; i < dim; i++)
        mat[i][i] = c;
}
// target >= threshold
mEdge make_hybrid_with_large_target(QubitCount q, GateMatrix g, Qubit target,
                                    Qubit threshold) {

    assert(threshold < q && target >= threshold);
    size_t dim = 1 << threshold;

    std::array<Complex **, 4> mats;

    for (auto i = 0; i < 4; i++) {
        mats[i] = new Complex *[dim];
        for (auto j = 0; j < dim; j++)
            mats[i][j] = new Complex[dim];

        set_diagonal(mats[i], dim, g[i]);
    }

    mEdge e;
    std::array<mEdge, 4> edges;
    if (threshold == target) {
        for (auto i = 0; i < 4; i++) {
            mEdge tmp;
            tmp.w = {1.0, 0.0};
            tmp.mat = mats[i];
            tmp.dim = dim;
            tmp.q = threshold - 1;
            edges[i] = tmp;
        }

    } else {
        for (auto i = 0; i < 4; i++) {
            mEdge tmp;
            tmp.w = {1.0, 0.0};
            tmp.mat = mats[i];
            tmp.dim = dim;
            tmp.q = threshold - 1;

            edges[i] =
                makeMEdge(threshold, {tmp, mEdge::zero, mEdge::zero, tmp});
        }

        for (auto z = threshold + 1; z < target; z++) {
            for (auto i = 0; i < 4; i++) {
                edges[i] = makeMEdge(
                    z, {edges[i], mEdge::zero, mEdge::zero, edges[i]});
            }
        }
    }

    e = makeMEdge(target, edges);

    for (auto z = target + 1; z < q; z++) {
        e = makeMEdge(z, {e, mEdge::zero, mEdge::zero, e});
    }

    return e;
}
mEdge makeHybridGate(QubitCount q, GateMatrix g, Qubit target,
                     Qubit threshold) {
    assert(target < q);

    if (target < threshold)
        return make_hybrid_with_small_target(q, g, target, threshold);
    else
        return make_hybrid_with_large_target(q, g, target, threshold);
}

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

    mEdge result;

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
            x = lhs;
        }
        if (rv == current_var && !rhs.isTerminal()) {
            y = rnode->getEdge(i);
            y.w = rhs.w * y.w;
        } else {
            y = rhs;
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
                x = lcopy;
            }

            if (rv == current_var && !rhs.isTerminal()) {
                y = rnode->getEdge((k << 1) | col);
            } else {
                y = rcopy;
            }

            product[k] = mm_multiply2(x, y, current_var - 1);
        }
        edges[i] = mm_add2(product[0], product[1], current_var - 1);
    }

    result = makeMEdge(current_var, edges);
    _mCache.set(lhs.n, rhs.n, result);

    result.w = result.w * lhs.w * rhs.w;
    if (result.w.isApproximatelyZero())
        return mEdge::zero;
    else
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
        return rhs;
    } else if (rhs.w.isApproximatelyZero()) {
        return lhs;
    }

    if (current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, vNode::terminal};
    }

    vEdge result;

    result = _aCache.find(lhs, rhs);
    if (result.n != nullptr) {
        return result;
    }

    vEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    vNode *lnode = lhs.getNode();
    vNode *rnode = rhs.getNode();
    std::array<vEdge, 2> edges;

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

        edges[i] = vv_add2(x, y, current_var - 1);
    }

    result = makeVEdge(current_var, edges);
    _aCache.set(lhs, rhs, result);

    return result;
}

vEdge vv_add(const vEdge &lhs, const vEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w + rhs.w, vNode::terminal};
    }

    // assume lhs and rhs are the same size vector.
    assert(lhs.getVar() == rhs.getVar());
    return vv_add2(lhs, rhs, lhs.getVar());
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

static void printVector_sparse2(const vEdge &edge, std::size_t row,
                                const std_complex &w, uint64_t left,
                                std::map<int, std_complex> &m) {

    std_complex wp = edge.w * w;

    const double thr = 0.001;

    if (edge.isTerminal() && left == 0) {
        if (norm(wp) > thr)
            m[row] = wp;
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
    if (dim != nullptr)
        *dim = d;
    return vector;
}
struct VectorGuard {
    VectorGuard(std_complex *m, std::size_t dim) : _m(m), _dim(dim) {}
    ~VectorGuard() { delete[] _m; }

    std_complex *_m;
    const std::size_t _dim;
};

VectorXcf vEdge::getEigenVector() {
    std::size_t dim;
    auto v = getVector(&dim);
    VectorGuard g(v, dim);
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
        std::array<vEdge, 2> product;
        for (auto k = 0; k < 2; k++) {
            if (lv == current_var && !lhs.isTerminal()) {
                x = lnode->getEdge((i << 1) | k);
            } else {
                x = lcopy;
            }

            if (rv == current_var && !rhs.isTerminal()) {
                y = rnode->getEdge(k);
            } else {
                y = rcopy;
            }

            product[k] = mv_multiply2(x, y, current_var - 1);
        }

        edges[i] = vv_add2(product[0], product[1], current_var - 1);
    }

    result = makeVEdge(current_var, edges);
    _mCache.set(lhs.n, rhs.n, result);
    result.w = result.w * lhs.w * rhs.w;
    if (result.w.isApproximatelyZero()) {
        return vEdge::zero;
    } else {
        return result;
    }
}

struct mEdgeRefGuard {
    mEdgeRefGuard(mEdge &ee) : e(ee) { e.incRef(); }
    ~mEdgeRefGuard() { e.decRef(); };
    mEdge e;
};
struct vEdgeRefGuard {
    vEdgeRefGuard(vEdge &ee) : e(ee) { e.incRef(); }
    ~vEdgeRefGuard() { e.decRef(); };
    vEdge e;
};

vEdge mv_multiply(mEdge lhs, vEdge rhs) {

    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {lhs.w * rhs.w, vNode::terminal};
    }

    // assume lhs and rhs are the same length.
    assert(lhs.getVar() == rhs.getVar());
    vEdge v = mv_multiply2(lhs, rhs, lhs.getVar());
    return v;
}

std::string measureAll(vEdge &rootEdge, const bool collapse,
                       std::mt19937_64 &mt, double epsilon) {
    if (std::abs(rootEdge.w.mag2() - 1.0L) > epsilon) {
        if (rootEdge.w.isApproximatelyZero()) {
            throw std::runtime_error("led to a 0-vector");
        }
    }

    vEdge cur = rootEdge;
    const auto nqubits = static_cast<QubitCount>(rootEdge.getVar() + 1);
    std::string result(nqubits, '0');
    std::uniform_real_distribution<double> dist(0.0, 1.0L);

    for (Qubit i = rootEdge.getVar(); i >= 0; --i) {
        double p0 = cur.n->getEdge(0).w.mag2();
        double p1 = cur.n->getEdge(1).w.mag2();
        double tmp = p0 + p1;

        if (std::abs(tmp - 1.0L) > epsilon) {
            throw std::runtime_error("Added probabilities differ from 1 by " +
                                     std::to_string(std::abs(tmp - 1.0L)));
        }

        p0 /= tmp;

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

mEdge makeSwap(QubitCount q, Qubit target0, Qubit target1) {
    Controls c1{Control{target0, Control::Type::pos}};
    mEdge e1 = makeGate(q, Xmat, target1, c1);

    Controls c2{Control{target1, Control::Type::pos}};
    mEdge e2 = makeGate(q, Xmat, target0, c2);

    mEdge e3 = mm_multiply(e2, e1);
    e3 = mm_multiply(e1, e3);

    return e3;
}

void mEdge::decRef() {
    if (isTerminal())
        return;
    assert(n->ref > 0);
    n->ref--;

    if (n->ref == 0) {
        for (auto i = 0; i < 4; i++)
            n->getEdge(i).decRef();
    }
}
void vEdge::decRef() {
    if (isTerminal())
        return;

    assert(n->ref > 0);
    n->ref--;

    if (n->ref == 0) {
        for (auto i = 0; i < 2; i++)
            n->getEdge(i).decRef();
    }
}
void mEdge::incRef() {
    if (isTerminal())
        return;

    if (n->ref < MAX_REF) {
        n->ref++;
    }

    if (n->ref == 1) {
        for (auto i = 0; i < 4; i++)
            n->getEdge(i).incRef();
    }
}
void vEdge::incRef() {
    if (isTerminal())
        return;

    if (n->ref < MAX_REF) {
        n->ref++;
    }

    if (n->ref == 1) {
        for (auto i = 0; i < 2; i++)
            n->getEdge(i).incRef();
    }
}
void mEdge::check() {
    if (isTerminal())
        return;

    if (n->v == -2 || n->v == -3 || n->v == -4) {
        assert(false);
    }
    for (auto i = 0; i < 4; i++)
        n->getEdge(i).check();
}

mEdge RX(QubitCount qnum, int target, float angle) {
    std::complex<float> i1 = {std::cos(angle / 2), 0};
    std::complex<float> i2 = {0, -std::sin(angle / 2)};
    return makeGate(qnum, GateMatrix{i1, i2, i2, i1}, target);
}
mEdge RY(QubitCount qnum, int target, float angle) {
    std::complex<float> i1 = {std::cos(angle / 2), 0};
    std::complex<float> i2 = {-std::sin(angle / 2), 0};
    std::complex<float> i3 = {std::sin(angle / 2), 0};
    return makeGate(qnum, GateMatrix{i1, i2, i3, i1}, target);
}
mEdge RZ(QubitCount qnum, int target, float angle) {
    std::complex<float> i1 = {std::cos(angle / 2), -std::sin(angle / 2)};
    std::complex<float> i2 = {std::cos(angle / 2), std::sin(angle / 2)};
    return makeGate(qnum, GateMatrix{i1, cf_zero, cf_zero, i2}, target);
}

mEdge CX(QubitCount qnum, int target, int control) {
    std::complex<float> zero = {0, 0};
    std::complex<float> one = {1, 0};
    Controls controls;
    controls.emplace(Control{control, Control::Type::pos});
    return makeGate(qnum, GateMatrix{zero, one, one, zero}, target, controls);
}

mEdge Dense(QubitCount qnum, int target, float r0, float i0, float r1, float i1,
            float r2, float i2, float r3, float i3) {
    std::complex<float> c0 = {r0, i0};
    std::complex<float> c1 = {r1, i1};
    std::complex<float> c2 = {r2, i2};
    std::complex<float> c3 = {r3, i3};
    return makeGate(qnum, GateMatrix{c0, c1, c2, c3}, target);
}

mEdge Dense1211(QubitCount qnum, int target, int control0, int control1,
                float r0, float i0, float r1, float i1, float r2, float i2,
                float r3, float i3) {
    std::complex<float> c0 = {r0, i0};
    std::complex<float> c1 = {r1, i1};
    std::complex<float> c2 = {r2, i2};
    std::complex<float> c3 = {r3, i3};
    Controls controls;
    controls.emplace(Control{control0, Control::Type::pos});
    controls.emplace(Control{control1, Control::Type::pos});
    return makeGate(qnum, GateMatrix{c0, c1, c2, c3}, target, controls);
}
