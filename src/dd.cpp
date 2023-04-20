#include "dd.h"
#include "cache.hpp"
#include "common.h"
#include "complexnumbers.hpp"
#include "table.hpp"
#include <algorithm>
#include <bitset>
#include <map>
#include <queue>

#define SUBTASK_THRESHOLD 5

mNodeTable mUnique(40);
vNodeTable vUnique(40);

std::vector<mEdge> identityTable(40);

mNode mNode::terminalNode = mNode(-1, {}, nullptr, MAX_REF);
vNode vNode::terminalNode = vNode(-1, {}, nullptr, MAX_REF);

mEdge mEdge::one{.w = Complex::one, .n = mNode::terminal, .q = -1};
mEdge mEdge::zero{.w = Complex::zero, .n = mNode::terminal, .q = -1};
vEdge vEdge::one{.w = Complex::one, .n = vNode::terminal};
vEdge vEdge::zero{.w = Complex::zero, .n = vNode::terminal};

AddCache _aCache(40);
MulCache _mCache(40);

static int LIMIT = 10000;
const int MINUS = 3;

ComplexNumbers cn{};

using CTEntry = ComplexTable<>::Entry;

static void normalizeM(mEdge &e) {

    // check for all zero weights

    auto zero = std::array{e.n->children[0].w.approximatelyZero(),
                           e.n->children[1].w.approximatelyZero(),
                           e.n->children[2].w.approximatelyZero(),
                           e.n->children[3].w.approximatelyZero()};

    if (zero[0] && zero[1] && zero[2] && zero[3]) {
        mUnique.returnNode(e.n);
        e = mEdge::zero;
        return;
    }

    auto result = std::max_element(e.n->children.begin(), e.n->children.end(),
                                   [](const mEdge &lhs, const mEdge &rhs) {
                                       return ComplexNumbers::mag2(lhs.w) <
                                              ComplexNumbers::mag2(rhs.w);
                                   });

    Complex max_weight = result->w;
    const std::size_t idx = std::distance(e.n->children.begin(), result);

    for (int i = 0; i < 4; i++) {
        if (i == idx) {
            if (e.w.exactlyOne()) {
                e.w = max_weight;
            } else {
                ComplexNumbers::mul(e.w, e.w, max_weight);
            }
            e.n->children[i].w = Complex::one;
        } else {
            if (zero[i]) {
                if (e.n->children[i].w != Complex::zero) {
                    cn.returnToCache(e.n->children[i].w);
                }

                e.n->children[i] = mEdge::zero;
                continue;
            }
            /*
            if (!zero[i] && !e.n->children[i].w.exactlyOne()) {
                cn.returnToCache(e.n->children[i].w);
            }
            */

            if (e.n->children[i].w.approximatelyOne()) {
                e.n->children[i].w = Complex::one;
            }

            auto c = cn.getTemporary();
            ComplexNumbers::div(c, e.n->children[i].w, max_weight);
            e.n->children[i].w = cn.lookup(c);
        }
    }

    mNode *n = mUnique.lookup(e.n);
    e.n = n;
    assert(n->v >= -1);
}

static void normalizeV(vEdge &e) {
    auto zero = std::array{e.n->children[0].w.approximatelyZero(),
                           e.n->children[1].w.approximatelyZero()};

    if (zero[0] && zero[1]) {
        vUnique.returnNode(e.n);
        e = vEdge::zero;
        return;
    }

    auto tmp = std::sqrt(ComplexNumbers::mag2(e.n->children[0].w) + ComplexNumbers::mag2(e.n->children[1].w));
    auto sum = cn.lookup(tmp, 0.0);

    if(!sum.approximatelyOne() && !sum.approximatelyZero()){
        for (int i = 0; i < 2; i++) {
            auto c = cn.getTemporary();
            ComplexNumbers::div(c, e.n->children[i].w, sum);
            e.n->children[i].w = cn.lookup(c);
        }
    }
    auto c = cn.getTemporary();
    ComplexNumbers::mul(c, e.w, sum);
    e.w = cn.lookup(c);

    vNode *n = vUnique.lookup(e.n);
    e.n = n;
    assert(n->v >= -1);
}

mEdge makeMEdge(Qubit q, const std::array<mEdge, 4> &c) {

    mNode *node = mUnique.getNode();
    node->v = q;
    node->children = c;

    mEdge e = mEdge{.w = Complex::one, .n = node, .q = q};

    normalizeM(e);

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

    vEdge e = vEdge{.w = Complex::one, .n = node};

    normalizeV(e);

    assert(e.getVar() == q || e.isTerminal());

    return e;
}

Qubit vEdge::getVar() const { return n->v; }

bool mEdge::isTerminal() const { return n == mNode::terminal; }

bool vEdge::isTerminal() const { return n == vNode::terminal; }

static void fillMatrix(const mEdge &edge, size_t row, size_t col,
                       const std_complex &w, uint64_t dim, std_complex **m) {

    std_complex wp =
        std_complex{CTEntry::val(edge.w.r), CTEntry::val(edge.w.i)} * w;

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

std::complex<double> **mEdge::getMatrix(std::size_t *dim) const {
    assert(!this->isTerminal());

    Qubit q = this->getVar();
    std::size_t d = 1 << (q + 1);

    std::complex<double> **matrix = new std::complex<double> *[d];
    for (std::size_t i = 0; i < d; i++)
        matrix[i] = new std::complex<double>[d];

    fillMatrix(*this, 0, 0, {1.0, 0.0}, d, matrix);
    if (dim != nullptr)
        *dim = d;
    return matrix;
}

bool mEdge::compareNumerically(const mEdge &other) const noexcept {
    if (this->getVar() != other.getVar())
        return false;

    auto m1 = this->getMatrix(nullptr);
    auto m2 = other.getMatrix(nullptr);
    Qubit q = this->getVar();
    std::size_t dim = 1 << (q + 1);

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
        edges[i] = mEdge{cn.lookup(g[i].real(), g[i].imag()), mNode::terminal};

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

static Qubit rootVar(const mEdge &lhs, const mEdge &rhs) {
    assert(!(lhs.isTerminal() && rhs.isTerminal()));

    return (lhs.isTerminal() ||
            (!rhs.isTerminal() && (rhs.getVar() > lhs.getVar())))
               ? rhs.getVar()
               : lhs.getVar();
}

mEdge mm_add2(const mEdge &lhs, const mEdge &rhs, int32_t current_var) {
    if (lhs.w.approximatelyZero()) {
        auto r = rhs;
        r.w = cn.getCached(CTEntry::val(rhs.w.r), CTEntry::val(rhs.w.i));
        return r;
    } else if (rhs.w.approximatelyZero()) {
        auto r = lhs;
        r.w = cn.getCached(CTEntry::val(lhs.w.r), CTEntry::val(lhs.w.i));
        return r;
    }

    if (current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {cn.addCached(lhs.w, rhs.w), mNode::terminal};
    }

    mEdge result;

    result = _aCache.find(lhs, rhs);
    if (result.n != nullptr) {
        if (result.w.approximatelyZero()) {
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
            x.w = cn.mulCached(lhs.w, x.w);
        } else {
            x = lhs;
        }
        if (rv == current_var && !rhs.isTerminal()) {
            y = rnode->getEdge(i);
            y.w = cn.mulCached(rhs.w, y.w);
        } else {
            y = rhs;
        }

        edges[i] = mm_add2(x, y, current_var - 1);
        if (!x.isTerminal() && x.n->v == current_var && x.w != Complex::zero) {
            cn.returnToCache(x.w);
        }
        if (!y.isTerminal() && y.n->v == current_var && y.w != Complex::zero) {
            cn.returnToCache(y.w);
        }
    }

    result = makeMEdge(current_var, edges);
    _aCache.set(lhs, rhs, result);

    return result;
}

mEdge mm_add(const mEdge &lhs, const mEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {cn.addCached(lhs.w, rhs.w), mNode::terminal};
    }

    Qubit root = rootVar(lhs, rhs);
    return mm_add2(lhs, rhs, root);
}

mEdge mm_multiply2(const mEdge &lhs, const mEdge &rhs, int32_t current_var) {

    if (lhs.w.approximatelyZero() || rhs.w.approximatelyZero()) {
        return mEdge::zero;
    }

    if (current_var == -1) {

        assert(lhs.isTerminal() && rhs.isTerminal());
        return {cn.mulCached(lhs.w, rhs.w), mNode::terminal};
    }

    mEdge result;
    result = _mCache.find(lhs.n, rhs.n);
    if (result.n != nullptr) {
        if (result.w.approximatelyZero()) {
            return mEdge::zero;
        } else {
            result = mEdge{.w = cn.getCached(CTEntry::val(result.w.r),
                                             CTEntry::val(result.w.i)),
                           .n = result.n};
            ComplexNumbers::mul(result.w, result.w, lhs.w);
            ComplexNumbers::mul(result.w, result.w, rhs.w);
            if (result.w.approximatelyZero()) {
                cn.returnToCache(result.w);
                return mEdge::zero;
            }

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
    lcopy.w = Complex::one;
    rcopy.w = Complex::one;

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

    ComplexNumbers::mul(result.w, result.w, lhs.w);
    ComplexNumbers::mul(result.w, result.w, rhs.w);

    if (result.w.approximatelyZero())
        return mEdge::zero;
    else
        return result;
}

mEdge mm_multiply(const mEdge &lhs, const mEdge &rhs) {

    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {cn.mulCached(lhs.w, rhs.w), mNode::terminal};
    }

    Qubit root = rootVar(lhs, rhs);
    mEdge result = mm_multiply2(lhs, rhs, root);
    return result;
}

mEdge mm_kronecker2(const mEdge &lhs, const mEdge &rhs) {
    if (lhs.isTerminal()) {
        return {cn.mulCached(lhs.w, rhs.w), rhs.n};
    }

    std::array<mEdge, 4> edges;
    mEdge x;
    mNode *lnode = lhs.getNode();

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    for (auto i = 0; i < 4; i++) {
        x = lnode->getEdge(i);
        ComplexNumbers::mul(x.w, x.w, lhs.w);

        edges[i] = mm_kronecker2(x, rhs);
    }

    mEdge ret = makeMEdge(lv + rv + 1, edges);
    return ret;
}

mEdge mm_kronecker(const mEdge &lhs, const mEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {cn.mulCached(lhs.w, rhs.w), mNode::terminal};
    }

    return mm_kronecker2(lhs, rhs);
}

vEdge vv_add2(const vEdge &lhs, const vEdge &rhs, int32_t current_var) {
    if (lhs.w.approximatelyZero()) {
        return rhs;
    } else if (rhs.w.approximatelyZero()) {
        return lhs;
    }

    if (current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {cn.addCached(lhs.w, rhs.w), vNode::terminal};
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
            x.w = cn.mulCached(x.w, lhs.w);
        } else {
            x = lhs;
        }
        if (rv == current_var && !rhs.isTerminal()) {
            y = rnode->getEdge(i);
            y.w = cn.mulCached(y.w, rhs.w);
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
        return {cn.addCached(lhs.w, rhs.w), vNode::terminal};
    }

    // assume lhs and rhs are the same size vector.
    assert(lhs.getVar() == rhs.getVar());
    return vv_add2(lhs, rhs, lhs.getVar());
}

vEdge vv_kronecker2(const vEdge &lhs, const vEdge &rhs) {
    if (lhs.isTerminal()) {
        return {cn.mulCached(lhs.w, rhs.w), rhs.n};
    }
    if (rhs.isTerminal()) {
        return {cn.mulCached(lhs.w, rhs.w), lhs.n};
    }

    std::array<vEdge, 2> edges;
    vEdge x;
    vNode *lnode = lhs.getNode();

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    for (auto i = 0; i < 2; i++) {
        x = lnode->getEdge(i);
        ComplexNumbers::mul(x.w, x.w, lhs.w);
        edges[i] = vv_kronecker2(x, rhs);
    }

    vEdge ret = makeVEdge(lv + rv + 1, edges);
    return ret;
}

vEdge vv_kronecker(const vEdge &lhs, const vEdge &rhs) {
    if (lhs.isTerminal() && rhs.isTerminal()) {
        return {cn.mulCached(lhs.w, rhs.w), vNode::terminal};
    }

    return vv_kronecker2(lhs, rhs);
}
static void fillVector(const vEdge &edge, std::size_t row, const std_complex &w,
                       uint64_t dim, std_complex *m) {

    std_complex wp =
        std_complex{CTEntry::val(edge.w.r), CTEntry::val(edge.w.i)} * w;
    if (edge.isTerminal()) {
        for (auto i = row; i < row + dim; i++) {
            m[i] = wp;
        }
        return;
    }

    vNode *node = edge.getNode();
    fillVector(node->getEdge(0), row, wp, dim / 2, m);
    fillVector(node->getEdge(1), row + dim / 2, wp, dim / 2, m);
}

void vEdge::printVector() const {
    if (this->isTerminal()) {
        std::cout << this->w << std::endl;
        return;
    }
    Qubit q = this->getVar();
    std::size_t dim = 1 << (q + 1);

    std_complex *vector = new std_complex[dim];

    fillVector(*this, 0, {1.0, 0.0}, dim, vector);

    for (size_t i = 0; i < dim; i++) {
        std::cout << vector[i] << " ";
        std::cout << "\n";
    }
    std::cout << std::endl;

    delete[] vector;
}

std_complex *vEdge::getVector(std::size_t *dim) const {
    assert(!this->isTerminal());

    Qubit q = this->getVar();
    std::size_t d = 1 << (q + 1);

    std_complex *vector = new std_complex[d];
    fillVector(*this, 0, {1.0, 0.0}, d, vector);
    if (dim != nullptr)
        *dim = d;
    return vector;
}

vEdge mv_multiply2(const mEdge &lhs, const vEdge &rhs, int32_t current_var) {

    if (lhs.w.approximatelyZero() || rhs.w.approximatelyZero()) {
        return vEdge::zero;
    }

    if (current_var == -1) {

        assert(lhs.isTerminal() && rhs.isTerminal());
        return {cn.mulCached(lhs.w, rhs.w), vNode::terminal};
    }

    vEdge result;

    result = _mCache.find(lhs.n, rhs.n);
    if (result.n != nullptr) {
        if (result.w.approximatelyZero()) {
            return vEdge::zero;
        } else {
            if(lhs.w.approximatelyOne()){
                return result;
            }
            if(result.w.approximatelyOne()){
                result.w = lhs.w;
                return result;
            }

            ComplexNumbers::mul(result.w, result.w, lhs.w);
            ComplexNumbers::mul(result.w, result.w, rhs.w);
            if (result.w.approximatelyZero())
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
    lcopy.w = Complex::one;
    rcopy.w = Complex::one;

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
    result.w = cn.mulCached(result.w, lhs.w);
    result.w = cn.mulCached(result.w, rhs.w);
    if (result.w.approximatelyZero()) {
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
        return {cn.mulCached(lhs.w, rhs.w), vNode::terminal};
    }

    // assume lhs and rhs are the same length.
    // assert(lhs.getVar() == rhs.getVar());
    vEdge v = mv_multiply2(lhs, rhs, lhs.getVar());
    v.w = Complex::one;
    return v;
}

std::string measureAll(vEdge &rootEdge, const bool collapse,
                       std::mt19937_64 &mt, double epsilon) {
    if (std::abs(ComplexNumbers::mag2(rootEdge.w) - 1.0L) > epsilon) {
        if (rootEdge.w.approximatelyZero()) {
            throw std::runtime_error("led to a 0-vector");
        }
    }

    vEdge cur = rootEdge;
    const auto nqubits = static_cast<QubitCount>(rootEdge.getVar() + 1);
    std::string result(nqubits, '0');
    std::uniform_real_distribution<double> dist(0.0, 1.0L);
    for (Qubit i = rootEdge.getVar(); i >= 0; --i)
    {
        double p0 = ComplexNumbers::mag2(cur.n->getEdge(0).w);
        double p1 = ComplexNumbers::mag2(cur.n->getEdge(1).w);
        double tmp = p0 + p1;

        //if (std::abs(tmp - 1.0L) > epsilon) {
        //    throw std::runtime_error("Added probabilities differ from 1 by " +
        //                             std::to_string(std::abs(tmp - 1.0L)));
        //}
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

double assignProbabilities(const vEdge &edge,
                           std::unordered_map<vNode *, double> &probs) {
    auto it = probs.find(edge.n);
    if (it != probs.end()) {
        return ComplexNumbers::mag2(edge.w) * it->second;
    }
    double sum{1};
    if (!edge.isTerminal()) {
        sum = assignProbabilities(edge.n->children.at(0), probs) +
              assignProbabilities(edge.n->children.at(1), probs);
    }

    probs.insert({edge.n, sum});

    return ComplexNumbers::mag2(edge.w) * sum;
}

std::pair<double, double>
determineMeasurementProbabilities(const vEdge &rootEdge, const Qubit index,
                                  const bool assumeProbabilityNormalization) {
    std::map<vNode *, double> probsMone;
    std::set<vNode *> visited;
    std::queue<vNode *> q;

    probsMone[rootEdge.n] = ComplexNumbers::mag2(rootEdge.w);
    visited.insert(rootEdge.n);
    q.push(rootEdge.n);

    while (q.front()->v != index) {
        vNode *ptr = q.front();
        q.pop();
        const double prob = probsMone[ptr];

        if (!ptr->children.at(0).w.approximatelyZero()) {
            const double tmp1 =
                prob * ComplexNumbers::mag2(ptr->children.at(0).w);

            if (visited.find(ptr->children.at(0).n) != visited.end()) {
                probsMone[ptr->children.at(0).n] =
                    probsMone[ptr->children.at(0).n] + tmp1;
            } else {
                probsMone[ptr->children.at(0).n] = tmp1;
                visited.insert(ptr->children.at(0).n);
                q.push(ptr->children.at(0).n);
            }
        }

        if (!ptr->children.at(1).w.approximatelyZero()) {
            const double tmp1 =
                prob * ComplexNumbers::mag2(ptr->children.at(1).w);

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

    if (assumeProbabilityNormalization) {
        while (!q.empty()) {
            vNode *ptr = q.front();
            q.pop();

            if (!ptr->children.at(0).w.approximatelyZero()) {
                pzero += probsMone[ptr] *
                         ComplexNumbers::mag2(ptr->children.at(0).w);
            }

            if (!ptr->children.at(1).w.approximatelyZero()) {
                pone += probsMone[ptr] *
                        ComplexNumbers::mag2(ptr->children.at(1).w);
            }
        }
    } else {
        std::unordered_map<vNode *, double> probs;
        assignProbabilities(rootEdge, probs);

        while (!q.empty()) {
            vNode *ptr = q.front();
            q.pop();

            if (!ptr->children.at(0).w.approximatelyZero()) {
                pzero += probsMone[ptr] * probs[ptr->children.at(0).n] *
                         ComplexNumbers::mag2(ptr->children.at(0).w);
            }

            if (!ptr->children.at(1).w.approximatelyZero()) {
                pone += probsMone[ptr] * probs[ptr->children.at(1).n] *
                        ComplexNumbers::mag2(ptr->children.at(1).w);
            }
        }
    }
    return {pzero, pone};
}

char measureOneCollapsing(vEdge &rootEdge, const Qubit index,
                          const bool assumeProbabilityNormalization,
                          std::mt19937_64 &mt, double epsilon) {
    const auto &[pzero, pone] = determineMeasurementProbabilities(
        rootEdge, index, assumeProbabilityNormalization);
    const double sum = pzero + pone;
    // if (std::abs(sum - 1) > epsilon) {
    //     throw std::runtime_error(
    //         "Numerical instability occurred during measurement: |alpha|^2 "
    //         "+ "
    //         "|beta|^2 = " +
    //         std::to_string(pzero) + " + " + std::to_string(pone) + " = " +
    //         std::to_string(pzero + pone) + ", but should be 1!");
    // }
    GateMatrix measurementMatrix{cf_zero, cf_zero, cf_zero, cf_zero};

    std::uniform_real_distribution<double> dist(0.0, 1.0L);

    double threshold = dist(mt);
    double normalizationFactor; // NOLINT(cppcoreguidelines-init-variables)
                                // always assigned a value in the following
                                // block
    char result; // NOLINT(cppcoreguidelines-init-variables) always assigned
                 // a value in the following block

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
    c = std_complex{CTEntry::val(e.w.r), CTEntry::val(e.w.i)} * c;
    e.w = cn.getCached(c.real(), c.imag());
    rootEdge = e;

    return result;
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

void genDot2(vNode *node, std::vector<std::string> &result, int depth) {
    std::stringstream node_ss;
    node_ss << (uint64_t)node << " [label=\"q" << depth <<"(n=" << node->v << ")\"]";
    result.push_back(node_ss.str());
    for (int i = 0; i < node->children.size(); i++) {
        std::stringstream ss;
        ss << (uint64_t)node << " -> " << (uint64_t)node->children[i].n
           << " [label=\"" << i <<"(w="<< node->children[i].w << ")\"]";
        result.push_back(ss.str());
        if (!node->children[i].isTerminal())
            genDot2(node->children[i].n, result, depth + 1);
    }
}

std::string genDot(vEdge &rootEdge) {
    std::vector<std::string> result;
    result.push_back("0 [label=\"start\"]");
    std::stringstream edge_ss;
    edge_ss << 0 << " -> " << (uint64_t)rootEdge.n
           << " [label=\"(w="<< rootEdge.w << ")\"]";
    result.push_back(edge_ss.str());
    genDot2(rootEdge.n, result, 0);

    // vNode::terminal
    std::stringstream node_ss;
    node_ss << (uint64_t)vNode::terminal << " [label=\"Term\"]";
    result.push_back(node_ss.str());

    std::stringstream finalresult;
    finalresult << "digraph qdd {" << std::endl;
    for (std::string line : result) {
        finalresult << "  " << line << std::endl;
    }
    finalresult << "}" << std::endl;
    return finalresult.str();
}
