#include "algorithms/shor.hpp"
#include "graph.hpp"

mEdge Shor::u_a_emulate(unsigned long long a, int q) {
    QuantumCircuit qc(n_qubits , _nworkers, _reduce);
    Worker w(n_qubits);

    mEdge limit = makeIdent(required_bits - 1);

    mEdge  f = mEdge::one;
    std::array<mEdge, 4> edges{
            mEdge::zero,
            mEdge::zero,
            mEdge::zero,
            mEdge::zero};

    for (unsigned int p = 0; p < required_bits; ++p) {
        edges[0] = f;
        edges[1] = f;
        f        = makeMEdge(p, edges);
    }

    f = mm_multiply(&w, f, limit);    

    edges[1] = mEdge::zero;


    unsigned long t = a;

    for (unsigned int i = 0; i < required_bits; ++i) {
        mEdge active = mEdge::one;
        for (unsigned int p = 0; p < required_bits; ++p) {
            if (p == i) {
                edges[3] = active;
                edges[0] = mEdge::zero;
            } else {
                edges[0] = edges[3] = active;
            }
            active = makeMEdge(p, edges);
        }

        active.w          = {-1.0,0.0};
        qc.emplace_back(mm_add(&w, limit, active));
        qc.emplace_back(f);
        active.w          = {1.0,0.0};
        qc.emplace_back(active);
        qc.emplace_back(f);
// here
        dd::mEdge tmp = addConstMod(t);
        active        = dd->multiply(tmp, active);

        dd->decRef(f);
        f = dd->add(active, passive);
        dd->incRef(f);
        dd->garbageCollect();

        t = (2 * t) % n;
    }

    dd->decRef(limit);
    dd->decRef(f);

    dd::mEdge e = f;

    for (int i = 2 * required_bits - 1; i >= 0; --i) {
        if (i == q) {
            edges[1] = edges[2] = dd::mEdge::zero;
            edges[0]            = dd->makeIdent(0, n_qubits - i - 2);
            edges[3]            = e;
            e                   = dd->makeDDNode(n_qubits - 1 - i, edges, false);
        } else {
            edges[1] = edges[2] = dd::mEdge::zero;
            edges[0] = edges[3] = e;
            e                   = dd->makeDDNode(n_qubits - 1 - i, edges, false);
        }
    }

    dd::vEdge tmp = dd->multiply(e, rootEdge);
    dd->incRef(tmp);
    dd->decRef(rootEdge);
    rootEdge = tmp;

    dd->garbageCollect();
}
void Shor::run(){

    QubitCount n_qubits; 
    vEdge rootEdge;

    n_qubits = 3 * required_bits;
    QuantumCircuit qc(n_qubits , _nworkers, _reduce);
    qc.setInput(makeZeroState(n_qubits));

    qc.emplace_back(Xmat, 0);

    if(coprime_a == 0){
        std::uniform_int_distribution<unsigned int> distribution(1, n - 1); // range is inclusive
        do {
            coprime_a = distribution(mt);
        } while (gcd(coprime_a, n) != 1 || coprime_a == 1);
        
    }

    auto* as                  = new unsigned long long[2 * required_bits];
    as[2 * required_bits - 1] = coprime_a;
    unsigned long long new_a  = coprime_a;
    for (int i = 2 * required_bits - 2; i >= 0; i--) {
        new_a = new_a * new_a;
        new_a = new_a % n;
        as[i] = new_a;
    }

    for (unsigned int i = 0; i < 2 * required_bits; i++) {
        qc.emplace_back(Hmat, (n_qubits-1) - i);
    }
    const int mod = std::ceil(2 * required_bits / 6.0); // log_0.9(0.5) is about 6
    for (unsigned int i = 0; i < 2 * required_bits; i++) {
        u_a_emulate(as[i], i);
    }
                                                        


}
