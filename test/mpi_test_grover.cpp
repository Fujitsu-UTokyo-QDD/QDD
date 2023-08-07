#include "table.hpp"
#include "cache.hpp"
#include<iostream>
#include<fstream>

static unsigned long long CalculateIterations(const unsigned short n_qubits) {
    constexpr long double PI_4 =
        0.785398163397448309615660845819875721049292349843776455243L;
    if (n_qubits <= 3) {
        return 1;
    } else {
        return static_cast<unsigned long long>(
            std::floor(PI_4 * std::pow(2.L, n_qubits / 2.0L)));
    }
}

mEdge buildUnitary(const std::vector<mEdge> &g) {

    if (g.size() == 0) {
        return mEdge();
    }

    mEdge lhs = g[0];
    for (int i = 1; i < g.size(); i++) {
        // std::cout<<"i: "<<i<<std::endl;
        lhs = mm_multiply(lhs, g[i]);
    }
    return lhs;
}

static mEdge groverIteration(const std::string &oracle, QubitCount n_qubits) {

    std::vector<mEdge> g;
    QubitCount total_qubits = n_qubits + 1;

    // prepare oracle
    Controls controls;
    for (auto i = 0; i < n_qubits; i++) {
        controls.emplace(Control{i, oracle.at(i) == '1' ? Control::Type::pos
                                                        : Control::Type::neg});
    }

    mEdge o = makeGate(total_qubits, Zmat, n_qubits, controls);
    g.push_back(o);
    //    o.printMatrix();

    // std::cout<<"ocale"<<std::endl;
    // o.printMatrix();

    // prepare diffusioin
    mEdge d = makeIdent(n_qubits);
    g.push_back(d);
    // 1. H to data qubits
    for (auto i = 0; i < n_qubits; i++) {
        g.emplace_back(makeGate(total_qubits, Hmat, i));
    }
    // 2. X to data qubits
    for (auto i = 0; i < n_qubits; i++) {
        g.emplace_back(makeGate(total_qubits, Xmat, i));
    }

    // 3. H to the last data qubit
    g.emplace_back(makeGate(total_qubits, Hmat, n_qubits - 1));

    // 4. CX to the last data qubit
    Controls diff_controls;
    for (auto i = 0; i < n_qubits - 1; i++) {
        diff_controls.emplace(Control{i, Control::Type::pos});
    }
    g.emplace_back(makeGate(total_qubits, Xmat, n_qubits - 1, diff_controls));

    // 5. H to the last data qubit
    g.emplace_back(makeGate(total_qubits, Hmat, n_qubits - 1));

    // 6. X to data qubits
    for (auto i = 0; i < n_qubits; i++) {
        g.emplace_back(makeGate(total_qubits, Xmat, i));
    }

    // 7. H to all qubits

    for (auto i = 0; i < n_qubits; i++) {
        g.emplace_back(makeGate(total_qubits, Hmat, i));
    }

    return buildUnitary(g);
}


int vNode_to_vec2(vNode *node, std::vector<vContent> &table,
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
    int i0 = vNode_to_vec2(node->children[0].getNode(), table, map);
    int i1 = vNode_to_vec2(node->children[1].getNode(), table, map);
    vContent nodeData(node->v, node->children[0].w, node->children[1].w, i0,
                      i1);
    table.push_back(nodeData);
    map[node] = table.size() - 1;
    return table.size() - 1;
}

vNode *vec_to_vNode2(std::vector<vContent> &table, vNodeTable &uniqTable) {
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
        node = uniqTable.lookup(node);
        map[i] = node;
    }
    return map[table.size() - 1];
}

vEdge grover_MPI(QubitCount n_qubits, bmpi::communicator &world) {
    auto t1 = std::chrono::high_resolution_clock::now();
    std::size_t iterations = CalculateIterations(n_qubits);
    std::mt19937_64 mt;
    std::array<std::mt19937_64::result_type, std::mt19937_64::state_size>
        random_data{};
    std::random_device rd;
    std::generate(std::begin(random_data), std::end(random_data),
                  [&rd]() { return rd(); });
    std::seed_seq seeds(std::begin(random_data), std::end(random_data));
    mt.seed(100);
    // Generate random oracle
    std::uniform_int_distribution<int> dist(0, 1); // range is inclusive
    std::string oracle = std::string(n_qubits, '0');
    for (Qubit i = 0; i < n_qubits; i++) {
        if (dist(mt) == 1) {
            oracle[i] = '1';
        }
    }
    if(world.rank()==0)
        std::cout << "orcale: " << oracle << std::endl;

    QubitCount total_qubits = n_qubits + 1;

    mEdge full_iteration = groverIteration(oracle, n_qubits);

    // set it up
    std::vector<mEdge> setup_gates;
    setup_gates.emplace_back(makeGate(total_qubits, Xmat, n_qubits));
    // apply H to all qubits
    for (auto i = 0; i < n_qubits; i++) {
        setup_gates.emplace_back(makeGate(total_qubits, Hmat, i));
    }
    mEdge setup = buildUnitary(setup_gates);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()/1000000<<" seconds"<<std::endl;
    world.barrier();

    auto t3 = std::chrono::high_resolution_clock::now();
    vEdge state = makeZeroStateMPI(total_qubits, world);
    state = mv_multiply_MPI(setup, state, world);

    unsigned int j_pre = 0;
    if(world.rank()==0)
        std::cout << "iterations: " << iterations << std::endl;

    while ((iterations - j_pre) % 8 != 0){
        state = mv_multiply_MPI(full_iteration, state, world);
        j_pre++;
    }
    int i = 0;
    for (unsigned long long j = j_pre; j < iterations; j += 8) {
        state = mv_multiply_MPI(full_iteration, state, world);
        state = mv_multiply_MPI(full_iteration, state, world);
        state = mv_multiply_MPI(full_iteration, state, world);
        state = mv_multiply_MPI(full_iteration, state, world);
        state = mv_multiply_MPI(full_iteration, state, world);
        state = mv_multiply_MPI(full_iteration, state, world);
        state = mv_multiply_MPI(full_iteration, state, world);
        state = mv_multiply_MPI(full_iteration, state, world);
        if(j%100000<8 && j>100000){
            dump(world, state, j);
            std::cout << world.rank() << " gc" << std::endl;
            std::vector<vContent> v;
            std::unordered_map<vNode *, int> map;
            vNode_to_vec2(state.n, v, map);

            vNodeTable new_table(NQUBITS);
            vUnique = std::move(new_table);
            state.n = vec_to_vNode2(v, vUnique);

            AddCache newA(NQUBITS);
            MulCache newM(NQUBITS);
            _aCache = std::move(newA);
            _mCache = std::move(newM);
            dump(world, state, j);
            std::cout << "Weight adjustment x" << adjust_weight(world, state) << std::endl;
        }
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    ms = t4 - t3;
    std::cout<<ms.count()/1000000<<" seconds"<<std::endl;
    return state;
}

int main(int argc, char** argv){
    bmpi::environment env(argc, argv);
    bmpi::communicator world;
    assert(argc == 2);
    auto result = grover_MPI(std::atoi(argv[1]), world);
    std::cout << std::endl << genDot(result);
}
