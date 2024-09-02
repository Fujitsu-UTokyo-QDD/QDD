#include "algorithms/grover.hpp"

#include "dd.h"
#include "table.hpp"

static unsigned long long CalculateIterations(const unsigned short n_qubits) {
    constexpr long double PI_4 =
        0.785398163397448309615660845819875721049292349843776455243L;  // dd::PI_4
                                                                       // is of
                                                                       // type
                                                                       // fp and
                                                                       // hence
                                                                       // possibly
                                                                       // smaller
                                                                       // than
                                                                       // long
                                                                       // double
    if (n_qubits <= 3) {
        return 1;
    } else {
        return static_cast<unsigned long long>(
            std::floor(PI_4 * std::pow(2.L, n_qubits / 2.0L)));
    }
}

static mEdge groverIteration(Scheduler &s, const std::string &oracle,
                             QubitCount n_qubits) {
    std::vector<mEdge> g;
    QubitCount total_qubits = n_qubits + 1;

    // prepare oracle
    Controls controls;
    for (auto i = 0; i < n_qubits; i++) {
        controls.emplace(Control{
            i, oracle.at(i) == '1' ? Control::Type::pos : Control::Type::neg});
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

    // std::cout<<"diff"<<std::endl;
    // d.printMatrix();

    mEdge full = s.buildUnitary(g);
    // std::cout<<"full"<<std::endl;
    // full.printMatrix();
    // std::cout<<"full ite"<<std::endl;
    // full.printMatrix();
    return full;
}

vEdge grover(Scheduler &s, QubitCount n_qubits) {
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
    std::uniform_int_distribution<int> dist(0, 1);  // range is inclusive
    std::string oracle = std::string(n_qubits, '0');
    for (Qubit i = 0; i < n_qubits; i++) {
        if (dist(mt) == 1) {
            oracle[i] = '1';
        }
    }
    // oracle = "011";
    std::cout << "orcale: " << oracle << std::endl;

    QubitCount total_qubits = n_qubits + 1;

    mEdge full_iteration = groverIteration(s, oracle, n_qubits);
    // full_iteration.incRef();
    std::cout << "after full: " << std::endl;

    // set it up
    std::vector<mEdge> setup_gates;
    setup_gates.emplace_back(makeGate(total_qubits, Xmat, n_qubits));
    // apply H to all qubits
    for (auto i = 0; i < n_qubits; i++) {
        setup_gates.emplace_back(makeGate(total_qubits, Hmat, i));
    }
    mEdge setup = s.buildUnitary(setup_gates);
    // std::cout<<"setup"<<std::endl;
    // setup.printMatrix();
    s.addGate(setup);

    unsigned int j_pre = 0;
    std::cout << "iterations: " << iterations << std::endl;

    while ((iterations - j_pre) % 8 != 0) {
        s.addGate(full_iteration);
        j_pre++;
    }
    int i = 0;
    for (unsigned long long j = j_pre; j < iterations; j += 8) {
        // std::cout<<"ite: "<<j_pre<<std::endl;
        s.addGate(full_iteration);
        s.addGate(full_iteration);
        s.addGate(full_iteration);
        s.addGate(full_iteration);
        s.addGate(full_iteration);
        s.addGate(full_iteration);
        s.addGate(full_iteration);
        s.addGate(full_iteration);
    }

    vEdge result = s.buildCircuit(makeZeroState(total_qubits));

    return result;
}
