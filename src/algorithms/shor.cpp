#include "algorithms/shor.hpp"
#include "graph.hpp"

int Shor::inverse_mod(int a, int n) {
    int t    = 0;
    int newt = 1;
    int r    = n;
    int newr = a;
    while (newr != 0) {
        int quotient = r / newr;
        int h        = t;
        t            = newt;
        newt         = h - quotient * newt;
        h            = r;
        r            = newr;
        newr         = h - quotient * newr;
    }
    if (r > 1) {
        std::cerr << "ERROR: a=" << a << " with n=" << n << " is not invertible\n";
        std::exit(3);
    }
    if (t < 0) {
        t = t + n;
    }
    return t;
}

void Shor::add_phi(int a, int c1, int c2) {
    for (int i = required_bits; i >= 0; --i) {
        double       q   = 1;
        unsigned int fac = 0;
        for (int j = i; j >= 0; --j) {
            if ((a >> j) & 1u) {
                fac |= 1u;
            }
            fac *= 2;
            q *= 2;
        }

        Controls controls;
        if (c1 != std::numeric_limits<int>::min()) {
            controls.emplace(Control{static_cast<Qubit>((n_qubits - 1) - c1)});
        }
        if (c2 != std::numeric_limits<int>::min()) {
            controls.emplace(Control{static_cast<Qubit>((n_qubits - 1) - c2)});
        }

        float         q_r = QDDcos(fac, q);
        float         q_i = QDDsin(fac, q);
        GateMatrix Qm{cf_one, cf_zero, cf_zero, {q_r, q_i}};

        qc.emplace_back(Qm, (n_qubits - 1)-(1+2*required_bits - i), controls);
    }
}

void Shor::add_phi_inv(int a, int c1, int c2) {
    for (int i = required_bits; i >= 0; --i) {
        double       q   = 1;
        unsigned int fac = 0;
        for (int j = i; j >= 0; --j) {
            if ((a >> j) & 1u) {
                fac |= 1u;
            }
            fac *= 2;
            q *= 2;
        }
        Controls controls;
        if (c1 != std::numeric_limits<int>::min()) {
            controls.emplace(Control{static_cast<Qubit>((n_qubits - 1) - c1)});
        }
        if (c2 != std::numeric_limits<int>::min()) {
            controls.emplace(Control{static_cast<Qubit>((n_qubits - 1) - c2)});
        }

        float         q_r = QDDcos(fac, -q);
        float         q_i = QDDsin(fac, -q);
        GateMatrix Qm{cf_one, cf_zero, cf_zero, {q_r, q_i}};
        qc.emplace_back(Qm, (n_qubits - 1)-(1+2*required_bits - i), controls);
    }
}


void Shor::mod_add_phi(int a, int N, int c1, int c2) {
    add_phi(a, c1, c2);
    add_phi_inv(N, std::numeric_limits<int>::min(), std::numeric_limits<int>::min());

    qft_inv();

    qc.emplace_back(Xmat, (n_qubits-1)-(2*required_bits +2), Controls{Control{(Qubit)((n_qubits - 1)-(required_bits + 1)), Control::Type::pos}});

    qft();
    add_phi(N, 2 * required_bits + 2, std::numeric_limits<int>::min());
    add_phi_inv(a, c1, c2);
    qft_inv();



    qc.emplace_back(Xmat, (n_qubits-1)-(2*required_bits +2), Controls{Control{(Qubit)((n_qubits - 1)-(required_bits + 1)), Control::Type::neg}});
    qft();
    add_phi(a, c1, c2);
}

void Shor::mod_add_phi_inv(int a, int N, int c1, int c2) {
    add_phi_inv(a, c1, c2);
    qft_inv();

    qc.emplace_back(Xmat, (n_qubits-1)-(2*required_bits +2), Controls{Control{(Qubit)((n_qubits - 1)-(required_bits + 1)), Control::Type::neg}});

    qft();
    add_phi(a, c1, c2);
    add_phi_inv(N, 2 * required_bits + 2, std::numeric_limits<int>::min());
    qft_inv();


    qc.emplace_back(Xmat, (n_qubits-1)-(2*required_bits +2), Controls{Control{(Qubit)((n_qubits - 1)-(required_bits + 1)), Control::Type::pos}});

    qft();
    add_phi(N, std::numeric_limits<int>::min(), std::numeric_limits<int>::min());
    add_phi_inv(a, c1, c2);
}

//[start, end)
static void qft_rotations(QuantumCircuit& qc, Qubit start, Qubit end){
    if(start == end) 
        return;

    end -= 1;
    qc.emplace_back(Hmat, end);

    for(Qubit q = start; q < end; q++){
        Controls c{Control{q, Control::Type::pos}};
        float r = std::cos(std::numbers::pi/(std::pow(2, end - q))); 
        float i = std::sin(std::numbers::pi/(std::pow(2, end - q))); 
        GateMatrix g{cf_one, cf_zero, cf_zero, {r,i}}; 
        qc.emplace_back(g, end, c);
    }

    qft_rotations(qc, start, end);
    
}


//[start, end)
static void qft_swap(QuantumCircuit& qc, Qubit start, Qubit end){
    QubitCount total_qubits = qc.getQubits();
    Qubit e = end-1;
    for(Qubit q = start; q < (start+end)/2; q++){
        qc.emplace_gate(makeSwap(total_qubits, q, e-- ));
    }
}

static void full_qft(QuantumCircuit& qc, Qubit start, Qubit end){
    qft_rotations(qc, start, end);
    qft_swap(qc, start, end);
}

void Shor::qft() {
    for (unsigned int i = required_bits + 1; i < 2 * required_bits + 2; i++) {
        qc.emplace_back(Hmat, (n_qubits-1) - i);

        double q = 2;
        for (unsigned int j = i + 1; j < 2 * required_bits + 2; j++) {
            float         q_r = QDDcos(1, q);
            float         q_i = QDDsin(1, q);
            GateMatrix Qm{cf_one, cf_zero, cf_zero, {q_r, q_i}};
            Controls c;
            c.emplace(Control{(Qubit)((n_qubits - 1)- j), Control::Type::pos});
            qc.emplace_back(Qm, (n_qubits - 1) - i, c );
            q *= 2;
        }
    }
}

void Shor::qft_inv() {
    for (unsigned int i = 2 * required_bits + 1; i >= required_bits + 1; i--) {
        double q = 2;
        for (unsigned int j = i + 1; j < 2 * required_bits + 2; j++) {
            float         q_r = QDDcos(1, -q);
            float         q_i = QDDsin(1, -q);
            GateMatrix Qm{cf_one, cf_zero, cf_zero, {q_r, q_i}};
            Controls c;
            c.emplace(Control{(Qubit)((n_qubits - 1)- j), Control::Type::pos});
            qc.emplace_back(Qm, (n_qubits - 1) - i, c );
            q *= 2;
        }
        qc.emplace_back(Hmat, (n_qubits- 1)- i);
    }
}
void Shor::cmult(int a, int N, int c) {
    qft();

    int t = a;
    for (int i = required_bits; i >= 1; i--) {
        mod_add_phi(t, N, i, c);
        t = (2 * t) % N;
    }
    qft_inv();
}

void Shor::cmult_inv(int a, int N, int c) {
    qft();
    int t = a;
    for (int i = required_bits; i >= 1; i--) {
        mod_add_phi_inv(t, N, i, c);
        t = (2 * t) % N;
    }
    qft_inv();
}
void Shor::u_a(unsigned long long a, int N, int c) {
    cmult(a, N, c);
    for (unsigned int i = 0; i < required_bits; i++) {
        qc.emplace_back(Xmat, (n_qubits-1)-(i + 1), Controls{Control{(Qubit)((n_qubits - 1)-(required_bits + 2 + i)), Control::Type::pos}});

        qc.emplace_back(Xmat, (n_qubits-1)-(required_bits + 2 + i), Controls{Control{(Qubit)((n_qubits - 1)-(i + 1)), Control::Type::pos}, Control{(Qubit)(n_qubits-1-c), Control::Type::pos}});

        qc.emplace_back(Xmat, (n_qubits-1)-(i + 1), Controls{Control{(Qubit)((n_qubits - 1)-(required_bits + 2 + i)), Control::Type::pos}});
    }

    cmult_inv(inverse_mod(a, N), N, c);
}

void Shor::run(){

    auto t1 = std::chrono::high_resolution_clock::now();
    vEdge rootEdge;

    qc.setInput(makeZeroState(n_qubits));

    qc.emplace_back(Xmat, n_qubits - 1);

    std::cout<<"coprime_a = "<< coprime_a<<std::endl;

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
        u_a(as[i], n, 0);
    }

    for (unsigned int i = 0; i < 2 * required_bits; i++) {

        double q = 2;

        for (int j = i - 1; j >= 0; j--) {
            float         q_r = QDDcos(1, -q);
            float         q_i = QDDsin(1, -q);
            GateMatrix Qm{cf_one, cf_zero, cf_zero, {q_r, q_i}};
            qc.emplace_back(Qm, n_qubits -1 - i, Controls{Control{(Qubit)(n_qubits- 1- j), Control::Type::pos}});
            q *= 2;
        }

        qc.emplace_back(Hmat, n_qubits - 1 - i);
    }
    qc.buildCircuit();
    qc.dump_task_graph();

    vEdge result = qc.wait().vectorResult();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" micro s"<<std::endl;
    /*
    {
        std::string sample_reversed = measureAll(result, false, mt, 0.001L);
        std::string sample{sample_reversed.rbegin(), sample_reversed.rend()};
        sim_factors = post_processing(sample);
        if (sim_factors.first != 0 && sim_factors.second != 0) {
            sim_result =
                    std::string("SUCCESS(") + std::to_string(sim_factors.first) + "*" +
                    std::to_string(sim_factors.second) + ")";
        } else {
            sim_result = "FAILURE";
        }
    }
*/
    delete[] as;
                                                        
}
std::pair<unsigned int, unsigned int> Shor::post_processing(const std::string& sample) const {
    unsigned long long res = 0;
    if (verbose) {
        std::clog << "measurement: ";
    }
    for (unsigned int i = 0; i < 2 * required_bits; i++) {
        if (verbose) {
            std::clog << sample.at(required_bits + i);
        }
        res = (res << 1u) + (sample.at(required_bits + i) == '1');
    }

    if (verbose) {
        std::clog << " = " << res << "\n";
    }
    unsigned long long denom = 1ull << (2 * required_bits);

    bool               success = false;
    unsigned long long f1{0}, f2{0};
    if (res == 0) {
        if (verbose) {
            std::clog << "Factorization failed (measured 0)!" << std::endl;
        }
    } else {
        if (verbose) {
            std::clog << "Continued fraction expansion of " << res << "/" << denom << " = " << std::flush;
        }
        std::vector<unsigned long long> cf;

        unsigned long long old_res   = res;
        unsigned long long old_denom = denom;
        while (res != 0) {
            cf.push_back(denom / res);
            unsigned long long tmp = denom % res;
            denom                  = res;
            res                    = tmp;
        }

        if (verbose) {
            for (const auto i: cf) {
                std::clog << i << " ";
            }
            std::clog << "\n";
        }

        for (unsigned int i = 0; i < cf.size(); i++) {
            //determine candidate
            unsigned long long denominator = cf[i];
            unsigned long long numerator   = 1;

            for (int j = i - 1; j >= 0; j--) {
                unsigned long long tmp = numerator + cf[j] * denominator;
                numerator              = denominator;
                denominator            = tmp;
            }
            if (verbose) {
                std::clog << "  Candidate " << numerator << "/" << denominator << ": ";
            }
            if (denominator > n) {
                if (verbose) {
                    std::clog << " denominator too large (greater than " << n << ")!\n";
                }
                success = false;
                if (verbose) {
                    std::clog << "Factorization failed!\n";
                }
                break;
            } else {
                double delta = (double)old_res / (double)old_denom - (double)numerator / (double)denominator;
                if (std::abs(delta) < 1.0 / (2.0 * old_denom)) {
                    unsigned long long fact = 1;
                    while (denominator * fact < n && modpow(coprime_a, denominator * fact, n) != 1) {
                        fact++;
                    }
                    if (modpow(coprime_a, denominator * fact, n) == 1) {
                        if (verbose) {
                            std::clog << "found period: " << denominator << " * " << fact << " = "
                                      << (denominator * fact) << "\n";
                        }
                        if ((denominator * fact) & 1u) {
                            if (verbose) {
                                std::clog << "Factorization failed (period is odd)!\n";
                            }
                        } else {
                            f1 = modpow(coprime_a, (denominator * fact) >> 1u, n);
                            f2 = (f1 + 1) % n;
                            f1 = (f1 == 0) ? n - 1 : f1 - 1;
                            f1 = gcd(f1, n);
                            f2 = gcd(f2, n);

                            if (f1 == 1ull || f2 == 1ull) {
                                if (verbose) {
                                    std::clog << "Factorization failed: found trivial factors " << f1 << " and " << f2
                                              << "\n";
                                }
                            } else {
                                if (verbose) {
                                    std::clog << "Factorization succeeded! Non-trivial factors are: \n"
                                              << "  -- gcd(" << n << "^(" << (denominator * fact) << "/2)-1"
                                              << "," << n
                                              << ") = " << f1 << "\n"
                                              << "  -- gcd(" << n << "^(" << (denominator * fact) << "/2)+1"
                                              << "," << n
                                              << ") = " << f2 << "\n";
                                }
                                success = true;
                            }
                        }

                        break;
                    } else {
                        if (verbose) {
                            std::clog << "failed\n";
                        }
                    }
                } else {
                    if (verbose) {
                        std::clog << "delta is too big (" << delta << ")\n";
                    }
                }
            }
        }
    }
    if (success) {
        return {f1, f2};
    } else {
        return {0, 0};
    }
}
