#pragma once
#include "common.h"
#include "dd.h"
#include <vector>
#include <random>
#include <algorithm>
#include "task.h"


class Shor{
    public:
        Shor(int composite_number, int w, int gcfreq,  bool v = false ):
            n(composite_number), coprime_a(2), _nworkers(w), verbose(v), s(w, gcfreq), 
            required_bits(std::ceil(std::log2(composite_number))), n_qubits(2 * required_bits + 3),  approximate(false) {
                if(n%2 == 0){
                    std::cout<<"only support factorizing odd numbers"<<std::endl;
                    exit(1);
                }
                std::cout<<"total qubits: "<<n_qubits<<std::endl;
                std::array<std::mt19937_64::result_type, std::mt19937_64::state_size> random_data{};
                std::random_device                                                    rd;
                std::generate(std::begin(random_data), std::end(random_data), [&rd]() { return rd(); });
                std::seed_seq seeds(std::begin(random_data), std::end(random_data));
                mt.seed(seeds);

            };

        static unsigned long long modpow(unsigned long long base, unsigned long long exp, unsigned long long modulus) {
            base %= modulus;
            unsigned long long result = 1ull;
            while (exp > 0) {
                if (exp & 1ull) result = (result * base) % modulus;
                base = (base * base) % modulus;
                exp >>= 1ull;
            }
            return result;
        }

        static unsigned long long gcd(unsigned long long a, unsigned long long b) {
            unsigned long long c;
            while (a != 0) {
                c = a;
                a = b % a;
                b = c;
            }
            return b;
        }

        static float QDDcos(float fac, float div) {
            return std::cos((PI * fac) / div);
        }

        static float QDDsin(float fac, float div) {
            return std::sin((PI * fac) / div);
        }


        void cmult_inv(int a, int N, int c);

        void cmult( int a, int N, int c);

        void mod_add_phi_inv( int a, int N, int c1, int c2);

        void mod_add_phi( int a, int N, int c1, int c2);

        void qft_inv();

        void qft();

        void add_phi_inv(int a, int c1, int c2);

        void add_phi( int a, int c1, int c2);

        static int inverse_mod(int a, int n);

        void u_a( unsigned long long a, int N, int c) ;


        std::vector<unsigned long long> ts;

        mEdge addConst(unsigned long long a);

        mEdge addConstMod(unsigned long long a);

        mEdge limitTo(unsigned long long a);

        std::pair<unsigned int, unsigned int> post_processing(const std::string& sample) const;

        void run();
    private:
        /// composite number to be factored
        const unsigned int n;
        /// coprime number to `n`. Setting this to zero will randomly generate a suitable number
        unsigned int       coprime_a;
        const unsigned int required_bits;
        QubitCount     n_qubits{};

        std::string                   sim_result = "did not start";
        std::pair<unsigned, unsigned> sim_factors{0, 0};

        std::string                   polr_result = "did not start";
        std::pair<unsigned, unsigned> polr_factors{0, 0};

        const bool verbose;

        const bool         approximate;
        unsigned long long approximation_runs{0};
        long double        final_fidelity{1.0L};
        double             step_fidelity{0.9};
        std::mt19937_64 mt;

        int _nworkers;


        Scheduler s;



};
