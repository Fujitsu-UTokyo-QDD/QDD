#include "algorithms/grover.hpp"

static unsigned long long CalculateIterations(const unsigned short n_qubits) {
    constexpr long double PI_4 = 0.785398163397448309615660845819875721049292349843776455243L; // dd::PI_4 is of type fp and hence possibly smaller than long double
    if (n_qubits <= 3) {
        return 1;
    } else {
        return static_cast<unsigned long long>(std::floor(PI_4 * std::pow(2.L, n_qubits / 2.0L)));
    }
}


mEdge Grover::makeFullIteration(){
    
    QuantumCircuit qc(n_qubits + n_anciallae, _nworkers, _reduce);
    QubitCount total_qubits = n_qubits + n_anciallae;

    //prepare oracle
    Controls controls;
    for(auto i = 0; i < n_qubits; i++){
        controls.emplace(Control{i, oracle.at(i) == '1'? Control::Type::pos: Control::Type::neg});
    }
    qc.emplace_back(Zmat, n_qubits, controls);

     

    //prepare diffusioin
    // 1. H to data qubits
    for(auto i = 0 ; i < n_qubits; i++){
        qc.emplace_back(Hmat, i);
    }

    //2. X to data qubits
    for(auto i = 0; i < n_qubits; i++){
        qc.emplace_back(Xmat,i);   
    }


    //3. H to the last data qubit
    qc.emplace_back(Hmat, n_qubits - 1);


    //4. CX to the last data qubit
    Controls diff_controls;
    for(auto i = 0; i < n_qubits - 1; i++){
        diff_controls.emplace(Control{i, Control::Type::pos});
    }
    qc.emplace_back(Xmat, n_qubits -1 , diff_controls);


    //5. H to the last data qubit
    qc.emplace_back(Hmat, n_qubits - 1);

    //6. X to data qubits
    for(auto i = 0; i < n_qubits; i++){
        qc.emplace_back(Xmat,i);   
    }


    //7. H to all qubits
    for(auto i = 0 ; i < total_qubits; i++){
        qc.emplace_back(Hmat, i);
    }


    qc.buildCircuit();

    return qc.wait().matrixResult();
    

}

Grover::Grover(QubitCount q, int workers, int reduce, std::size_t seed): n_qubits(q), _nworkers(workers), _reduce(reduce), seed(seed){
    
    iterations = CalculateIterations(n_qubits);
    std::cout<<iterations<<" iterations"<<std::endl;
    std::array<std::mt19937_64::result_type, std::mt19937_64::state_size> random_data{};
    std::random_device                                                    rd;
    std::generate(std::begin(random_data), std::end(random_data), [&rd]() { return rd(); });
    std::seed_seq seeds(std::begin(random_data), std::end(random_data));
    mt.seed(seeds);
    //Generate random oracle
    std::uniform_int_distribution<int> dist(0, 1); // range is inclusive
    oracle = std::string(n_qubits, '0');
    for (Qubit i = 0; i < n_qubits; i++) {
        if (dist(mt) == 1) {
            oracle[i] = '1';
        }
    }
    std::cout<<"orcale: "<<oracle<<std::endl;

}




void Grover::full_grover(){
    auto t1 = std::chrono::high_resolution_clock::now();
    mEdge full_iteration = makeFullIteration();

    QuantumCircuit qc(n_qubits + n_anciallae, _nworkers, _reduce);
    QubitCount total_qubits = n_qubits + n_anciallae;

    //create init state  and set it up
    qc.setInput(makeZeroState(total_qubits));
    //prepare the oracle qubit to be 1
    qc.emplace_back(Xmat, n_qubits);
    //apply H to all qubits
    for(auto i = 0; i < total_qubits; i++) qc.emplace_back(Hmat, i);


    for(auto j = 0;j < iterations; j++){
        qc.emplace_gate(full_iteration);
    
    }
    
    qc.buildCircuit();
    vEdge result = qc.wait().vectorResult();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" micro s"<<std::endl;
    //auto s = measureAll(result, false, mt, 0.001L);
    //std::cout<<s<<std::endl;

}
