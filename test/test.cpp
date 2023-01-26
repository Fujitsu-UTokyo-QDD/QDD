#include <iostream>
#include "common.h"
#include "dd.h"
#include <chrono>
#include <cstdlib>
#include <random>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "gtest/gtest.h"
#include "graph.hpp"



using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using Eigen::MatrixXcf;
using Eigen::Matrix2cf;
using Eigen::VectorXcf;


#define EIGEN_MATRIX(GATE) \
    static Matrix2cf GATE = (Matrix2cf()<<::GATE[0], ::GATE[1], ::GATE[2], ::GATE[3]).finished();

namespace qdd_test{

    EIGEN_MATRIX(Imat)
    EIGEN_MATRIX(Hmat)
    EIGEN_MATRIX(Xmat)
    EIGEN_MATRIX(Ymat)
    EIGEN_MATRIX(Zmat)
    EIGEN_MATRIX(Smat)
    EIGEN_MATRIX(Sdagmat)
    EIGEN_MATRIX(Tmat)
    EIGEN_MATRIX(SXmat)
    EIGEN_MATRIX(SXdagmat)
    EIGEN_MATRIX(Vmat)
    EIGEN_MATRIX(Vdagmat)

}

/*

mEdge run(){
    QuantumCircuit qc(10,8);
    qc.emplace_back(Hmat, 1);
    qc.emplace_back(Xmat, 2);
    qc.emplace_back(Ymat, 3);
    qc.emplace_back(Vdagmat, 4);
    qc.emplace_back(Hmat, 5);
    qc.emplace_back(Sdagmat, 6);
    qc.emplace_back(Hmat, 7);
    qc.emplace_back(Xmat, 8);
    qc.emplace_back(Ymat, 9);
    qc.emplace_back(Vdagmat, 1);
    qc.emplace_back(Hmat, 2);
    qc.emplace_back(Sdagmat, 3);
    qc.buildCircuit();
//    std::cout<<"Task graph: "<<std::endl;
//    qc.dump_task_graph();
    auto t1 = std::chrono::high_resolution_clock::now();
    mEdge result = qc.wait().matrixResult();
    auto t2 = std::chrono::high_resolution_clock::now();
    duration<double, std::micro> ms = t2 - t1;
    std::cout<<ms.count()<<" micro s"<<std::endl;
    return result;

}


int main(int argc, char* argv[]){
    
    mEdge r1 = run();
    mEdge r2 = run();
    assert(r1.compareNumerically(r2));
    
}
*/



static testing::AssertionResult matrixEqual(const mEdge& e, const MatrixXcf& m){
    std::size_t dim;
    std_complex** mm = e.getMatrix(&dim);
    if(dim != m.rows()){
        return testing::AssertionFailure()<<"mEdge dim: "<<dim<<" , Matrix row: "<< m.rows();
    }

    if(dim != m.cols()){
        return testing::AssertionFailure()<<"mEdge dim: "<<dim<<" , Matrix col: "<< m.cols();
    }


    for(auto i = 0; i < dim; i++){
        for(auto j = 0; j < dim; j++){
            if(mm[i][j] != m(i,j)){
                return testing::AssertionFailure()<<"mEdge["<<i<<"]["<<j<<"]= "<<mm[i][j]<<", Matrix["<<i<<"]["<<j<<"] = "<<m(i,j);
            }
        }
    }

    return testing::AssertionSuccess();

}

static testing::AssertionResult vectorEqual(const vEdge& e, const VectorXcf& v){
    std::size_t dim;
    std_complex* vv = e.getVector(&dim);
    if(dim != v.rows()){
        return testing::AssertionFailure()<<"vEdge dim: "<<dim<<" , Vector dim: "<< v.size();
    }


    for(auto i = 0; i < dim; i++){
        if(vv[i] != v(i)){
            return testing::AssertionFailure()<<"vEdge["<<i<<"] = "<<vv[i]<<", Vector["<<i<<"] = "<<v(i);
        }
    }

    return testing::AssertionSuccess();

}

static MatrixXcf makeEigenGate(QubitCount q, const Matrix2cf& g, Qubit target){
    
    MatrixXcf k(1,1); 
    k(0,0) = std::complex<float>(1.0,0.0);

    MatrixXcf id = Matrix2cf::Identity();

    Qubit i = 0;
    for(; i < target; i++){
        k = Eigen::kroneckerProduct(id, k).eval();    
    }

    k = Eigen::kroneckerProduct(g, k).eval();

    for(i = i + 1; i < q; i++){
        k = Eigen::kroneckerProduct(id,k).eval(); 
    }

    return k;


}



TEST(QddTest, ArithmeticTest){
    {
        //Identity Matrix
        Qubit q = 9;
        std::size_t dim = (1<<(q+1));
        mEdge ident = makeIdent(q);
        MatrixXcf m = MatrixXcf::Identity(dim,dim);
        ASSERT_TRUE(matrixEqual(ident,m));
    }

    {
        //makeGate
        mEdge e = makeGate(2, Hmat, 0);
        MatrixXcf m = makeEigenGate(2,qdd_test::Hmat, 0);
        ASSERT_TRUE(matrixEqual(e, m));
        e = makeGate(8, Sdagmat, 2);
        m = makeEigenGate(8,qdd_test::Sdagmat, 2);
        ASSERT_TRUE(matrixEqual(e, m));
        
    }

    {
        //makeGate with controls
        Controls controls;
        controls.emplace(Control{2, Control::Type::pos});
        mEdge e = makeGate(3, Hmat, 0, controls);

        MatrixXcf one = Matrix2cf::Zero();
        one(1,1) = {1.0,0.0};
        MatrixXcf zero = Matrix2cf::Zero();
        zero(0,0) = {1.0,0.0};
        MatrixXcf ident = Matrix2cf::Identity();
        MatrixXcf h = makeEigenGate(1, qdd_test::Hmat, 0);

        MatrixXcf r1 = Eigen::kroneckerProduct(one, ident);
        r1 = Eigen::kroneckerProduct(r1, h).eval();
        MatrixXcf r2 = Eigen::kroneckerProduct(zero, ident).eval();
        r2 = Eigen::kroneckerProduct(r2, ident).eval();
        MatrixXcf result =  r1 + r2;
        ASSERT_TRUE(matrixEqual(e, result));
        

    }

    mEdge le, re, result_e;
    le = makeGate(8, Vdagmat, 2);
    re = makeGate(8, Hmat, 4);
    MatrixXcf lm, rm, result_m;
    lm = makeEigenGate(8, qdd_test::Vdagmat,2);
    rm = makeEigenGate(8, qdd_test::Hmat,4);

    Worker w(10);

    //mm_add
    {
        result_e = mm_add(&w, le, re);
        result_m = lm + rm;
        ASSERT_TRUE(matrixEqual(result_e, result_m));
    }

    //mm_multiply
    {
        result_e = mm_multiply(&w, le, re);
        result_m = lm * rm;
        ASSERT_TRUE(matrixEqual(result_e, result_m));
    }

    //mv_multiply
    {
        vEdge v = makeZeroState(8); 
        vEdge result_v = mv_multiply(&w, le, v);

        VectorXcf vec =  VectorXcf::Zero(1<<8);
        vec(0) = std::complex{1.0,0.0};
        VectorXcf result_vec = lm * vec;
        ASSERT_TRUE(vectorEqual(result_v, result_vec));
    }



}
TEST(QddTest, CircuitTest){
    vEdge v;
    VectorXcf vec;
    {
        //compute using graph.hpp
        vEdge input = makeZeroState(5);
        QuantumCircuit qc(5,8, 10, input);
        qc.emplace_back(Hmat, 1);
        qc.emplace_back(Xmat, 2);
        qc.emplace_back(Ymat, 3);
        qc.emplace_back(Vdagmat, 4);
        qc.emplace_back(Hmat, 3);
        qc.emplace_back(Sdagmat, 2);
        qc.buildCircuit();
        v = qc.wait().vectorResult();

        //compute using Eigen
        VectorXcf vinput =  VectorXcf::Zero(1<<5);
        vinput(0) = std::complex{1.0,0.0};
        MatrixXcf m1 = makeEigenGate(5,qdd_test::Hmat, 1);
        MatrixXcf m2 = makeEigenGate(5,qdd_test::Xmat, 2);
        MatrixXcf m3 = makeEigenGate(5,qdd_test::Ymat, 3);
        MatrixXcf m4 = makeEigenGate(5,qdd_test::Vdagmat, 4);
        MatrixXcf m5 = makeEigenGate(5,qdd_test::Hmat, 3);
        MatrixXcf m6 = makeEigenGate(5,qdd_test::Sdagmat, 2);
        vec = m6 * m5 * m4 * m3 * m2 * m1 * vinput;
    }
    ASSERT_TRUE(vectorEqual(v, vec));
}


TEST(QddTest, GroverTest){

    QubitCount n_qubits = 2;
    QubitCount n_oracles = 1;
    QubitCount total_qubit = n_qubits + n_oracles;
    
    std::string oracle = "00"; 


    QuantumCircuit qc(total_qubit,2, 10);
    qc.setInput(makeZeroState(total_qubit));


    //prepare the oracle qubit to be 1
    qc.emplace_back(Xmat, n_qubits);


    //apply H to all qubits
    for(auto i = 0; i < total_qubit; i++) qc.emplace_back(Hmat, i);
    
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
    for(auto i = 0 ; i < n_qubits; i++){
        qc.emplace_back(Hmat, i);
    }

    qc.buildCircuit();
    vEdge result = qc.wait().vectorResult();
    result.printVector();


}
