#include "gtest/gtest.h"

#include "common.h"
#include "dd.h"

bool isNearlyEqual(std_complex lhs, std::complex<double> rhs){
    // Here, tolerance is larger than dd.h
    double TOL = 0.000001;
    if (lhs.r == rhs.real() && lhs.i == rhs.imag())
            return true;
    if (std::abs(lhs.r - rhs.real()) <= TOL &&
        std::abs(lhs.i - rhs.imag()) <= TOL)
        return true;

    return false;
}

TEST(QddTest, GateTest){
    {
        mEdge m = makeGate(1, Xmat, 0);
        size_t dim;
        std_complex**  mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[0][1] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[1][0] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[1][1] == (std::complex<double>(0.0, 0.0)));
    }
    {
        mEdge m = makeGate(1, Ymat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, -1.0)));
        ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, 1.0)));
        ASSERT_TRUE(mat[1][1] == (std::complex<double>(0.0, 0.0)));
    }
    {
        mEdge m = makeGate(1, Zmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[1][1] == (std::complex<double>(-1.0, 0.0)));
    }
    {
        mEdge m = makeGate(1, Hmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        double val = 1.0 / std::sqrt(2);
        ASSERT_TRUE(isNearlyEqual(mat[0][0], {val, 0.0}));
        ASSERT_TRUE(isNearlyEqual(mat[0][1], {val, 0.0}));
        ASSERT_TRUE(isNearlyEqual(mat[1][0], {val, 0.0}));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], {-val, 0.0}));
    }
    {
        mEdge m = makeGate(2, Xmat, 0, Controls{Control{(Qubit)1, Control::Type::pos}});
        size_t dim;
        std_complex**  mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[0][2] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[0][3] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[1][1] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[1][2] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[1][3] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[2][0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[2][1] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[2][2] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[2][3] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[3][0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[3][1] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[3][2] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[3][3] == (std::complex<double>(0.0, 0.0)));
    }
    {
        mEdge m = makeGate(1, Smat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == std::complex<double>(1.0, 0.0));
        ASSERT_TRUE(mat[0][1] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(mat[1][0] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], std::exp(std::complex<double>(0,PI/2))));
    }
    {
        mEdge m = makeGate(1, Sdagmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == std::complex<double>(1.0, 0.0));
        ASSERT_TRUE(mat[0][1] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(mat[1][0] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], std::exp(std::complex<double>(0,-PI/2))));
    }
    {
        mEdge m = makeGate(1, Tmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == std::complex<double>(1.0, 0.0));
        ASSERT_TRUE(mat[0][1] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(mat[1][0] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], std::exp(std::complex<double>(0,PI/4))));
    }
    {
        mEdge m = makeGate(1, Tdagmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == std::complex<double>(1.0, 0.0));
        ASSERT_TRUE(mat[0][1] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(mat[1][0] == std::complex<double>(0.0, 0.0));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], std::exp(std::complex<double>(0,-PI/4))));
    }
    {
        mEdge m = makeGate(1, SXmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(isNearlyEqual(mat[0][0], {0.5,0.5}));
        ASSERT_TRUE(isNearlyEqual(mat[0][1], {0.5,-0.5}));
        ASSERT_TRUE(isNearlyEqual(mat[1][0], {0.5,-0.5}));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], {0.5,0.5}));
    }
    {
        mEdge m = makeGate(1, SXdagmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(isNearlyEqual(mat[0][0], {0.5,-0.5}));
        ASSERT_TRUE(isNearlyEqual(mat[0][1], {0.5,0.5}));
        ASSERT_TRUE(isNearlyEqual(mat[1][0], {0.5,0.5}));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], {0.5,-0.5}));
    }
    // V dag
    {
        mEdge m = makeGate(1, Vmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(isNearlyEqual(mat[0][0], {std::cos(PI/4),0}));
        ASSERT_TRUE(isNearlyEqual(mat[0][1], {0,-std::sin(PI/4)}));
        ASSERT_TRUE(isNearlyEqual(mat[1][0], {0,-std::sin(PI/4)}));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], {std::cos(PI/4),0}));
    }
    {
        mEdge m = makeGate(1, Vdagmat, 0);
        size_t dim;
        std_complex** mat = m.getMatrix(&dim);
        ASSERT_TRUE(isNearlyEqual(mat[0][0], {std::cos(-PI/4),0}));
        ASSERT_TRUE(isNearlyEqual(mat[0][1], {0,-std::sin(-PI/4)}));
        ASSERT_TRUE(isNearlyEqual(mat[1][0], {0,-std::sin(-PI/4)}));
        ASSERT_TRUE(isNearlyEqual(mat[1][1], {std::cos(-PI/4),0}));
    }
    {
        for (int i = 0; i<16; i++){
            float angle = PI * 2 * i / 16;
            mEdge m = RX(1, 0, angle);
            size_t dim;
            std_complex** mat = m.getMatrix(&dim);
            ASSERT_TRUE(isNearlyEqual(mat[0][0], {std::cos(angle / 2), 0}));
            ASSERT_TRUE(isNearlyEqual(mat[0][1], {0,-std::sin(angle/2)}));
            ASSERT_TRUE(isNearlyEqual(mat[1][0], {0,-std::sin(angle/2)}));
            ASSERT_TRUE(isNearlyEqual(mat[1][1], {std::cos(angle/2),0}));
        }
    }
    {
        for (int i = 0; i<16; i++){
            float angle = PI * 2 * i / 16;
            mEdge m = RY(1, 0, angle);
            size_t dim;
            std_complex** mat = m.getMatrix(&dim);
            ASSERT_TRUE(isNearlyEqual(mat[0][0], {std::cos(angle / 2), 0}));
            ASSERT_TRUE(isNearlyEqual(mat[0][1], {-std::sin(angle/2),0}));
            ASSERT_TRUE(isNearlyEqual(mat[1][0], {std::sin(angle/2),0}));
            ASSERT_TRUE(isNearlyEqual(mat[1][1], {std::cos(angle/2),0}));
        }
    }
    {
        for (int i = 0; i<16; i++){
            float angle = PI * 2 * i / 16;
            mEdge m = RZ(1, 0, angle);
            size_t dim;
            std_complex** mat = m.getMatrix(&dim);
            ASSERT_TRUE(isNearlyEqual(mat[0][0], std::exp(-std::complex<double>(0,angle/2))));
            ASSERT_TRUE(isNearlyEqual(mat[0][1], {0,0}));
            ASSERT_TRUE(isNearlyEqual(mat[1][0], {0,0}));
            ASSERT_TRUE(isNearlyEqual(mat[1][1], std::exp(std::complex<double>(0,angle/2))));
        }
    }
}

TEST(QddTest, InitialStateTest){
    {
        // ZeroState (1)
        vEdge v = makeZeroState(1);
        size_t dim;
        std_complex *vec = v.getVector(&dim);
        ASSERT_TRUE(vec[0] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(vec[1] == (std::complex<double>(0.0, 0.0)));
    }
    {
        // OneState (1)
        vEdge v = makeOneState(1);
        size_t dim;
        std_complex *vec = v.getVector(&dim);
        ASSERT_TRUE(vec[0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(vec[1] == (std::complex<double>(1.0, 0.0)));
    }
    {
        // ZeroState (1)
        vEdge v = makeZeroState(2);
        size_t dim;
        std_complex *vec = v.getVector(&dim);
        ASSERT_TRUE(vec[0] == (std::complex<double>(1.0, 0.0)));
        for(int i=1; i<dim; i++)
            ASSERT_TRUE(vec[i] == (std::complex<double>(0.0, 0.0)));
    }
    {
        // OneState (1)
        vEdge v = makeOneState(2);
        size_t dim;
        std_complex *vec = v.getVector(&dim);
        for(int i=0; i<dim-1; i++)
            ASSERT_TRUE(vec[i] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(vec[3] == (std::complex<double>(1.0, 0.0)));
    }
}

TEST(QddTest, AddTest){
    {
        // Matrix + Matrix (2x2)
        mEdge m1 = makeGate(1, Xmat, 0);
        mEdge m2 = makeGate(1, Zmat ,0);
        mEdge m3 = mm_add(m1, m2);
        size_t dim;
        std_complex**  mat = m3.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[0][1] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[1][0] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[1][1] == (std::complex<double>(-1.0, 0.0)));
    }

    {
        // Vector + Vector (2x1)
        vEdge zero = makeZeroState(1);
        vEdge one = makeOneState(1);
        zero.w = {0.6,0};
        one.w = {0,-0.8};
        vEdge v = vv_add(zero, one);
        size_t dim;
        std_complex* vec = v.getVector(&dim);
        ASSERT_TRUE(vec[0].isApproximatelyEqual({0.6,0.0}));
        ASSERT_TRUE(vec[1].isApproximatelyEqual({0.0,-0.8}));
    }
}

TEST(QddTest, MulTest){
    {
        // Matrix * Matrix (2x2)
        mEdge m1 = makeGate(1, Xmat, 0);
        mEdge m2 = makeGate(1, Zmat ,0);
        mEdge m3 = mm_multiply(m1, m2);
        size_t dim;
        std_complex**  mat = m3.getMatrix(&dim);
        ASSERT_TRUE(mat[0][0] == (std::complex<double>(0.0, 0.0)));
        ASSERT_TRUE(mat[0][1] == (std::complex<double>(-1.0, 0.0)));
        ASSERT_TRUE(mat[1][0] == (std::complex<double>(1.0, 0.0)));
        ASSERT_TRUE(mat[1][1] == (std::complex<double>(0.0, 0.0)));
    }
    {
        // Matrix * Vector (2)
        mEdge m1 = makeGate(1, SXmat, 0);
        vEdge zero = makeZeroState(1);
        vEdge v = mv_multiply(m1, zero);
        size_t dim;
        std_complex* vec = v.getVector(&dim);
        ASSERT_TRUE(vec[0] == (std::complex<double>(0.5, 0.5)));
        ASSERT_TRUE(vec[1] == (std::complex<double>(0.5, -0.5)));
    }
}

TEST(QddTest, MeasureTest){
    std::mt19937_64 mt(0);
    {
        vEdge state = makeZeroState(2);
        state = mv_multiply(makeGate(2,Hmat,0),state);
        state = mv_multiply(CX(2, 1, 0), state);
        std::map<std::string, int> resultmap;
        for (int i = 0; i < 100; i++){
            std::string result = measureAll(state, false, mt);
            if(resultmap.find(result) != resultmap.end()){
                resultmap[result] += 1;
            }else{
                resultmap[result] = 1;
            }
        }
        ASSERT_TRUE(resultmap["11"] > 40 && resultmap["11"] < 60);
        ASSERT_TRUE(resultmap["00"] > 40 && resultmap["00"] < 60);
    }

    {
        vEdge state = makeZeroState(2);
        state = mv_multiply(makeGate(2,Xmat,0),state);
        std::map<std::string, int> resultmap;
        for (int i = 0; i < 100; i++){
            std::string result = measureAll(state, false, mt);
            if(resultmap.find(result) != resultmap.end()){
                resultmap[result] += 1;
            }else{
                resultmap[result] = 1;
            }
        }
        ASSERT_TRUE(resultmap["01"] == 100);
    }

    {
        bool isZeroChecked = false, isOneChecked = false;
        while(isZeroChecked == false || isOneChecked == false){
            vEdge state = makeZeroState(3);
            state = mv_multiply(makeGate(3,Hmat,0),state);
            state = mv_multiply(CX(3, 1, 0), state);
            state = mv_multiply(CX(3, 2, 0), state);
            const auto res = measureOneCollapsing(state, 0, mt);
            if(res == '0'){
                for (int i = 0; i < 10; i++){
                    std::string allres = measureAll(state, false, mt);
                    ASSERT_TRUE(allres == "000");
                }
                isZeroChecked = true;
            }else if (res == '1'){
                for (int i = 0; i < 10; i++){
                    std::string allres = measureAll(state, false, mt);
                    ASSERT_TRUE(allres == "111");
                }
                isOneChecked = true;
            }else{
                FAIL();
            }
        }
    }

    {
        bool isZeroChecked = false, isOneChecked = false;
        std::map<char, int> resultmap{{'0',0},{'1',0}};
        for(int i=0;i<100;i++){
            vEdge state = makeZeroState(3);
            state = mv_multiply(RX(3,0,PI/3),state);
            state = mv_multiply(CX(3, 1, 0), state);
            state = mv_multiply(CX(3, 2, 0), state);
            const auto res = measureOneCollapsing(state, 0, mt);
            ASSERT_TRUE(resultmap.find(res) != resultmap.end());
            resultmap[res] += 1;
        }
        ASSERT_TRUE(resultmap['0'] < 80 && resultmap['0'] > 70);
    }
    
}

TEST(QddTest, DotTest){
    {
        vEdge state = makeZeroState(2);
        state = mv_multiply(makeGate(2,Hmat,0),state);
        state = mv_multiply(CX(2, 1, 0), state);
        std::string dot = genDot(state);
    }
}

TEST(QddTest, GCTest){
    {
        vEdge state = makeZeroState(2);
        state = mv_multiply(makeGate(2,Hmat,0),state);
        state = mv_multiply(CX(2, 1, 0), state);
        state = gc(state, true);
        state = mv_multiply(CX(2, 0, 1), state);
    }
    {
        mEdge h = makeGate(3,Hmat,0);
        mEdge cx = CX(3, 1, 0);
        mEdge mul = mm_multiply(h,cx);
        mul = gc_mat(mul, true);
        mEdge h2 = makeGate(3,Hmat,2);
        mEdge cx2 = CX(3, 2, 0);
        mEdge mul2 = mm_multiply(h2,cx2);
        mEdge mul3 = mm_multiply(mul2,mul);
    }
}