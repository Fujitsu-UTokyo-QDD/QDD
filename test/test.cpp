#include "gtest/gtest.h"

#include "common.h"
#include "dd.h"

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
        double val = 0.5 / std::sqrt(2);
        ASSERT_TRUE(mat[0][0] == (std::complex<double>(val, 0.0)));
        ASSERT_TRUE(mat[0][1] == (std::complex<double>(val, 0.0)));
        ASSERT_TRUE(mat[1][0] == (std::complex<double>(val, 0.0)));
        ASSERT_TRUE(mat[1][1] == (std::complex<double>(-val, 0.0)));
    }
   {
       mEdge m = makeGate(1, Imat, 0);
       size_t dim;
       std_complex** mat = m.getMatrix(&dim);
       ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
       ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][1] == (std::complex<double>(1.0, 0.0)));
   }
   {
       mEdge m = makeGate(1, Smat, 0);
       size_t dim;
       std_complex** mat = m.getMatrix(&dim);
       double s = std::sqrt(0.5);
       ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
       ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, -s)));
       ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, s)));
       ASSERT_TRUE(mat[1][1] == (std::complex<double>(1.0, 0.0)));
   }
   {
       mEdge m = makeGate(1, Sdagmat, 0);
       size_t dim;
       std_complex** mat = m.getMatrix(&dim);
       ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
       ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][1] == (std::complex<double>(0.0, -1.0)));
   }
   {
       mEdge m = makeGate(1, Tmat, 0);
       size_t dim;
       std_complex** mat = m.getMatrix(&dim);
       double val = 0.5 * (1.0 + std::sqrt(2));
       ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
       ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][1] == (std::complex<double>(val, val)));
   }
   {
       mEdge m = makeGate(1, Tdagmat, 0);
       size_t dim;
       std_complex** mat = m.getMatrix(&dim);
       double val = 1.0 / std::sqrt(2);
       ASSERT_TRUE(mat[0][0] == (std::complex<double>(1.0, 0.0)));
       ASSERT_TRUE(mat[0][1] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][0] == (std::complex<double>(0.0, 0.0)));
       ASSERT_TRUE(mat[1][1] == (std::complex<double>(val, -val)));
   }
   {
       mEdge m = makeGate(1, SXmat, 0);
       size_t dim;
       std_complex** mat = m.getMatrix(&dim);
       double val = 1.0 / std::sqrt(2);
       ASSERT_TRUE(mat[0][0] == (std::complex<double>(val, 0.0)));
       ASSERT_TRUE(mat[0][1] == (std::complex<double>(val, 0.0)));
       ASSERT_TRUE(mat[1][0] == (std::complex<double>(val, 0.0)));
       ASSERT_TRUE(mat[1][1] == (std::complex<double>(-val, 0.0)));
   }
   {
       m = makeGate(1, SXdagmat, 0);
       mat = m.getMatrix(&dim);
       val = 1.0 / std::sqrt(2);
       ASSERT_TRUE(mat[0][0] == (std::complex<double>(val, 0.0)));
       ASSERT_TRUE(mat[0][1] == (std::complex<double>(-val, 0.0)));
       ASSERT_TRUE(mat[1][0] == (std::complex<double>(val, 0.0)));
       ASSERT_TRUE(mat[1][1] == (std::complex<double>(val, 0.0)));
   }
   {
       mEdge m = makeGate(1, Vdagmat, 0);
       size_t dim;
       std_complex** mat = m.getMatrix(&dim);
       ASSERT_TRUE(mat[0][0] == cf_SQRT2_2);
       ASSERT_TRUE(mat[0][1] == cf_iSQRT2_2);
       ASSERT_TRUE(mat[1][0] == cf_iSQRT2_2);
       ASSERT_TRUE(mat[1][1] == cf_SQRT2_2);
   }
   {
        // todo Vmat
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
        state.printVector();
        std::map<std::string, int> resultmap;
        for (int i = 0; i < 100; i++){
            std::string result = measureAll(state, false, mt);
            if(resultmap.contains(result)){
                resultmap[result] += 1;
            }else{
                resultmap[result] = 1;
            }
        }
        for(auto itr: resultmap)
            std::cout << itr.first << ":" << itr.second << ", ";
        std::cout << std::endl;
        ASSERT_TRUE(resultmap["11"] > 40 && resultmap["11"] < 60);
        ASSERT_TRUE(resultmap["00"] > 40 && resultmap["00"] < 60);
    }

    {
        vEdge state = makeZeroState(2);
        state = mv_multiply(makeGate(2,Xmat,0),state);
        state.printVector();
        std::map<std::string, int> resultmap;
        for (int i = 0; i < 100; i++){
            std::string result = measureAll(state, false, mt);
            if(resultmap.contains(result)){
                resultmap[result] += 1;
            }else{
                resultmap[result] = 1;
            }
        }
        for(auto itr: resultmap)
            std::cout << itr.first << ":" << itr.second << ", ";
        std::cout << std::endl;
        ASSERT_TRUE(resultmap["01"] == 100);
    }

    {
        std::cout << "measureOneCollapsing test" << std::endl;
        bool isZeroChecked = false, isOneChecked = false;
        while(isZeroChecked == false || isOneChecked == false){
            vEdge state = makeZeroState(3);
            state = mv_multiply(makeGate(3,Hmat,0),state);
            state = mv_multiply(CX(3, 1, 0), state);
            state = mv_multiply(CX(3, 2, 0), state);
            const auto res = measureOneCollapsing(state, 0, true, mt);
            if(res == '0'){
                state.printVector();
                for (int i = 0; i < 10; i++){
                    std::string allres = measureAll(state, false, mt);
                    ASSERT_TRUE(allres == "000");
                }
                isZeroChecked = true;
            }else if (res == '1'){
                state.printVector();
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
        std::cout << "measureOneCollapsing test" << std::endl;
        bool isZeroChecked = false, isOneChecked = false;
        std::map<char, int> resultmap{{'0',0},{'1',0}};
        for(int i=0;i<100;i++){
            vEdge state = makeZeroState(3);
            state = mv_multiply(RX(3,0,std::numbers::pi/3),state);
            state = mv_multiply(CX(3, 1, 0), state);
            state = mv_multiply(CX(3, 2, 0), state);
            const auto res = measureOneCollapsing(state, 0, true, mt);
            ASSERT_TRUE(resultmap.contains(res));
            resultmap[res] += 1;
        }
        for(auto itr: resultmap)
            std::cout << itr.first << ":" << itr.second << ", ";
        std::cout << std::endl;
        ASSERT_TRUE(resultmap['0'] < 80 && resultmap['0'] > 70);
    }
    
}

TEST(QddTest, DotTest){
    {
        vEdge state = makeZeroState(2);
        state = mv_multiply(makeGate(2,Hmat,0),state);
        state = mv_multiply(CX(2, 1, 0), state);
        std::cout << genDot(state) << std::endl;
    }
}
