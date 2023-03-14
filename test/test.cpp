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
        // Please write a test for Y gate...
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
