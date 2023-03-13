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
}
