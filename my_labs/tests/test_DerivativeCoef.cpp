#include "../third_party/googletest/googletest/include/gtest/gtest.h"
#include "../src/labs/DerivativeCoef/DerivativeCoef.h"

TEST(calcDerivativeCoef, TEST_EXAMPLE)
{
    std::array<double, 2> constexpr a = {-1, 1};
    DerivativeCoef<double, 2> const coefs = calcDerivativeCoef<double, 2, 2>(a);
    EXPECT_EQ(coefs.centralCoef, -2);
    std::array<double, 2> constexpr expectedOthercoefs = {1, 1};
    EXPECT_EQ(coefs.otherCoefs, expectedOthercoefs);
}

