#include "../third_party/googletest/googletest/include/gtest/gtest.h"
#include "../src/labs/integrate/integrate.h"
#include <cmath>

double f1(double const &x) {return std::sin(x);}
double f2(double const &x) {return x*x;}
double f3(double const &x) {return x*std::sin(x);}

TEST(integrate, TEST_sin)
{
    auto const Int = integrate<2>(f1, 0, 10, 0.0000001);
    ASSERT_NEAR(Int, 1.83907152907645, 1e-7);
}

TEST(integrate, TEST_xx)
{
    auto const Int = integrate<4>(f2, 0, 10, 0.0000001);
    ASSERT_NEAR(Int, 333.33333333, 1e-5);
}

TEST(integrate, TEST_xsin)
{
    auto const Int = integrate<6>(f3, 0, 10, 0.0000001);
    ASSERT_NEAR(Int, 7.84669417987515, 1e-5);
}

TEST(integrateRichardsonExtrapolation, TEST_xx)
{
    auto const Int = integrateRichardsonExtrapolation<6>(f3, 0, 10, 0.0000001);
    ASSERT_NEAR(Int, 7.84669417987515, 1e-7);
}

TEST(integrate, TEST_RUNGE_xx)
{
    auto const Int = integrate<3>(f3, 0, 10, 1e-15, 1000);
    ASSERT_NEAR(Int, 7.84669417987515, 1e-14);
}
