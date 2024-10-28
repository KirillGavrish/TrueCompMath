#include "../src/labs/CubicSpline/CubicSpline.h"
#include "../third_party/googletest/googletest/include/gtest/gtest.h"

TEST(CubicSpline, TEST)
{
    std::vector<double> points = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0,7, 0.8, 0.9, 1};
    std::vector<double> values = {
        1,
        1.105170918075647,
        1.221402758160169,
        1.349858807576003,
        1.491824697641270,
        1.648721270700128,
        1.822118800390508,
        2.013752707470476,
        2.225540928492467,
        2.459603111156949,
        2.718281828459045
    };
    auto const interpolated = CubicSpline<double, double>(points, values, 0, 0);
    ASSERT_NEAR(interpolated.interpolate(0.51), 1.6652911949458, 1e-2);
}