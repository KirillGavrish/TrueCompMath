#include "labs/CubicSpline/CubicSpline.h"
#include <cmath>
#include <iostream>

int main()
{
    double constexpr l = 10;

    std::vector<u32> N(6); N[0] = 5;
    for (u32 i = 1; i < N.size(); ++i)
        N[i] = N[i-1]*2;

    std::vector<std::vector<double>> points(N.size());
    for (u32 i = 0; i < N.size(); ++i)
    {
        points[i] = std::vector<double>(N[i]);
        for (u32 j = 0; j < N[i]; ++j)
            points[i][j] = j * l / (N[i]-1);
    }
    std::vector<double> errors(N.size());
    for (u32 i = 0; i < N.size(); ++i)
    {
        auto values = std::vector<double>(N[i]);
        for (u32 j = 0; j < N[i]; ++j)
            values[j] = std::exp(points[i][j]);
        auto interpolated = CubicSpline<double, double>(points[i], values, 1, std::exp(10));
        std::vector<double> web(1000);
        for (u32 j = 0; j < web.size(); ++j)
            web[j] = j * l / (999);
        double maxError = 0;
        for (auto p : web)
        {
            auto err = std::abs(interpolated.interpolate(p) - std::exp(p));
            if (err > maxError) maxError = err;
        }
        errors[i] = maxError;
    }
    for (auto n : N)
        std::cout << std::log(n) << ", ";
    std::cout << std::endl;
    for (auto err : errors)
        std::cout << std::log(err) << ", ";
}