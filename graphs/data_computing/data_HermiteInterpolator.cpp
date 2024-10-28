#include "labs/HermiteInterpolator/HermiteInterpolator.h"
#include <fstream>
#include <iostream>
#include <cmath>

int main()
{
    u32 constexpr N = 5;
    std::array<double, 6> l; l[0] = 2;
    for (u32 i = 1; i < l.size(); ++i)
        l[i] = l[i-1]/2;

    std::array<std::array<double, N>, l.size()> points;
    for (u32 i = 0; i < l.size(); ++i)
    {
        for (u32 j = 0; j < N; ++j)
            points[i][j] = j * l[i] / (N-1);
    }
    std::array<double, N> values;
    std::array<double, l.size()> errors;
    for (u32 i = 0; i < l.size(); ++i)
    {
        for (u32 j = 0; j < N; ++j)
            values[j] = std::exp(points[i][j]);
        auto interpolated = HermiteInterpolator<double, double, N>(points[i], values, values);
        std::array<double, 1000> web;
        for (u32 j = 0; j < web.size(); ++j)
            web[j] = j * l[i] / (web.size() - 1);
        double maxError = 0;
        for (auto p : web)
        {
            auto err = std::abs(interpolated.interpolate(p) - std::exp(p));
            if (err > maxError) maxError = err;
        }
        errors[i] = maxError;
    }
    for (u32 i = 0; i < l.size(); ++i)
        std::cout << l[i] << ", ";
    std::cout << std::endl;
    for (u32 i = 0; i < errors.size(); ++i)
        std::cout << errors[i] << ", ";
}
