#include "labs/integrate/integrate.h"
#include <cmath>
#include <iostream>
#include <vector>
double func(double const &x) {return std::sin(x);}

int main()
{
    double constexpr l = 10;
    double constexpr realI = 1.83907152907645;
    u32 constexpr N = 5;
    u32 constexpr xNum = 17;
    std::vector<double> errors(xNum);
    u32 constexpr RUNGE_FLAG = 1u;
    if constexpr (RUNGE_FLAG == 1u)
    {
        double eps_i = realI/10;
        std::vector<double> epsilon(xNum);
        for (u32 i = 0; i < xNum; ++i)
        {
            auto compI = integrate<N>(func, 0, l, eps_i, 1e10);
            errors[i] = std::abs(compI - realI);
            epsilon[i] = eps_i;
            eps_i /= 10;
        }
        for (auto eps : epsilon)
            std::cout << std::log(eps) << ", ";
        std::cout << std::endl;
        for (auto err : errors)
            std::cout << std::log(err) << ", ";
        std::cout << errors[1];
    }
    else
    {
        double hi = l;
        std::vector<double> h(xNum);
        for (u32 i = 0; i < xNum; ++i)
        {
            auto compI = integrate<N>(func, 0, l, hi);
            errors[i] = std::abs(compI - realI);
            h[i] = hi;
            hi /= 2;
        }
        for (auto step : h)
            std::cout << std::log(step) << ", ";
        std::cout << std::endl;
        for (auto err : errors)
            std::cout << std::log(err) << ", ";
    }
}
