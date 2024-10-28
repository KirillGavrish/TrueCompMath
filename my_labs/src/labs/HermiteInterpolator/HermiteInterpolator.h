#ifndef INTERPOLANT_H
#define INTERPOLANT_H
#include <cstdint>
#include <array>

using u32 = std::uint32_t;

template<typename xType, typename yType, u32 N>
class HermiteInterpolator {
    std::array<xType, 2 * N> z;
    std::array<yType, 2 * N> coefs;
public:
    HermiteInterpolator(std::array<xType, N> const &, std::array<yType, N> const &, std::array<yType, N> const &) noexcept;
    yType interpolate(xType const &) const noexcept;
};

template<typename xType, typename yType, u32 N>
HermiteInterpolator<xType, yType, N>::HermiteInterpolator(
    std::array<xType, N> const &points,
    std::array<yType, N> const &values,
    std::array<yType, N> const &deriv
    ) noexcept
{
    for (int i = 0; i < N; ++i)
    {
        z[2 * i] = points[i];
        z[2 * i + 1] = points[i];

        coefs[2 * i] = values[i];
        coefs[2 * i + 1] = values[i];
    }
    for (u32 i = 1; i < 2 * N; ++i)
    {
        for (u32 j = 2 * N - 1; j >= i; --j)
        {
            if (z[j] == z[j - i])
                coefs[j] = deriv[j-i];
            else
                coefs[j] = (coefs[j] - coefs[j - 1]) / (z[j] - z[j - i]);
        }
    }
}

template<typename xType, typename yType, u32 N>
yType HermiteInterpolator<xType, yType, N>::interpolate(xType const &x) const noexcept
{
    yType interpolated, tmp;
    for(u32 i = 0; i < 2 * N; ++i)
    {
        tmp = coefs[i];
        for(u32 j = 0; j < i; j++)
            tmp *= (x - z[j]);
        interpolated += tmp;
    }
    return interpolated;
}

#endif //INTERPOLANT_H
