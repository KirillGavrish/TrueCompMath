#include "Eigen/Dense"
#include<cstdint>
#include<array>

using u32 = std::uint32_t;

template<typename RealType, u32 N>
struct DerivativeCoef
{
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;

    DerivativeCoef(RealType const &, Eigen::Vector<RealType, N> const &);
};

template<typename RealType, u32 N>
DerivativeCoef<RealType, N>::DerivativeCoef(RealType const &otherСentralCoef, Eigen::Vector<RealType, N> const &otherOtherCoefs)
        : centralCoef(otherСentralCoef)
{
    for (u32 i = 0; i < otherCoefs.size(); ++i)
        otherCoefs[i] = otherOtherCoefs(i);
}

template<typename RealType, u32 N, u32 L>
DerivativeCoef<RealType, N> calcDerivativeCoef(std::array<RealType, N> const &points) noexcept
{
    Eigen::Matrix<RealType, N + 1, N + 1> A = Eigen::Matrix<RealType, N + 1, N + 1>::Zero();
    A.row(0) = Eigen::Vector<RealType, N + 1>::Ones();

    for (u32 i = 1; i < A.rows(); ++i)
    {
        for (u32 j = 1; j < A.cols(); ++j)
                 A(i, j) = A(i-1,j) * points[j-1] / i;
    }
    Eigen::Vector<RealType, N + 1> b = Eigen::Vector<RealType, N + 1>::Zero(); b(L)= 1;
    Eigen::Vector<RealType, N + 1> x = A.lu().solve(b);
    return {x(0), x.segment(1, N)};
}


