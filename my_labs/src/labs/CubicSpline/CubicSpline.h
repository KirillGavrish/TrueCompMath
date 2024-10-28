#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H
#include <vector>
#include <type_traits>
#include <cstdint>

using u32 = std::uint32_t;

template<typename T>
class ThreeDiagonalMatrix
{
private:
    std::vector<T> a_;
    std::vector<T> b_;
    std::vector<T> c_;
public:
    ThreeDiagonalMatrix(std::vector<T> const &a, std::vector<T> const &b, std::vector<T> const &c)
        : a_(a), b_(b), c_(c)
    {}

    std::vector<T> const &a() const;
    std::vector<T> const &b() const;
    std::vector<T> const &c() const;

    T const &a(std::size_t const &) const;
    T const &b(std::size_t const &) const;
    T const &c(std::size_t const &) const;
};

template<typename T>
std::vector<T> const &ThreeDiagonalMatrix<T>::a() const {return a_;}
template<typename T>
std::vector<T> const &ThreeDiagonalMatrix<T>::b() const {return b_;}
template<typename T>
std::vector<T> const &ThreeDiagonalMatrix<T>::c() const {return c_;}

template<typename T>
T const &ThreeDiagonalMatrix<T>::a(std::size_t const &i) const {return a_[i];}
template<typename T>
T const &ThreeDiagonalMatrix<T>::b(std::size_t const &i) const {return b_[i];}
template<typename T>
T const &ThreeDiagonalMatrix<T>::c(std::size_t const &i) const {return c_[i];}

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solveThreeDiagonalMatrix(ThreeDiagonalMatrix<mType> const &A, std::vector<cType> const &b)
{
    std::vector<mType> p(b.size());
    std::vector<DivisType<cType, mType>> q(b.size());
    std::vector<DivisType<cType, mType>> x(b.size());

    p[0] = -A.c(0) / A.b(0);
    q[0] = b[0] / A.b(0);
    for (u32 i = 1; i < b.size() - 1; ++i)
    {
        p[i] = -(A.c(i) / (A.a(i - 1) * p[i - 1] + A.b(i)));
        q[i] = (b[i] - A.a(i - 1) * q[i - 1]) / (A.a(i - 1) * p[i - 1] + A.b(i));
    }
    q[b.size() - 1] = (b[b.size() - 1] - A.a(b.size() - 2) * q[b.size() - 2]) /
                      (A.a(b.size() - 2) * p[b.size() - 2] + A.b(b.size() - 1));
    x[b.size() - 1] = q[b.size() - 1];
    for (std::size_t i = 1; i < b.size(); ++i)
        x[b.size() - 1 - i] = p[b.size() - 1 - i] * x[b.size() - i] + q[b.size() - 1 - i];
    return x;
}

template<typename xType, typename yType>
class CubicSpline
{
    std::vector<xType> points_;
    u32 size_;

    using DeltaXType = DiffType<xType>;
    using DerivType = DivisType<DiffType<yType>, DeltaXType>;
    using Deriv2Type = DivisType<DiffType<DerivType>, DeltaXType>;

    struct Coeffs
    {
        yType a;
        yType b;
        yType c;
        yType d;
    };
    std::vector<Coeffs> coeffs_;
public:
    CubicSpline(std::vector<xType> const &, std::vector<yType> const &, Deriv2Type const &, Deriv2Type const &);
    yType interpolate(xType const &) const noexcept;
};

template<typename xType, typename yType>
CubicSpline<xType, yType>::CubicSpline(std::vector<xType> const &points,
    std::vector<yType> const &values,
    Deriv2Type const &first,
    Deriv2Type const &second
    ) : points_(points)
{
    size_ = points.size();
    coeffs_.resize(size_);

    coeffs_[0].c = first;
    coeffs_[size_ - 1].c = second;

    std::vector<yType> u(size_ - 2);
    std::vector<xType> a(size_ - 3), b(size_ - 2, 2), c(size_ - 3);

    std::vector<DeltaXType> h(size_);
    std::vector<DiffType<yType>> sepDifferences(size_);

    for(int i = 1; i < size_; i++){
        h[i] = points_[i] - points_[i - 1];
        sepDifferences[i] = (values[i] - values[i - 1]) / h[i];
    }

    u[0] = 6 * (sepDifferences[2] - sepDifferences[1]) / (points_[2] - points_[0]);
    for(int i = 0; i < size_ - 3; i++){
        a[i] = h[i + 2] / (h[i + 2] + h[i + 1]);
        c[i] = h[i + 1] / (h[i + 2] + h[i + 1]);
        u[i + 1] = 6 * (sepDifferences[i + 3] - sepDifferences[i + 2]) / (points_[i + 3] - points_[i + 1]);
    }

    ThreeDiagonalMatrix<xType> threeDiagM(a, b, c);
    std::vector<Deriv2Type> solution = solveThreeDiagonalMatrix(threeDiagM, u);

    for(int i = 1; i < size_ - 1; i++){
        coeffs_[i].c = solution[i - 1];
        coeffs_[i].d = (coeffs_[i].c - coeffs_[i - 1].c) / h[i];
        coeffs_[i].a = values[i];
        coeffs_[i].b = h[i] / 3 * (coeffs_[i].c + coeffs_[i - 1].c / 2) + sepDifferences[i];
    }

    coeffs_[size_ - 1].a = values[size_ - 1];
    coeffs_[size_ - 1].b = h[size_ - 1] / 3 * (coeffs_[size_ - 1].c + coeffs_[size_ - 2].c / 2) + sepDifferences[size_ - 1];
    coeffs_[size_ - 1].d = (coeffs_[size_ - 1].c - coeffs_[size_ - 2].c) / h[size_ - 1];
}

template<typename xType, typename yType>
yType CubicSpline<xType, yType>::interpolate(xType const &x) const noexcept
{
    u32 n;
    for (u32 i = 1; i < size_; ++i)
    {
        if ((points_[i - 1] <= x) && (x <= points_[i]))
        {
            n = i;
            break;
        }
    }
    DiffType<xType> dx = (x - points_[n]);
    return coeffs_[n].a + coeffs_[n].b * dx + coeffs_[n].c / 2 * dx * dx + coeffs_[n].d / 6 * dx * dx * dx;
}

#endif //CUBICSPLINE_H
