#include "labs/DerivativeCoef/DerivativeCoef.h"
#include <fstream>
#include <iostream>
#include <numbers>

template <u32 N, u32 L, typename Callable>
double calcDerivative(double const h, double const x, Callable const &func) {
	std::array<double, N-1> points;
	for (u32 i = 0; i < N-1; ++i)
		points[i] = i+1;

	auto coefs = calcDerivativeCoef<double, N-1, L>(points);
	double res = func(x) * coefs.centralCoef;
	for (u32 i = 0; i < N-1; ++i)
		res += func(x + h * points[i]) * coefs.otherCoefs[i];
	return res / std::pow(h, L);
}

double Exp(double const &x)
{
	return std::exp(x);
}

template <u32 Nh, u32 N, u32 L, typename Callable>
std::array<double, Nh> compDerivatives(std::array<double, Nh> const &h, double const x, Callable const &func)
{
	std::array<double, Nh> computed;
	for (u32 i = 0; i < computed.size(); ++i)
		computed[i] = calcDerivative<N, L>(h[i], x, func);
	return computed;
}

int main()
{
	double const x = 1;
	u32 constexpr N = 3;
	u32 const L = 2;
	std::array<double, 16> h; h[0] = 1;
	for (u32 i = 1; i < h.size(); ++i)
		h[i] = h[i-1]/10;

	auto y= compDerivatives<h.size(), N, L>(h, x, Exp);
	for (u32 i = 0; i < y.size(); ++i)
		std::cout << std::log(h[i]) << ", ";
	std::cout << std::endl;
	for (u32 i = 0; i < h.size(); ++i)
		std::cout << std::log(std::abs(y[i] - std::numbers::e_v<double>)) << ", ";
}


