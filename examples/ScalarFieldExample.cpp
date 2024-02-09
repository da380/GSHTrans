#include <GSHTrans/All>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

int main() {
  using namespace GSHTrans;
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Complex;
  using MRange = All;
  using NRange = All;
  using Vector = std::vector<Scalar>;

  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 4;
  auto nMax = 2;

  auto grid = std::make_shared<Grid>(lMax, lMax);

  auto f = [](auto theta, auto phi) {
    auto ii = Complex(0, 1);
    auto fac = std::sqrt(static_cast<Real>(3) / static_cast<Real>(8)) *
               std::numbers::inv_sqrtpi_v<Real>;
    return fac * std::sin(theta) * std::exp(ii * phi);
  };

  auto data = std::vector<Complex>(grid->ComponentSize());
  auto y = std::ranges::views::all(data) | FormCanonicalComponentView(grid, 0);
  auto x = ComplexCanonicalComponent(grid, 0);
  y = 2 * y + 1;
  x = (y + y) / 2;

  std::cout << Integrate(x * conj(y)) << std::endl;
}
