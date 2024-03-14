#include <GSHTrans/All>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
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

  auto grid = Grid(lMax, lMax);

  auto f = [](auto theta, auto phi) {
    constexpr auto ii = Complex(0, 1);
    auto fac = std::sqrt(static_cast<Real>(3) / static_cast<Real>(8)) *
               std::numbers::inv_sqrtpi_v<Real>;
    return fac * std::sin(theta) * std::exp(ii * phi);
  };

  auto data = std::vector<Complex>(grid.ComponentSize());
  auto y = grid.InterpolateFunction(f) | FormCanonicalComponentView(grid, 0);
  auto x = CanonicalComponent<Grid, ComplexValued>(grid, 0, f);

  x(2, 1) = 2;
  auto z = 2 * x;

  std::cout << Integrate(z * conj(y)) << std::endl;

  std::cout << z(1, 2) << std::endl;
}
