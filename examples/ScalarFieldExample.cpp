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

  // auto grid = std::make_shared<Grid>(lMax, nMax);

  auto grid = std::make_shared<Grid>(lMax, lMax);

  auto data = std::vector<Real>(grid->ComponentSize());

  auto v = std::ranges::views::all(data) | CanonicalComponent(grid);

  v = InterpolateCanonicalComponent(
      grid, [](auto theta, auto phi) -> Real { return 1; });

  std::cout << Integrate(v + 7) << std::endl;
}
