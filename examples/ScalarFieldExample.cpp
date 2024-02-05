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

  auto grid = std::make_shared<Grid>(lMax, nMax);

  for (auto [theta, phi] : grid->Points())
    std::cout << theta << " " << phi << std::endl;

  auto data = std::vector<Complex>(grid->ComponentSize(), 1);

  auto v = std::ranges::views::all(data) | Views::CanonicalComponent(grid);

  auto w = 2 * v + 1;

  // for (auto val : w) std::cout << val << std::endl;
}
