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

using namespace GSHTrans;

int main() {
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Real;
  using MRange = All;
  using NRange = All;
  using Vector = std::vector<Scalar>;

  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 8;
  auto nMax = 2;

  auto grid = Grid(lMax, nMax);

  auto size = grid.ComponentSize();
  auto data = std::vector<Scalar>(size);

  auto u = ComponentField(grid, data, 1);

  u = 1;

  auto v = u * u + 1;

  std::cout << v.Rank() << std::endl;
  for (auto val : v) std::cout << val << std::endl;
}