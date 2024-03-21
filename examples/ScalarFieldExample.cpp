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
  using Scalar = Real;
  using MRange = All;
  using NRange = All;
  using Vector = std::vector<Scalar>;

  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 8;
  auto nMax = 1;

  auto grid = Grid(lMax, nMax);

  auto x = CanonicalComponent<Grid, ComplexValued>(grid, 0, 1);

  x(2, 1) = 2;
  auto z = x * 2.0;

  for (auto val : z) std::cout << val << std::endl;
}