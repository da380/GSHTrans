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
  auto nMax = 0;

  auto grid = Grid(lMax, nMax);

  auto x =
      Testing::Component<Grid, Testing::Coefficient, Testing::RealValued>(grid);

  for (auto val : x) std::cout << val << std::endl;

  // auto x = Component<Grid, RealValued>(grid);

  // std::cout << x.UpperIndex() << std::endl;
}