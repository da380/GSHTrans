#include <GSHTrans/All>
#include <array>
#include <iostream>
#include <ranges>
#include <tuple>

using namespace GSHTrans;
int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Complex;
  using MRange = All;
  using NRange = All;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 8;
  auto nMax = 2;

  auto grid = Grid(lMax, nMax);

  auto u = VectorField<Grid, RealValued>(grid);

  auto v = u[-1];

  // std::cout << u[-1, 0, 0] << std::endl;

  std::cout << v[1, 2] << std::endl;
}