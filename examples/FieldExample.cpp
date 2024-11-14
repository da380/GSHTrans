#include <GSHTrans/All>
#include <array>
#include <format>
#include <iostream>
#include <print>
#include <ranges>
#include <tuple>

using namespace GSHTrans;
int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using MRange = All;
  using NRange = All;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 4;
  auto nMax = 0;

  auto grid = Grid(lMax, nMax);

  auto u = RealCanonicalComponentField<0, Grid>(
      grid, [](auto theta, auto phi) { return 2; });
  auto v = RealCanonicalComponentField<0, Grid>(
      grid, [](auto theta, auto phi) { return 5; });

  auto ulm = CanonicalComponentExpansion<0, Grid, ComplexValued>(grid);

  std::cout << ulm << std::endl;
}
