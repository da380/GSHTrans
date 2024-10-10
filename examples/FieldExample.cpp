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
  auto nMax = 2;

  auto grid = Grid(lMax, nMax);

  auto u = RealCanonicalComponentField(1, grid,
                                       [](auto theta, auto phi) { return 1; });
  auto v = RealCanonicalComponentField(1, grid,
                                       [](auto theta, auto phi) { return 5; });

  std::cout << u + 2 * v;

  //  auto ulm = RealScalarFieldExpansion(grid, [](auto l, auto m) { return 1;
  //  });

  // for (auto [l, m] : ulm.Indices()) std::cout << l << " " << m << std::endl;
}
