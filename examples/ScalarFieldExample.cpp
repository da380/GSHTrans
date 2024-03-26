#include <GSHTrans/All>
#include <iostream>

using namespace GSHTrans;

int main() {
  using Real = double;
  using Complex = std::complex<Real>;
  using MRange = All;
  using NRange = All;

  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 8;
  auto nMax = 2;

  auto grid = Grid(lMax, nMax);

  auto size = grid.ComponentSize();
  auto data = std::vector<Complex>(size);

  auto u = Field(data, grid);

  u = 1;

  for (auto val : real(u) + u) std::cout << val << std::endl;
}