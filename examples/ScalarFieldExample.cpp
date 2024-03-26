#include <GSHTrans/All>
#include <iostream>

using namespace GSHTrans;

int main() {
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Real;
  using MRange = All;
  using NRange = All;

  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 8;
  auto nMax = 2;

  auto grid = Grid(lMax, nMax);

  auto [u, dataU] = AllocateRealField(grid);
  auto [v, dataV] = AllocateComplexField(grid);

  u = 1;
  v = [](auto theta, auto phi) { return 2; };

  for (auto val : real(u + v)) std::cout << val << std::endl;
}