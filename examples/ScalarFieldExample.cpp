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

  auto size = grid.FieldSize();
  auto data0 = std::vector<Real>(size);
  auto data1 = std::vector<Real>(size);
  auto data2 = std::vector<Real>(size);

  auto u = RealVectorField(data0, data1, data2, grid);

  u(0) = 2;

  auto w = u + u;

  for (auto val : w(0)) std::cout << val << std::endl;
}