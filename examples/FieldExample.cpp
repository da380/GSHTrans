#include <GSHTrans/Core>
#include <GSHTrans/Field>
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

  auto lMax = 4;
  auto nMax = 2;

  auto grid = Grid(lMax, nMax);

  auto A = RealMatrixField(grid);

  auto B = RealMatrixField(grid);

  auto C = Tr(B);

  C.Print();
}