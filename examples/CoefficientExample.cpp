
#include <GSHTrans/Coefficient>
#include <GSHTrans/Core>
#include <GSHTrans/Field>
#include <iostream>

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

  auto u = ScalarFieldCoefficient<Grid, ComplexValued>(grid);

  u[2, 1] = Complex{1, 2};

  auto v = -u;

  std::cout << v << std::endl;
}