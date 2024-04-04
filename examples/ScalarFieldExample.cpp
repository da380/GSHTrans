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

  auto u = ScalarField<Grid, RealValued>(grid, 1);

  auto v = u * Complex(1);

  auto w = conj(v);

  w.Print();
}