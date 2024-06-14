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

  auto data = std::vector<Real>(grid.FieldSize());

  auto v = RealVectorField(grid, [](auto theta, auto phi) {
    return CanonicalVector<Real>{1, 2, 3};
  });

  auto f =
      RealAbstractScalarField(grid, [](auto theta, auto phi) { return 2; });

  auto i = Complex{0, 1};
  std::cout << f * v + v / 2 << std::endl;
}
