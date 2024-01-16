#include <GSHTrans/All>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

int main() {
  using namespace GSHTrans;
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, R2C, Ortho>;
  using Scalar = Grid::scalar_type;

  auto lMax = 8;
  auto nMax = 2;

  auto grid = std::make_shared<Grid>(lMax, nMax);

  ScalarField x(grid);

  /*
  ScalarField x(grid, [](auto theta, auto phi) -> Scalar { return 1; });
  ScalarField y(grid, [](auto theta, auto phi) -> Scalar { return 1; });
  std::cout << x + y;
  */
}
