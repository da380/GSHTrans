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

  auto f = ScalarField(grid);

  auto g = ScalarFieldView(f, 2);

  auto h = ScalarFieldView(g, 0.5);

  std::cout << f(2, 2) << " " << g(2, 2) << " " << h(2, 2) << std::endl;
}
