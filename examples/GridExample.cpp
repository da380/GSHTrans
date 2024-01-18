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
  using Type = C2C;
  using MRange = Type::IndexRange;
  using Grid = GaussLegendreGrid<Real, Type>;
  using Scalar = Grid::scalar_type;

  auto lMax = 64;
  auto nMax = 2;

  /*

  auto n = -2;
  auto grid = Grid(lMax, nMax);
  auto f = grid.FieldVector<Scalar>();
  auto flm = grid.RandomCoefficientVector(n);
  auto glm = grid.CoefficientVector(n);

  grid.InverseTransformation(n, flm.begin(), f.begin());
  grid.ForwardTransformation(n, f.begin(), glm.begin());

  auto flmView = GSHView<Complex, MRange>(lMax, lMax, n, flm.begin());
  auto glmView = GSHView<Complex, MRange>(lMax, lMax, n, glm.begin());

  for (auto [l, m] : glmView.Indices()) {
    auto err = std::abs(flmView(l)(m) - glmView(l)(m));
    if (err > 1e-10)
      std::cout << l << " " << m << " " << flmView(l)(m) << glmView(l)(m)
                << std::endl;
  }

  */
}
