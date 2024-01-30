#include <GSHTrans/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
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
  using Scalar = Complex;
  using MRange = All;
  using NRange = All;
  using Vector = std::vector<Scalar>;

  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = 4;
  auto nMax = 2;

  auto grid = Grid(lMax, nMax);

  auto f = CanonicalComponent<Grid, RealValued>(grid, 1);

  auto g = RealCanonicalComponent(grid, [](auto theta, auto phi) { return 1; });

  auto h = ComplexCanonicalComponent(grid, 4);

  f = 2 * g;

  auto n = 0;

  auto flm = RealCanonicalCoefficient(lMax, n, grid);

  for (auto& val : flm) val = 1;
}
