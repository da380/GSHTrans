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

  auto grid = std::make_shared<Grid>(lMax, nMax);

  auto size = grid->ComponentSize();
  auto data = std::vector<Complex>(size, 1);

  //  auto v = std::ranges::views::all(data) | Views::CanonicalComponent(grid);
  auto v = CanonicalComponentView(data, grid);
  auto w = 2 * v * (v + v / 3);

  auto z = MakeRealCanonicalComponent(grid);

  for (auto& val : z) val = 1;

  for (auto val : z) std::cout << val << std::endl;
}
