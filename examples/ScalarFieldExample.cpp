#include <GSHTrans/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
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
  using Grid = GaussLegendreGrid<Real, R2C>;
  using Scalar = Grid::scalar_type;

  auto lMax = 4;
  auto nMax = 2;

  auto f = ScalarField<Grid>(lMax);

  auto g = ScalarFieldView(f);

  g(1, 3) = 2;

  /*
  auto x = std::vector<double>(10, 1);
  auto y = std::vector<double>(10, 2);

  auto z = boost::combine(x, y) | boost::adaptors::transformed([](auto pair) {
             auto a = boost::get<0>(pair);
             auto b = boost::get<1>(pair);
             return a + b;
           });

  for (auto value : z) std::cout << value << std::endl;

  */
}
