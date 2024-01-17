#include <GSHTrans/All>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>
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
  /*

  using namespace GSHTrans;
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, R2C>;
  using Scalar = Grid::scalar_type;

  auto lMax = 4;
  auto nMax = 2;

  auto f = ScalarField<Grid>(lMax);

  for (auto& value : f) value = 2;

  for (auto value : f) std::cout << value << std::endl;

  */

  auto x = std::vector<double>(10, 1);
  auto y = std::vector<double>(10, 2);

  auto z = boost::combine(x, y) | boost::adaptors::transformed([](auto pair) {
             return boost::get<0>(pair) + boost::get<1>(pair);
           });

  for (auto value : z) std::cout << value << std::endl;
}
