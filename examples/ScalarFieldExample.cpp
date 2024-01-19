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
  using Type = R2C;

  /*

  using Grid = GaussLegendreGrid<Real, Type>;

  using Scalar = Grid::scalar_type;

  auto lMax = 4;
  auto nMax = 2;

  auto grid = std::make_shared<Grid>(lMax, nMax);

  auto f = ScalarField(grid);

  auto data = std::vector<Scalar>(grid->FieldDimension());
  auto g = ScalarFieldView(grid, data);

  auto h = ScalarFieldUnaryFunction(f, [](auto x) { return x + 1; });

  auto k =
      ScalarFieldBinaryFunction(f, g, [](auto x, auto y) { return x + y; });

  std::cout << f(1, 1) << " " << g(1, 1) << " " << h(1, 1) << " " << k(1, 1)
            << std::endl;

  */

  /*

  auto h = ScalarFieldUnaryFunction(f, [](auto x) { return x + 1; });

  auto i = ScalarField<Grid>(h);

  auto j = ScalarFieldUnaryFunction(h, [](auto x) { return -x; });

  auto k =
      ScalarFieldBinaryFunction(i, j, [](auto x, auto y) { return x + y; });

  std::cout << f(1, 1) << " " << h(1, 1) << " " << j(1, 1) << " " << k(1, 1)
            << std::endl;

  */

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
