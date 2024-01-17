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

  auto start_time = std::chrono::high_resolution_clock::now();

  auto lMax = 8;
  auto nMax = 2;

  auto grid = GaussLegendreGrid<Real, C2C>(lMax, nMax);

  auto end_time = std::chrono::high_resolution_clock::now();

  std::chrono::duration<Real> duration = end_time - start_time;

  //  std::cout << duration.count() << std::endl;

  auto f = grid.Interpolate([](Real theta, Real phi) -> Complex {
    constexpr auto ii = std::complex<Real>(0, 1);
    using std::sin;
    using std::cos;
    using std::exp;
    return sin(theta) * (5 * cos(theta) * cos(theta) - 1) * exp(ii * phi);
  });
  auto n = 0;
  auto flm = grid.CoefficientVector(n);

  grid.ForwardTransformation(n, f.begin(), flm.begin());

  auto g = grid.FunctionVector<Complex>();

  grid.InverseTransformation(lMax, n, flm.begin(), g.begin());

  auto fIter = f.begin();
  auto gIter = g.begin();
  while (fIter != f.end()) std::cout << *fIter++ - *gIter++ << std::endl;
}
