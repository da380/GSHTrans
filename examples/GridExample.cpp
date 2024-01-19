#include <FFTWpp/All>
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

using namespace GSHTrans;

using Int = std::ptrdiff_t;

auto RandomDegree() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> d(3, 7);
  return std::pow(2, d(gen));
}

auto RandomUpperIndex(Int nMin, Int nMax) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> d(nMin, nMax);
  return d(gen);
}

int main() {
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Complex;
  using MRange = All;
  using NRange = All;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = RandomDegree();
  auto nMax = 4;
  auto grid = Grid(lMax, nMax);
  auto n = RandomUpperIndex(-nMax, nMax);

  // Make a random coefficient.
  auto size = GSHIndices<All>(lMax, lMax, n).size();
  auto flm = FFTWpp::vector<Complex>(size);
  {
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<Real> d{0., 1.};
    std::ranges::generate(flm, [&gen, &d]() {
      return Complex{d(gen), d(gen)};
    });
    auto flmView = GSHView<Complex, MRange>(lMax, lMax, n, flm.begin());
    flmView(lMax)(lMax) = 0;
  }

  auto f = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
                                  grid.NumberOfCoLatitudes());

  auto glm = FFTWpp::vector<Complex>(size);

  grid.InverseTransformation(lMax, n, flm.begin(), f.begin());
  grid.ForwardTransformation(lMax, n, f.begin(), glm.begin());

  std::ranges::transform(flm, glm, flm.begin(),
                         [](auto f, auto g) { return f - g; });

  auto err = *std::ranges::max_element(
      flm, [](auto a, auto b) { return std::abs(a) < std::abs(b); });

  std::cout << std::abs(err) << std::endl;
}
