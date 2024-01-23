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

Int RandomDegree() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> d(4, 128);
  return d(gen);
}

template <IndexRange NRange>
Int RandomUpperIndex(Int nMax) {
  std::random_device rd;
  std::mt19937 gen(rd());
  if constexpr (std::same_as<NRange, All>) {
    std::uniform_int_distribution<Int> d(-nMax, nMax);
    return d(gen);
  } else {
    std::uniform_int_distribution<Int> d(0, nMax);
    return d(gen);
  }
}

int main() {
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Real;
  using MRange = All;
  using NRange = All;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMax = RandomDegree();
  auto nMax = std::min(lMax, Int(4));
  auto grid = Grid(lMax, nMax);
  auto n = RandomUpperIndex<NRange>(nMax);

  auto size = grid.CoefficientSize<Scalar>(lMax, n);
  auto flm = FFTWpp::vector<Complex>(size);
  if constexpr (ComplexFloatingPoint<Scalar>) {
    grid.RandomComplexCoefficient(lMax, n, flm);
  } else {
    grid.RandomRealCoefficient(lMax, n, flm);
  }

  auto f = FFTWpp::vector<Scalar>(grid.ComponentSize());
  auto glm = FFTWpp::vector<Complex>(size);

  grid.InverseTransformation(lMax, n, flm, f);
  grid.ForwardTransformation(lMax, n, f, glm);

  std::ranges::transform(flm, glm, flm.begin(),
                         [](auto f, auto g) { return f - g; });

  auto err = *std::ranges::max_element(
      flm, [](auto a, auto b) { return std::abs(a) < std::abs(b); });

  std::cout << lMax << " " << n << " " << std::abs(err) << std::endl;
}
