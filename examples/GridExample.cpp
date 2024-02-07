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

Int RandomDegree(Int lMin, Int lMax) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> d(lMin, lMax);
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

  {
    auto lMaxGrid = RandomDegree(4, 256);
    auto lMax = RandomDegree(4, lMaxGrid);
    auto nMax = std::min(lMax, Int(4));
    lMaxGrid = 256;
    lMax = 256;
    nMax = 2;

    auto grid = Grid(lMaxGrid, nMax);
    auto n = RandomUpperIndex<NRange>(nMax);

    Int size;
    if constexpr (ComplexFloatingPoint<Scalar>) {
      size = grid.ComplexCoefficientSize(lMax, n);
    } else {
      size = grid.RealCoefficientSize(lMax, n);
    }
    auto flm = FFTWpp::vector<Complex>(size);
    if constexpr (ComplexFloatingPoint<Scalar>) {
      grid.RandomComplexCoefficient(lMax, n, flm);
    } else {
      grid.RandomRealCoefficient(lMax, n, flm);
    }

    auto f = FFTWpp::vector<Scalar>(grid.ComponentSize());
    auto glm = FFTWpp::vector<Complex>(size);

    grid.InverseTransformation(lMax, n, std::ranges::views::all(flm), f);

    grid.ForwardTransformation(
        lMax, n,
        std::ranges::views::all(f) |
            std::ranges::views::transform([](auto x) { return x; }),
        glm);

    auto errors = std::ranges::views::zip_transform(
        [](auto x, auto y) { return std::abs(x - y); }, flm, glm);

    auto err = *std::ranges::max_element(errors);

    std::cout << lMaxGrid << " " << lMax << " " << n << " " << std::abs(err)
              << std::endl;
  }

  FFTWpp::CleanUp();
}
