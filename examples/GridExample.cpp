
#include <GSHTrans/All>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
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
  using Scalar = Complex;
  using MRange = All;
  using NRange = All;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  {
    auto lMaxGrid = RandomDegree(4, 256);
    auto lMax = RandomDegree(4, lMaxGrid);
    auto nMax = std::min(lMax, Int(4));

    auto grid = Grid(lMaxGrid, nMax);
    auto n = RandomUpperIndex<NRange>(nMax);

    Int size;
    if constexpr (ComplexFloatingPoint<Scalar>) {
      size = grid.CoefficientSize(lMax, n);
    } else {
      size = grid.CoefficientSizeNonNegative(lMax, n);
    }
    auto flm = FFTWpp::vector<Complex>(size);
    if constexpr (ComplexFloatingPoint<Scalar>) {
      grid.RandomComplexCoefficient(lMax, n, flm);
    } else {
      grid.RandomRealCoefficient(lMax, n, flm);
    }

    auto f = FFTWpp::vector<Scalar>(grid.FieldSize());
    auto glm = FFTWpp::vector<Complex>(size);

    grid.InverseTransformation(lMax, n, flm, f);

    grid.ForwardTransformation(lMax, n, f, glm);

    auto error = std::ranges::max(std::ranges::views::zip_transform(
        [](auto x, auto y) { return std::abs(x - y); }, flm, glm));

    ;

    std::cout << lMaxGrid << " " << lMax << " " << n << " " << error
              << std::endl;
  }

  FFTWpp::CleanUp();
}
