#ifndef CHECK_COEFF_2_COEFF_GUARD_H
#define CHECK_COEFF_2_COEFF_GUARD_H

#include <FFTWpp/Core>
#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <concepts>
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

template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange,
          IndexRange NRange>
auto Coeff2Coeff() {
  using Real = FFTWpp::RemoveComplex<Scalar>;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMaxGrid = RandomDegree(4, 256);
  auto lMax = RandomDegree(4, lMaxGrid);
  auto nMax = std::min(lMax, Int(4));
  auto grid = Grid(lMaxGrid, nMax);

  auto n = RandomUpperIndex<NRange>(nMax);

  auto getSize = [](Int lMax, Int n) {
    if constexpr (RealFloatingPoint<Scalar>) {
      return GSHIndices<NonNegative>(lMax, lMax, n).size();
    } else {
      return GSHIndices<All>(lMax, lMax, n).size();
    }
  };

  auto size = getSize(lMax, n);
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

  std::ranges::transform(flm, glm, flm.begin(),
                         [](auto f, auto g) { return f - g; });

  return std::ranges::any_of(flm, [](auto f) {
    constexpr auto eps = 50000 * std::numeric_limits<Real>::epsilon();
    return std::abs(f) > eps;
  });
}
#endif  // CHECK_COEFF_2_COEFF_GUARD_H
