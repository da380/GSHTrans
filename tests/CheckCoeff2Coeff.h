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

#include "TestRandom.h"

using namespace GSHTrans;

using Int = std::ptrdiff_t;

// Round trip: random coefficients -> field -> coefficients. Returns true on
// failure. The generator is supplied by the caller so that a failing run can
// be reproduced from the seed it reports.
template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange,
          IndexRange NRange>
auto Coeff2Coeff(GSHTransTest::Generator& gen) {
  using Real = RemoveComplex<Scalar>;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  auto lMaxGrid = GSHTransTest::RandomDegree(gen, 4, 256);
  auto lMax = GSHTransTest::RandomDegree(gen, 4, lMaxGrid);
  auto nMax = std::min(lMax, Int(4));
  auto grid = Grid(lMaxGrid, nMax);

  // Real-valued fields exist only at upper index zero.
  auto n = RealFloatingPoint<Scalar>
               ? Int{0}
               : GSHTransTest::RandomUpperIndex<NRange>(gen, nMax);

  auto getSize = [](Int lMax, Int n) {
    if constexpr (RealFloatingPoint<Scalar>) {
      return GSHIndices<NonNegative>(lMax, lMax, n).Size();
    } else {
      return GSHIndices<All>(lMax, lMax, n).Size();
    }
  };

  auto size = getSize(lMax, n);
  auto flm = FFTWpp::vector<Complex>(size);

  if constexpr (ComplexFloatingPoint<Scalar>) {
    GSHTransTest::RandomComplexCoefficient(grid, lMax, n, flm, gen);
  } else {
    GSHTransTest::RandomRealCoefficient(grid, lMax, flm, gen);
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
