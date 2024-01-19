#ifndef CHECK_COEFF_2_COEFF_GUARD_H
#define CHECK_COEFF_2_COEFF_GUARD_H

#include <FFTWpp/All>
#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <execution>
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

auto RandomUpperIndex(Int nMax) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> d(0, nMax);
  return d(gen);
}

template <RealFloatingPoint Real, OrderIndexRange MRange, IndexRange NRange,
          RealOrComplexFloatingPoint Scalar>
auto Coeff2Coeff() {
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  constexpr auto realTransform = RealFloatingPoint<Scalar>;
  constexpr auto complexTransform = ComplexFloatingPoint<Scalar>;

  auto lMax = RandomDegree();
  auto nMax = 4;

  auto n = RandomUpperIndex(nMax);
  auto grid = Grid(lMax, nMax);

  // Make a random coefficient.
  auto size = Int();
  if constexpr (realTransform) {
    size = GSHIndices<NonNegative>(lMax, lMax, n).size();
  } else {
    size = GSHIndices<All>(lMax, lMax, n).size();
  }
  auto flm = FFTWpp::vector<Complex>(size);
  {
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<Real> d{0., 1.};
    std::ranges::generate(flm, [&gen, &d]() {
      return Complex{d(gen), d(gen)};
    });
    auto flmView = GSHView<Complex, MRange>(lMax, lMax, n, flm.begin());
    if constexpr (complexTransform) {
      flmView(lMax)(lMax) = 0;
    } else {
      for (auto l : flmView.Degrees()) {
        flmView(l)(0) = 0;
      }
      flmView(lMax)(lMax).imag(0);
    }
  }

  // Allocate other vectors
  auto f = FFTWpp::vector<Scalar>(grid.NumberOfCoLatitudes() *
                                  grid.NumberOfLongitudes());
  auto glm = FFTWpp::vector<Complex>(size);

  grid.InverseTransformation(n, flm.begin(), f.begin());
  grid.ForwardTransformation(n, f.begin(), glm.begin());

  std::ranges::transform(flm, glm, flm.begin(),
                         [](auto f, auto g) { return f - g; });

  return std::ranges::any_of(flm, [](auto f) {
    constexpr auto eps = 10000 * std::numeric_limits<Real>::epsilon();
    return std::abs(f) < eps;
  });
}
#endif  // CHECK_COEFF_2_COEFF_GUARD_H
