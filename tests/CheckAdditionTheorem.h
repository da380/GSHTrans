#ifndef CHECK_ADDITION_THEOREM_ZERO_GUARD
#define CHECK_ADDITION_THEOREM_ZERO_GUARD

#include <Eigen/Core>
#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

template <std::floating_point Real>
int CheckAdditionTheorem() {
  using namespace GSHTrans;

  // Set the degree, order and upper index
  int lMax = 40;
  int mMax = lMax;
  int nMax = lMax;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Real> dist1{static_cast<Real>(0),
                                             std::numbers::pi_v<Real>};
  auto theta = dist1(gen);

  auto d = Wigner<Real, All, All, FourPi>(lMax, lMax, lMax, theta);

  constexpr auto eps = 1000 * std::numeric_limits<Real>::epsilon();

  for (auto n : d.UpperIndices()) {
    auto d1 = d(n)();
    for (auto np : d.UpperIndices()) {
      auto d2 = d(np)();
      auto lMin = std::max(std::abs(n), std::abs(np));
      for (auto l = lMin; l <= lMax; l++) {
        auto sum = std::inner_product(d1(l).begin(), d1(l).end(), d2(l).begin(),
                                      Real{0});
        if (n == np) --sum;
        if (std::abs(sum) > eps) return 1;
      }
    }
  }

  return 0;
}

#endif
