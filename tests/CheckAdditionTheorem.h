#ifndef CHECK_ADDITION_THEOREM_ZERO_GUARD
#define CHECK_ADDITION_THEOREM_ZERO_GUARD

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
  std::uniform_real_distribution<Real> dist{static_cast<Real>(0),
                                            std::numbers::pi_v<Real>};
  auto theta = dist(gen);

  // Define small number for comparison.
  constexpr auto eps = 1000 * std::numeric_limits<Real>::epsilon();
  for (auto n = -nMax; n <= nMax; n++) {
    Wigner<Real, All, FourPi> d1(lMax, mMax, n, theta);
    for (auto np = 0; np <= nMax; np++) {
      Wigner<Real, All, FourPi> d2(lMax, mMax, np, theta);
      auto lstart = std::max(std::abs(n), std::abs(np));
      for (auto l = lstart; l <= lMax; l++) {
        auto sum = std::inner_product(d1.beginForDegree(l), d1.endForDegree(l),
                                      d2.beginForDegree(l), Real{0});
        if (n == np) --sum;
        if (std::abs(sum) > eps) return 1;
      }
    }
  }
  return 0;
}

#endif
