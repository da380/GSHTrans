#pragma once

#include <GSHTrans/GSHTrans.hpp>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <random>

#include "TestRandom.hpp"

template <std::floating_point Real>
int CheckAdditionTheorem() {
  using namespace GSHTrans;

  // Set the degree, order and upper index
  int lMax = 40;

  // A seeded angle, so that a failure can be run again: see TestRandom.hpp.
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);
  std::uniform_real_distribution<Real> dist1{static_cast<Real>(0),
                                             std::numbers::pi_v<Real>};
  auto theta = dist1(gen);

  auto d = Wigner<Real, All, All>(lMax, lMax, lMax, theta);

  constexpr auto eps = 1000 * std::numeric_limits<Real>::epsilon();

  for (auto n : d.UpperIndices()) {
    auto d1 = d[n];
    for (auto np : d.UpperIndices()) {
      auto d2 = d[np];
      auto lMin = std::max(std::abs(n), std::abs(np));
      for (auto l = lMin; l <= lMax; l++) {
        // The stored values carry the orthonormal factor sqrt((2l+1)/(4 pi)),
        // the only normalisation the library offers, so
        // sum_m d^l_{nm} d^l_{n'm} = delta_{nn'} (2l+1)/(4 pi). Undo that
        // factor and the check is the addition theorem as written.
        const auto normalisation =
            4 * std::numbers::pi_v<Real> / static_cast<Real>(2 * l + 1);
        auto sum =
            normalisation * std::inner_product(d1[l].begin(), d1[l].end(),
                                               d2[l].begin(), Real{0});
        if (n == np) --sum;
        if (std::abs(sum) > eps) return 1;
      }
    }
  }

  return 0;
}
