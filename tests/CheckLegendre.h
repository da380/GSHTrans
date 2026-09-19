#ifndef CHECK_LEGENDRE_GUARD
#define CHECK_LEGENDRE_GUARD

#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

#include "TestRandom.h"

// The upper-index-zero row against the standard library's spherical Legendre
// function, which is an independent implementation of the same thing.
//
// Three things this once got wrong, all of them making it easier to pass.
// The "relative" error divided by `norm`, which was declared inside the
// condition of an `if` and so was the *bool* of the comparison -- one -- and
// the check was absolute. The standard function takes unsigned arguments, so
// every negative order reached it as an enormous degree, came back zero, fell
// under the size threshold and was skipped: half the table was never looked
// at. And the angle came from std::random_device, so a failure could not have
// been reproduced. The error is now measured against the largest value at the
// degree, negative orders are checked through X_{l,-m} = (-1)^m X_{lm}, and
// the angles are seeded and reported.
template <std::floating_point Real>
int CheckLegendre() {
  using namespace GSHTrans;

  const int lMax = 300;
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);
  auto dist = std::uniform_real_distribution<Real>{0, std::numbers::pi_v<Real>};

  // A few drawn, and the two ends of the range a draw will not find.
  auto angles = std::vector<Real>{static_cast<Real>(1e-3),
                                  std::numbers::pi_v<Real> - Real(1e-3)};
  for (auto i = 0; i < 3; i++) angles.push_back(dist(gen));

  constexpr auto eps = 100000 * std::numeric_limits<Real>::epsilon();

  for (auto theta : angles) {
    const auto d = Wigner(lMax, lMax, 0, theta);
    for (auto l : d.Degrees(0)) {
      auto expected = std::vector<Real>{};
      auto scale = Real{0};
      for (auto m : d[l].Orders()) {
        const auto magnitude =
            std::sph_legendre(static_cast<unsigned>(l),
                              static_cast<unsigned>(std::abs(m)), theta);
        expected.push_back(m < 0 && (m % 2 != 0) ? -magnitude : magnitude);
        scale = std::max(scale, std::abs(magnitude));
      }

      auto k = std::size_t{0};
      for (auto m : d[l].Orders()) {
        const auto difference = std::abs(d[l][m] - expected[k++]);
        if (!(difference <= eps * scale)) {
          std::cerr << "CheckLegendre: l = " << l << ", m = " << m
                    << ", theta = " << theta << ", off by " << difference
                    << " against a scale of " << scale << "; seed " << seed
                    << "\n";
          return 1;
        }
      }
    }
  }

  return 0;
}

#endif
