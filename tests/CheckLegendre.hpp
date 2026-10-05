#pragma once

#include <GSHTrans/GSHTrans.hpp>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

#include "TestRandom.hpp"

// The upper-index-zero row against the standard library's spherical Legendre
// function, which is an independent implementation of the same thing.
//
// Three details keep the check honest. The error is relative to the largest
// value at the degree, so that small entries are neither skipped nor held to
// an absolute bound. The standard function takes unsigned arguments, so a
// negative order passed straight through would arrive as an enormous degree
// and come back zero; negative orders are instead checked through
// X_{l,-m} = (-1)^m X_{lm}. And the angles are seeded and reported, so that
// a failure can be reproduced.
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
