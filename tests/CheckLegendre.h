#ifndef CHECK_LEGENDRE_GUARD
#define CHECK_LEGENDRE_GUARD

#include <GSHTrans/All>
#include <cmath>
#include <concepts>
#include <execution>
#include <limits>
#include <numbers>
#include <random>

template <std::floating_point Real>
int CheckLegendre() {
  using namespace GSHTrans;

  // Set the degree, order and upper index
  int lMax = 300;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Real> dist{0., std::numbers::pi_v<Real>};
  auto theta = dist(gen);

  // Construct the normalised Wigner values
  auto d = Testing::Wigner(lMax, lMax, 0, theta);

  // Define small numbers for comparison.
  constexpr auto eps = 100000 * std::numeric_limits<Real>::epsilon();
  constexpr auto tiny = 1000 * std::numeric_limits<Real>::min();

  // Compare values to std library function
  for (auto l : d.Degrees()) {
    for (auto m : d(l).Orders()) {
      Real plm = d(l)(m);
      Real plmSTD = std::sph_legendre(l, m, theta);
      if (auto norm = std::abs(plmSTD) > tiny) {
        Real diff = std::abs(plm - plmSTD) / norm;
        if (diff > eps) return 1;
      }
    }
  }

  return 0;
}

#endif
