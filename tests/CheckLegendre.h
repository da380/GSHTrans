#ifndef CHECK_LEGENDRE_GUARD
#define CHECK_LEGENDRE_GUARD

#include <GSHTrans/All>
#include <cmath>
#include <concepts>
#include <execution>
#include <limits>
#include <numbers>
#include <random>

template <std::floating_point Float>
int CheckLegendre() {
  using namespace GSHTrans;

  // Set the degree, order and upper index
  int lMax = 300;
  int mMax = lMax;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Float> dist{0., std::numbers::pi_v<Float>};
  auto theta = dist(gen);

  // Construct the normalised Wigner values
  LegendreArray<Float, NonNegative> d(lMax, mMax, theta);

  // Define small numbers for comparison.
  constexpr auto eps = 100000 * std::numeric_limits<Float>::epsilon();
  constexpr auto tiny = 1000 * std::numeric_limits<Float>::min();

  // Compare values to std library function
  for (int l = 0; l <= lMax; l++) {
    for (int m = 0; m <= l; m++) {
      Float plm = d(l, m);
      Float plmSTD = std::sph_legendre(l, m, theta);
      if (auto norm = std::abs(plmSTD) > tiny) {
        Float diff = std::abs(plm - plmSTD) / norm;
        if (diff > eps) return 1;
      }
    }
  }

  return 0;
}

#endif
