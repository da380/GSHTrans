#ifndef CHECK_WIGNER_GUARD
#define CHECK_WIGNER_GUARD

#include <GSHTrans/All>
#include <cmath>
#include <concepts>
#include <execution>
#include <limits>
#include <numbers>
#include <random>

template <std::floating_point Float>
int CheckWigner() {
  using namespace GSHTrans;

  // Set the degree and order
  int lMax = 200;
  int mMax = lMax;

  // Set a random upper index
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_int_distribution idist{-lMax, lMax};
  int nMax = idist(gen);

  // Pick a random angle
  std::uniform_real_distribution<Float> rdist{0.0, std::numbers::pi_v<Float>};
  auto theta = rdist(gen);

  // Set the normalisation
  Normalisation norm;
  if (idist(gen) > 0) {
    norm = Normalisation::Ortho;
  } else {
    norm = Normalisation::FourPi;
  }

  // Construct the WignerArrayN object
  WignerArrayN d(lMax, mMax, nMax, theta, norm);

  // Define small numbers for comparison.
  constexpr auto eps = 1.0e-8;
  constexpr auto tiny = 1.0e-8;

  // Loop over degrees.
  for (int l = std::abs(nMax); l <= lMax; l++) {
    // Create WignerArrayLN object.
    WignerArrayLN d2(l, nMax, theta, norm);

    // Compare values at each order.
    for (int m = -l; m <= l; m++) {
      Float v1 = d(l, m);
      Float v2 = d2(m);
      Float norm = std::max(std::abs(v1), std::abs(v2));
      if (norm > tiny) {
        Float diff = std::abs(v1 - v2) / norm;
        if (diff > eps) return 1;
      }
    }
  }

  return 0;
}

#endif
