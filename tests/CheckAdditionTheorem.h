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

template <std::floating_point Float>
int CheckAdditionTheorem() {
  using namespace GSHTrans;

  // Set the degree, order and upper index
  int L = 30;
  int M = L;
  int N = L;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Float> dist{static_cast<Float>(0),
                                             std::numbers::pi_v<Float>};
  auto theta = dist(gen);

  // Define small number for comparison.
  constexpr auto eps = 1000 * std::numeric_limits<Float>::epsilon();

  for (int n = -N; n <= N; n++) {
    WignerN<Float, AllOrders> d1(L, M, n, theta, Normalisation::FourPi);
    for (int np = 0; np <= N; np++) {
      WignerN<Float, AllOrders> d2(L, M, np, theta, Normalisation::FourPi);

      auto lstart = std::max(std::abs(n), std::abs(np));
      for (int l = lstart; l <= L; l++) {
        Float sum = 0.0;
        for (int m = -l; m <= l; m++) {
          sum += d1(l, m) * d2(l, m);
        }
        if (n == np) --sum;
        if (std::abs(sum) > eps) return 1;
      }
    }
  }

  return 0;
}

#endif
