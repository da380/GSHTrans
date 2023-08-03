#ifndef CHECK_ADDITION_THEOREM_ZERO_GUARD
#define CHECK_ADDITION_THEOREM_ZERO_GUARD

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

#include "GSHTrans/Core"

template <std::floating_point Float>
int CheckAdditionTheorem() {
  using namespace GSHTrans;

  // Set the degree, order and upper index
  int L = 50;
  int M = L;
  int N = L;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Float> dist{static_cast<Float>(0),
                                             std::numbers::pi_v<Float>};
  auto theta = dist(gen);

  using Wigner = GSHTrans::Wigner<Float, AllOrders, FourPiNormalised>;
  std::vector<std::unique_ptr<Wigner>> d(2 * L + 1);

  for (int n = -L; n <= L; n++) {
    d.push_back(GSHTrans::MakeWigner<Float, AllOrders, FourPiNormalised>(
        L, M, n, theta));
  }

  /*

  // Construct the unnormalised Wigner values
  Wigner d(L, M, N, theta);

  // Define small number for comparison.
  constexpr auto eps = 1000 * std::numeric_limits<Float>::epsilon();

  for (int l : d.Degrees()) {
    for (int n : d.UpperIndices(l)) {
      for (int np : d.UpperIndices(l)) {
        Float sum = static_cast<Float>(0);
        for (int m : d.Orders(l)) {
          sum += d(l, m, n) * d(l, m, np);
        }
        if (n == np) {
          if (std::abs(sum - static_cast<Float>(1)) > eps) return 1;
        } else {
          if (abs(sum) > eps) return 1;
        }
      }
    }
  }

  */
  return 0;
}

#endif
