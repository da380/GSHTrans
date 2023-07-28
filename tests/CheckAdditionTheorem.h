#ifndef CHECK_ADDITION_THEOREM_ZERO_GUARD
#define CHECK_ADDITION_THEOREM_ZERO_GUARD

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>
#include <random>

#include "GSHTrans/Core"

template <std::floating_point Float>
int CheckAdditionTheorem() {
  using namespace GSHTrans;

  // Set the degree, order and upper index
  int L = 100;
  int M = L;
  int N = L;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Float> dist{0., std::numbers::pi_v<Float>};
  auto theta = dist(gen);

  // Construct the unnormalised Wigner values
  Wigner d(L, M, N, theta);

  // Define small number for comparison.
  constexpr auto eps = 1000 * std::numeric_limits<Float>::epsilon();

  for (int l : d.Degrees()) {
    for (int n : d.UpperIndices(l, true)) {
      auto start = d.cbegin(n, l);
      auto finish = d.cend(n, l);
      auto sum = std::accumulate(start, finish, static_cast<Float>(0),
                                 [](Float x, Float y) { return x + y * y; });
      auto diff = std::abs(sum - static_cast<Float>(1));
      if (diff > eps) {
        return 1;
      }
    }
  }
  return 0;
}

#endif
