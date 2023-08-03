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
  int L = 10;
  int M = L;
  int N = 0;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Float> dist{static_cast<Float>(0),
                                             std::numbers::pi_v<Float>};
  auto theta = dist(gen);

  // Define small number for comparison.
  constexpr auto eps = 1000 * std::numeric_limits<Float>::epsilon();

  for (int n = -N; n <= N; n++) {
    Wigner<Float, AllOrders> d1(L, M, n, theta, Normalisation::FourPi);
    for (int np = -N; np <= N; np++) {
      Wigner<Float, AllOrders> d2(L, M, np, theta, Normalisation::FourPi);

      auto lstart = std::max(std::abs(n), std::abs(np));
      for (int l = lstart; l <= L; l++) {
        Float sum = std::inner_product(d1.begin(l), d1.end(l), d2.begin(l),
                                       static_cast<Float>(0.0));
        //        if (n == np) --sum;
        std::cout << n << " " << np << " " << l << " " << sum << std::endl;
      }
    }
  }

  return 0;
}

#endif
