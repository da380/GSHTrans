

#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

int main() {
  using namespace GSHTrans;

  using Real = double;

  int lMax = 40;
  int mMax = lMax;
  int nMax = lMax;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Real> dist{static_cast<Real>(0),
                                            std::numbers::pi_v<Real>};
  auto theta = dist(gen);
  auto d = Wigner<Real, All, All, FourPi>(lMax, lMax, lMax, theta);

  constexpr auto eps = 1000 * std::numeric_limits<Real>::epsilon();

  for (auto n : d.UpperIndices()) {
    auto d1 = d(n, 0);
    for (auto np : d.UpperIndices()) {
      auto d2 = d(np, 0);
      auto lMin = std::max(std::abs(n), std::abs(np));
      for (auto l = lMin; l <= lMax; l++) {
        auto sum = std::inner_product(d1(l).begin(), d1(l).end(), d2(l).begin(),
                                      Real{0});
        if (n == np) --sum;
        std::cout << sum << std::endl;
      }
    }
  }
}
