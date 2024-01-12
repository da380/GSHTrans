

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

  // Set the degree, order and upper index
  int lMax = 3;
  int mMax = lMax;
  int nMax = lMax;

  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  constexpr auto pi = std::numbers::pi_v<Real>;
  auto dist1 = std::uniform_real_distribution<Real>(static_cast<Real>(0), pi);
  const auto theta1 = dist1(gen);
  auto dist2 = std::uniform_real_distribution<Real>(theta1, pi);
  const auto theta3 = dist2(gen);
  const auto theta2 = theta3 - theta1;

  auto d = Wigner<Real, All, All, FourPi>(lMax, lMax, lMax,
                                          std::vector{theta1, theta2, theta3});

  auto mat = d.Matrix(0, 2);

  for (auto l : d.Degrees()) {
    auto m1 = d.Matrix(0, l);
    auto m2 = d.Matrix(1, l);
    auto m3 = d.Matrix(2, l);
    auto err = (m3 - m2 * m1).lpNorm<Eigen::Infinity>();
    std::cout << err << std::endl;
  }

  /*



  constexpr auto eps = 1000 * std::numeric_limits<Real>::epsilon();

  for (auto n : d.UpperIndices()) {
    auto d1 = d(n)();
    for (auto np : d.UpperIndices()) {
      auto d2 = d(np)();
      auto lMin = std::max(std::abs(n), std::abs(np));
      for (auto l = lMin; l <= lMax; l++) {
        auto sum = std::inner_product(d1(l).begin(), d1(l).end(), d2(l).begin(),
                                      Real{0});
        if (n == np) --sum;
        std::cout << sum << std::endl;
      }
    }
  }

  */
}
