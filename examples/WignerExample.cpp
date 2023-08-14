

#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <concepts>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

int main() {
  using namespace GSHTrans;

  using Float = double;

  // Set the degree, order and upper index
  int L = 3;
  int M = L;
  int N = L;

  // Pick a random angle
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_real_distribution<Float> dist{static_cast<Float>(0),
                                             std::numbers::pi_v<Float>};
  Float theta = 0.6;

  std::cout.setf(std::ios_base::scientific);
  std::cout.setf(std::ios_base::showpos);
  std::cout.precision(7);

  std::vector<std::unique_ptr<Wigner<Float, AllOrders>>> d;
  for (int n = -N; n <= N; n++) {
    d.push_back(std::make_unique<Wigner<Float, AllOrders>>(
        L, M, n, theta, Normalisation::FourPi));
  }

  for (int l = L; l <= L; l++) {
    auto nRange = std::min(l, N);
    for (int n = -nRange; n <= -1; n++) {
      std::cout << n << "  " << l << " | ";
      auto mRange = std::min(l, M);
      for (int m = -mRange; m <= mRange; m++) {
        std::cout << d[n + N]->operator()(l, m) << "   ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
