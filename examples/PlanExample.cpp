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

  Plan<double, R2C> plan(10, 0);

  auto f = [](auto theta, auto phi) { return 1; };

  auto integral = plan.Integrate(f);

  std::cout << integral / (4 * std::numbers::pi) << std::endl;

  for (auto theta : plan.CoLatitudes()) std::cout << theta << std::endl;
}
