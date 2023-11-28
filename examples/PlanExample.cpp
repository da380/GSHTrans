

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

  PlanComplex<double> plan(10, 0);

  auto f = [](auto theta, auto phi) { return 1; };

  auto integral = plan.Integrate(f);

  std::cout << integral / (4 * std::numbers::pi) << std::endl;
}
