#include <GSHTrans/All>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

int main() {
  using namespace GSHTrans;

  auto start_time = std::chrono::high_resolution_clock::now();

  Plan<double, C2C> plan(256, 2);

  auto end_time = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> duration = end_time - start_time;

  std::cout << duration.count() << std::endl;

  auto f = [](auto theta, auto phi) { return 1; };

  auto integral = plan.Integrate(f);

  std::cout << integral / (4 * std::numbers::pi) << std::endl;
}
