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
  using Real = double;
  using Complex = std::complex<double>;

  auto start_time = std::chrono::high_resolution_clock::now();

  auto lMax = 512 * 2;
  auto nMax = 2;

  auto grid = GaussLegendreGrid<Real, C2C>(lMax, nMax);

  auto end_time = std::chrono::high_resolution_clock::now();

  std::chrono::duration<Real> duration = end_time - start_time;

  std::cout << duration.count() << std::endl;

  /*
auto f = [](Real theta, Real phi) -> Real { return 1.0; };

auto integral = plan.Integrate(f);

//  std::cout << integral / (4 * std::numbers::pi) << std::endl;

auto nPhi = plan.NumberOfLongitudes();
auto nTheta = plan.NumberOfCoLatitudes();

std::vector<Complex> x(nTheta * nPhi);
std::vector<Complex> y((lMax + 1) * (lMax + 1));

plan.Execute(x, y, 0);
  */
}
