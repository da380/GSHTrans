#include <FFTWpp/All>
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

using namespace GSHTrans;

using Int = std::ptrdiff_t;

auto RandomDegree() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> d(3, 7);
  return std::pow(2, d(gen));
}

auto RandomUpperIndex(Int nMin, Int nMax) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> d(nMin, nMax);
  return d(gen);
}

int main() {
  using Real = double;
  using Complex = std::complex<Real>;
  using MRange = All;
  using NRange = All;
  using Scalar = Complex;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  constexpr auto realTransform = RealFloatingPoint<Scalar>;
  constexpr auto complexTransform = ComplexFloatingPoint<Scalar>;

  //  auto lMax = RandomDegree();
  auto lMax = 4;
  auto nMax = 0;
  auto grid = Grid(lMax, nMax);
  //  auto n = RandomUpperIndex(grid.MinUpperIndex(), nMax);
  auto n = 0;

  // Make a random coefficient.
  auto size = Int();
  if constexpr (realTransform) {
    size = GSHIndices<NonNegative>(lMax, lMax, n).size();
  } else {
    size = GSHIndices<All>(lMax, lMax, n).size();
  }
  auto flm = FFTWpp::vector<Complex>(size);
  {
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<Real> d{0., 1.};
    std::ranges::generate(flm, [&gen, &d]() {
      return Complex{d(gen), d(gen)};
    });
    auto flmView = GSHView<Complex, MRange>(lMax, lMax, n, flm.begin());
    if constexpr (complexTransform) {
      flmView(lMax)(lMax) = 0;
    } else {
      for (auto l : flmView.Degrees()) {
        flmView(l)(0) = 0;
      }
      flmView(lMax)(lMax).imag(0);
    }
  }

  // Allocate other vectors
  auto f = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
                                  grid.NumberOfCoLatitudes());

  auto glm = FFTWpp::vector<Complex>(size);

  grid.InverseTransformation(n, flm.begin(), f.begin());

  //  for (auto val : flm) std::cout << val << std::endl;

  /*

  grid.ForwardTransformation(n, f.begin(), glm.begin());



  std::ranges::transform(flm, glm, flm.begin(),
                         [](auto f, auto g) { return f - g; });

  auto check = std::ranges::all_of(flm, [](auto f) {
    constexpr auto eps = 10000 * std::numeric_limits<Real>::epsilon();
    return std::abs(f) < eps;
  });

  */
}
