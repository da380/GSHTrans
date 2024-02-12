#include <FFTWpp/All>
#include <GSHTrans/All>
#include <iostream>
#include <ranges>
#include <vector>

int main() {
  using namespace GSHTrans;

  using Int = std::ptrdiff_t;
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Real;
  using Vector = FFTWpp::vector<Scalar>;
  using MRange = All;

  Int l = 20;
  Int mMax = 4;

  auto Size = [](auto l, auto mMax) {
    if constexpr (std::same_as<MRange, All>) {
      return 2 * std::min(l, mMax) + 1;
    } else {
      return std::min(l, mMax) + 1;
    }
  };

  auto size = Size(l, mMax);
  auto data = Vector(size);
  auto view = CoefficientSubView(l, mMax, data, MRange{});
}