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
  using MRange = NonNegative;

  Int lMax = 5;
  Int mMax = 3;
  Int n = 4;

  auto indices = Testing::GSHIndices<MRange>(lMax, mMax, n);

  auto count = 0;
  for (auto l : indices.Degrees()) {
    auto [offset, subindices] = indices.Index(l);
    for (auto m : subindices.Orders()) {
      std::cout << l << " " << m << " "
                << offset + subindices.Index(m) - count++ << std::endl;
    }
  }
}