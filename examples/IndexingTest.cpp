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

  Int lMax = 5;
  Int mMax = 3;
  Int n = 1;
}