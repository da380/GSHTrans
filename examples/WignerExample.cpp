
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
  constexpr int lMax = 256;
  constexpr int mMax = lMax;
  constexpr int nMax = lMax;

  constexpr auto theta = std::array<Real, 1>{0.1};

  auto d = Wigner<Real, Ortho, All, All, Single, ColumnMajor>(lMax, mMax, nMax,
                                                              theta);

  /*
  for (auto l : d.Degrees(0)) {
    for (auto m : d[l].Orders()) {
      std::cout << l << " " << m << " " << d[l][m] << std::endl;
    }
  }
*/
}
