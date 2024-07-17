
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
  int lMax = 2;
  int mMax = 2;
  int nMax = 0;

  auto theta = std::vector<Real>(4, 0.1);

  auto d = WignerNew<double, Ortho, All, All, Single, ColumnMajor>(lMax, mMax,
                                                                   nMax, theta);

  /*
  auto d = Wigner<double, Ortho, All, Single, Single, ColumnMajor>(lMax, mMax,
                                                                   nMax, theta);

  for (auto l : d.Degrees()) {
    for (auto m : d(l).Orders()) {
      std::cout << l << " " << m << " " << d(l)(m) << std::endl;
    }
  }

*/
}
