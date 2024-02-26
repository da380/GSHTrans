
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
  int lMax = 4;
  int mMax = 4;
  int nMax = 0;

  auto theta = double(1);

  auto d = Wigner<double, Ortho, All, Single, Single, ColumnMajor>(lMax, mMax,
                                                                   nMax, theta);

  for (auto l : d.Degrees()) {
    for (auto m : d(l).Orders()) {
      std::cout << l << " " << m << " " << d(l)(m) << std::endl;
    }
  }
}
