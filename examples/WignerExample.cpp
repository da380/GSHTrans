

#include <Eigen/Core>
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
  using Float = double;

  int lMax = 10;
  int mMax = lMax;
  int nMax = lMax;
  Float theta = 0.2;

  // Construct the Wigner array.
  WignerArray<Float, All, All> d(lMax, mMax, nMax, theta);

  std::cout << d(3, -2, -2) << std::endl;

  auto D = WignerMatrix(3, theta, Normalisation::Ortho);
  std::cout << D << std::endl;
}
