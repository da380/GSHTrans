
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
  int lMax = 10;
  int mMax = 5;
  int nMax = 2;

  auto theta = double(1);

  auto d = Wigner<double, NonNegative, All, Ortho>(lMax, mMax, mMax, theta);
}
