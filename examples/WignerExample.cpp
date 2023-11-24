

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

  int lMax = 100;
  int mMax = 50;
  int n = 0;

  double theta = 2.1;

  auto p1 = WignerN<double, All>(lMax, mMax, n, theta);
  auto p2 = WignerN<double, NonNegative>(lMax, mMax, n, theta);

  auto p3 =
      WignerN<double, All>(lMax, mMax, n, std::vector<double>{theta, theta});

  std::cout << p1(4, 2) << " " << p2(4, 2) << " " << p3(0, 4, 2) << std::endl;

  auto range = p1.RangeForAngleAndDegree(0, 2);
  for (auto val : range) std::cout << val << std::endl;
}
