

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

  auto size = WignerNStorage<All>(lMax, mMax, n);

  std::vector<double> data1(size);
  std::vector<double> data2(size);

  double theta = 2.1;

  auto p1 = WignerN<double, All>(data2, lMax, mMax, n, theta);

  auto p2 = WignerN<double, All>(data2, lMax, mMax, n, theta);

  auto p3 = WignerN<double, All>(lMax, mMax, n, theta);

  std::cout << p1(4, 2) << " " << p2(4, 2) << " " << p3(4, 2) << std::endl;
}
