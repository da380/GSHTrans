

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

  auto size1 = WignerNStorage<All>(lMax, mMax, n);
  std::vector<double> data1(size1);

  auto size2 = WignerNStorage<NonNegative>(lMax, mMax, n);
  std::vector<double> data2(size2);

  double theta = 2.1;

  WignerN<double, All> p1(data1, lMax, mMax, n, theta);
  std::cout << p1(4, 2) << std::endl;

  WignerN<double, NonNegative> p2(data2, lMax, mMax, n, theta);
  std::cout << p2(4, 2) << std::endl;

  std::cout << p1(4, 2) << " " << p2(4, 2) << std::endl;
}
