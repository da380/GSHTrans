

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

  std::vector<double> data(size);

  double th = 1.9;

  auto p = WignerN(data.begin(), data.end(), lMax, mMax, n, th);

  //  auto q = WignerNVector(lMax, mMax, n, th);

  //  std::cout << q(3, 2) << std::endl;
}
