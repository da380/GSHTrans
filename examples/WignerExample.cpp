

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
  int n = 5;

  auto size = WingerNStorage<All>(lMax, mMax, n);

  std::vector<double> data(size);

  double th = 1.9;

  auto p = WignerN(data.begin(), data.end(), lMax, mMax, n, th);

  auto p1 = WignerArrayN(lMax, mMax, n, th);

  for (int l = n; l < lMax; l++) {
    auto mS = std::min(mMax, l);
    for (int m = -mS; m <= mS; m++) {
      std::cout << p(l, m) << " " << p1(l, m) << std::endl;
    }
  }
}
