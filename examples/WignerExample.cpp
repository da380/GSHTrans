

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

  auto p1 = Wigner<double, All>(lMax, mMax, n, theta);
  auto p2 = Wigner<double, NonNegative>(lMax, mMax, n, theta);

  auto p3 = Wigner<double, All>(lMax, mMax, n,
                                std::vector<double>{theta, theta, theta});

  auto p4 = AssociatedLegendre<double, All>(lMax, mMax, theta);

  auto p5 = Legendre<double>(lMax, theta);

  std::cout << p1(4, 2) << " " << p2(4, 2) << " " << p3(2, 4, 2) << " "
            << p4(4, 2) << std::endl;

  std::cout << p1(0, 4, 0) << " " << p5(0, 4) << std::endl;
}
