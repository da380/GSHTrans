

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

  int lMax = 5;
  int mMax = 3;
  int nMax = 2;
  auto theta = std::vector<Real>(2, 1.0);

  auto d = WignerNew<Real, All, NonNegative, Ortho>(lMax, mMax, nMax, theta);

  auto i = 0;
  for (auto n : d.UpperIndices()) {
    for (auto iTheta : d.AngleIndices()) {
      auto e = d(n, iTheta);
      for (auto l : e.Degrees()) {
        auto f = e(l);
        for (auto m : f.Orders()) {
          std::cout << n << " " << iTheta << " " << l << " " << m << " "
                    << f(m) - i++ << std::endl;
        }
      }
    }
    std::cout << "-------------------------------\n";
  }

  /*

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

  */
}
