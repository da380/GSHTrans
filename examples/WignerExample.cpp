
#include <cmath>
#include <concepts>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <tuple>

#include "GSHTrans/Core"

int main() {
  using Float = double;

  using std::cout;
  using std::endl;
  cout.setf(std::ios_base::scientific);
  cout.setf(std::ios_base::showpos);
  cout.precision(8);
  using namespace GSHTrans;

  int L = 5;
  int M = L;
  int n = 0;
  Float theta = 0.6;
  Wigner<Float, AllOrders, FullyNormalised> d(L, M, n, theta);

  auto p = MakeWigner<Float, AllOrders, FullyNormalised>(L, M, n, theta);
}
