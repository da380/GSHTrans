
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
  int N = 0;
  Float theta = 0.6;
  Wigner d(L, M, N, theta);

  for (int l : d.Degrees()) {
    cout << endl;
    for (int n : d.UpperIndices(l)) {
      cout << n << " " << l << " |  ";
      for (int m : d.Orders(l)) {
        cout << d(l, m, n) << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
}
