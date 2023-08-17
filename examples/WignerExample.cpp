

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

  using Float = double;

  int n = 100;
  Float theta1 = 0;
  Float theta2 = std::numbers::pi_v<Float>;
  Float dtheta = (theta2 - theta1) / static_cast<Float>(n - 1);

  int l = 2;
  int N = -2;

  std::ofstream file("WignerExample.out");
  for (int i = 0; i < n; i++) {
    Float theta = theta1 + i * dtheta;
    WignerLN d(l, N, theta, Normalisation::Ortho);
    WignerN d2(l, l, N, theta, Normalisation::Ortho);
    file << theta;
    for (int m = -l; m <= l; m++) {
      file << " " << d(m) - d2(l, m);
    }
    file << std::endl;
  }
}
