#include <GSHTrans/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
#include <chrono>
#include <cmath>
#include <concepts>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

int main() {
  using namespace GSHTrans;
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Complex;
  using MRange = All;
  using NRange = All;

  using Grid = GaussLegendreGrid<Real, MRange, NRange>;
}
