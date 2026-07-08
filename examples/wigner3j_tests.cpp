#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "wigner3j.hpp"

// ============================================================================
// Test helpers
// ============================================================================

template <typename T>
bool approx_equal(T a, T b,
                  T tol = T{100} * std::numeric_limits<T>::epsilon()) {
  const T scale = std::max<T>(T{1}, std::max(std::abs(a), std::abs(b)));
  return std::abs(a - b) <= tol * scale;
}

// Sum of squares of the whole grid. For any triple satisfying the triangle
// inequality this equals 1 (completeness); otherwise it is 0.
template <typename T>
T grid_norm_sq(const Wigner3j<T>& w) {
  T s{0};
  for (int m1 = -w.l1(); m1 <= w.l1(); ++m1) {
    for (int m2 = -w.l2(); m2 <= w.l2(); ++m2) {
      const T v = w(m1, m2);
      s += v * v;
    }
  }
  return s;
}

// ============================================================================
// Tests
// ============================================================================

void test_completeness() {
  // Valid triples, including asymmetric ones, should each sum to 1.
  const int triples[][3] = {{2, 3, 2}, {1, 1, 1}, {3, 1, 3}, {2, 2, 2},
                            {4, 3, 5}, {3, 4, 3}, {0, 2, 2}, {5, 5, 5}};
  for (const auto& t : triples) {
    Wigner3j<double> w(t[0], t[1], t[2]);
    const double s = grid_norm_sq(w);
    if (!approx_equal(s, 1.0, 1e-10)) {
      throw std::runtime_error("Completeness failed for (" +
                               std::to_string(t[0]) + "," +
                               std::to_string(t[1]) + "," +
                               std::to_string(t[2]) +
                               "): sum_sq = " + std::to_string(s));
    }
  }
  std::cout << "[PASS] Completeness (sum of squares == 1 for valid triples)\n";
}

void test_selection_rules() {
  const int l1 = 3, is = 1, l2 = 3;
  Wigner3j<double> w(l1, is, l2);
  for (int m1 = -l1; m1 <= l1; ++m1) {
    for (int m2 = -l2; m2 <= l2; ++m2) {
      if (std::abs(m1 - m2) > is && !approx_equal(w(m1, m2), 0.0)) {
        throw std::runtime_error("Selection rule violated.");
      }
    }
  }
  std::cout << "[PASS] Selection Rules (forbidden states are 0)\n";
}

void test_parity_rule() {
  // l1 + is + l2 = 7 (odd) forces the m1 = m2 = 0 entry to vanish.
  Wigner3j<double> w(2, 3, 2);
  if (!approx_equal(w(0, 0), 0.0)) {
    throw std::runtime_error("Parity rule violated.");
  }
  std::cout << "[PASS] Parity Rule (odd l-sum gives A(0,0) = 0)\n";
}

void test_triangle_inequality() {
  // is = 10 lies outside |2-2| .. |2+2|: the whole grid must be zero.
  Wigner3j<double> w(2, 10, 2);
  if (w.satisfies_triangle()) {
    throw std::runtime_error("satisfies_triangle() should be false.");
  }
  if (!approx_equal(grid_norm_sq(w), 0.0)) {
    throw std::runtime_error("Invalid triple did not give a zero grid.");
  }
  std::cout << "[PASS] Triangle Inequality (invalid triple gives zero grid)\n";
}

void test_known_value() {
  // (1 1 1; -1 1 0) = -1/sqrt(6); the (-1)^m1 phase makes A(1,0) = +1/sqrt(6).
  Wigner3j<double> w(1, 1, 1);
  if (!approx_equal(w(1, 0), 1.0 / std::sqrt(6.0))) {
    throw std::runtime_error("Known value mismatch.");
  }
  std::cout << "[PASS] Known Value Calculation\n";
}

void test_bounds_checking() {
  Wigner3j<double> w(2, 3, 2);
  bool threw = false;
  try {
    (void)w.at(3, 0);  // m1 = 3 > l1 = 2
  } catch (const std::out_of_range&) {
    threw = true;
  }
  if (!threw) {
    throw std::runtime_error("at() failed to throw on out-of-range access.");
  }

  threw = false;
  try {
    Wigner3j<double> bad(-1, 0, 0);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  if (!threw) {
    throw std::runtime_error("Constructor failed to reject negative l.");
  }
  std::cout << "[PASS] Bounds & Argument Checking\n";
}

// ============================================================================
// Driver
// ============================================================================

int main() {
  std::cout << "Running Wigner 3-j tests...\n";
  std::cout << std::string(52, '-') << "\n";
  try {
    test_completeness();
    test_selection_rules();
    test_parity_rule();
    test_triangle_inequality();
    test_known_value();
    test_bounds_checking();
  } catch (const std::exception& e) {
    std::cerr << "[FAIL] " << e.what() << "\n";
    return 1;
  }
  std::cout << std::string(52, '-') << "\n";
  std::cout << "ALL TESTS PASSED.\n";
  return 0;
}
