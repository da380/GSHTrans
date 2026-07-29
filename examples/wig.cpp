#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

// ============================================================================
// Wigner3j Class Definition
// ============================================================================

/**
 * @brief Class encapsulating the calculation of Wigner 3-j symbols.
 * * Computes and stores the matrix A such that:
 * A(m1, m2) = (-1)^m1 * (  l1     is     l2 )
 * ( -m1   m1-m2    m2 )
 * * Memory is managed internally in standard C++ row-major order.
 * Values are accessed using physical quantum numbers (m1, m2).
 * * @tparam T A floating-point type (e.g., double, float).
 */
template <std::floating_point T = double>
class Wigner3j {
 private:
  int l1_, is_, l2_;
  int rows_, cols_;
  std::vector<T> data_;

  void compute() {
    std::fill(data_.begin(), data_.end(), static_cast<T>(0.0));

    // Triangle inequality check
    if (is_ < std::abs(l1_ - l2_) || is_ > (l1_ + l2_)) {
      return;  // Leaves matrix as all zeros
    }

    // Standard C++ row-major internal indexing
    auto A = [&](int r, int c) -> T& { return data_[r * cols_ + c]; };

    // 1. Compute the corner value
    auto r1 = static_cast<T>(1.0) /
              std::sqrt(static_cast<T>(2 * std::max(l1_, l2_) + 1));
    auto ldel = std::abs(l1_ - l2_);
    auto num = is_ - ldel;

    for (auto n = 1; n <= num; ++n) {
      auto isc = ldel + n;
      r1 *= std::sqrt(static_cast<T>(l1_ - isc + l2_ + 1) /
                      static_cast<T>(l1_ + l2_ + isc + 1));
    }
    A(0, 0) = r1;

    // 2. Compute first row
    auto num1 = std::min(is_ - l1_ + l2_, l2_);
    for (auto n = 1; n <= num1; ++n) {
      auto m2 = -l2_ + n;
      A(0, n) =
          -A(0, n - 1) *
          std::sqrt(static_cast<T>((l1_ + is_ + m2) * (is_ - l1_ - m2 + 1)) /
                    static_cast<T>((l2_ - m2 + 1) * (l2_ + m2)));
    }

    // 3. Compute first column
    auto num2 = is_ - l2_ + l1_;
    for (auto n = 1; n <= num2; ++n) {
      auto m1 = -l1_ + n;
      A(n, 0) =
          -A(n - 1, 0) *
          std::sqrt(static_cast<T>((l2_ + is_ + m1) * (-l2_ + is_ - m1 + 1)) /
                    static_cast<T>((l1_ + m1) * (l1_ - m1 + 1)));
    }

    // 4. Iterate south-east
    auto iss = is_ * (is_ + 1) - l1_ * (l1_ + 1) - l2_ * (l2_ + 1);
    auto numd = std::min(2 * is_ + 1, is_ + l1_);

    for (auto nd = 1; nd <= numd; ++nd) {
      auto m1b = std::max(-l1_, is_ - l2_ - nd + 1);
      auto m2b = std::max(-l2_, -is_ - l1_ + nd - 1);
      auto numit = std::min(-m2b, l1_ - m1b);

      for (auto nit = 1; nit <= numit; ++nit) {
        auto m1 = m1b + nit;
        auto m2 = m2b + nit;

        auto i = m1 + l1_;
        auto j = m2 + l2_;

        auto r1_val =
            (i >= 2 && j >= 2) ? A(i - 2, j - 2) : static_cast<T>(0.0);

        auto val = -r1_val *
                   std::sqrt(static_cast<T>((l1_ + m1 - 1) * (l1_ - m1 + 2) *
                                            (l2_ - m2 + 2) * (l2_ + m2 - 1)));
        val += static_cast<T>(iss + 2 * (m1 - 1) * (m2 - 1)) * A(i - 1, j - 1);
        val /= std::sqrt(static_cast<T>((l1_ + m1) * (l1_ - m1 + 1) *
                                        (l2_ - m2 + 1) * (l2_ + m2)));

        A(i, j) = val;
      }
    }

    // 5. Fill the rest of the array via symmetry
    auto sgn = ((l1_ + l2_ + is_) % 2 != 0) ? static_cast<T>(-1.0)
                                            : static_cast<T>(1.0);

    if (l2_ != 0) {
      for (auto i = 0; i < rows_; ++i) {
        for (auto m2 = 1; m2 <= l2_; ++m2) {
          A(i, l2_ + m2) = sgn * A(2 * l1_ - i, l2_ - m2);
        }
      }
    }

    // 6. Adjust signs for specific rows
    auto start_i = (l2_ + is_ + 1) % 2;
    for (auto i = start_i; i < rows_; i += 2) {
      for (auto j = 0; j < cols_; ++j) {
        A(i, j) = -A(i, j);
      }
    }
  }

 public:
  Wigner3j(int l1, int is, int l2)
      : l1_(l1),
        is_(is),
        l2_(l2),
        rows_(2 * l1 + 1),
        cols_(2 * l2 + 1),
        data_(rows_ * cols_, static_cast<T>(0.0)) {
    compute();
  }

  T operator()(int m1, int m2) const {
    if (m1 < -l1_ || m1 > l1_ || m2 < -l2_ || m2 > l2_) {
      throw std::out_of_range("Wigner3j: m1 or m2 out of valid bounds.");
    }
    auto i = m1 + l1_;
    auto j = m2 + l2_;
    return data_[i * cols_ + j];
  }

  int l1() const { return l1_; }
  int is() const { return is_; }
  int l2() const { return l2_; }
};

// ============================================================================
// Unit Tests
// ============================================================================

template <typename T>
bool approx_equal(T a, T b, T epsilon = 1e-12) {
  return std::abs(a - b) < epsilon;
}

void test_unit_norm() {
  int l1 = 2, is = 3, l2 = 2;
  Wigner3j<double> w3j(l1, is, l2);

  double sum_sq = 0.0;
  for (int m1 = -l1; m1 <= l1; ++m1) {
    for (int m2 = -l2; m2 <= l2; ++m2) {
      double val = w3j(m1, m2);
      sum_sq += (val * val);
    }
  }

  assert(approx_equal(sum_sq, 1.0) &&
         "Test Failed: Sum of squares must equal 1.0");
  std::cout << "[PASS] Unit Norm (Sum of squares = " << sum_sq << ")\n";
}

void test_selection_rules() {
  int l1 = 3, is = 1, l2 = 3;
  Wigner3j<double> w3j(l1, is, l2);

  for (int m1 = -l1; m1 <= l1; ++m1) {
    for (int m2 = -l2; m2 <= l2; ++m2) {
      if (std::abs(m1 - m2) > is) {
        assert(approx_equal(w3j(m1, m2), 0.0) &&
               "Test Failed: Physics violation not yielding 0.0");
      }
    }
  }
  std::cout << "[PASS] Selection Rules (Impossible states are 0.0)\n";
}

void test_parity_rule() {
  // Sum of degrees: 2 + 3 + 2 = 7 (Odd number)
  int l1 = 2, is = 3, l2 = 2;
  Wigner3j<double> w3j(l1, is, l2);

  assert(approx_equal(w3j(0, 0), 0.0) && "Test Failed: Parity rule violated");
  std::cout << "[PASS] Parity Rule (Odd L sum with m=0 yields 0.0)\n";
}

void test_triangle_inequality() {
  // is=10 is way outside |2-2| to |2+2|
  int l1 = 2, is = 10, l2 = 2;
  Wigner3j<double> w3j(l1, is, l2);

  double sum = 0.0;
  for (int m1 = -l1; m1 <= l1; ++m1) {
    for (int m2 = -l2; m2 <= l2; ++m2) {
      sum += std::abs(w3j(m1, m2));
    }
  }

  assert(approx_equal(sum, 0.0) &&
         "Test Failed: Out of bounds 'is' did not yield a zero matrix");
  std::cout
      << "[PASS] Triangle Inequality (Invalid degrees yield zero matrix)\n";
}

void test_known_value() {
  // For l1=1, is=1, l2=1, m1=1, m2=0
  // The standard physical 3-j symbol for (1, 1, 1 / -1, 1, 0) is actually
  // -1/sqrt(6). Woodhouse formulation applies (-1)^m1: (-1)^1 * (-1/sqrt(6)) =
  // +1/sqrt(6)
  Wigner3j<double> w3j(1, 1, 1);

  // Corrected to positive 1/sqrt(6)
  double expected = 1.0 / std::sqrt(6.0);
  double actual = w3j(1, 0);

  assert(approx_equal(actual, expected) &&
         "Test Failed: Known value does not match");
  std::cout << "[PASS] Known Value Calculation\n";
}

// ============================================================================
// Main Driver
// ============================================================================
int main() {
  std::cout << "Running Wigner 3-j Tests...\n";
  std::cout << std::string(50, '-') << "\n";

  test_unit_norm();
  test_selection_rules();
  test_parity_rule();
  test_triangle_inequality();
  test_known_value();

  std::cout << std::string(50, '-') << "\n";
  std::cout << "ALL TESTS PASSED SUCCESSFULLY.\n";

  return 0;
}