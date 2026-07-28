#ifndef WIGNER3J_HPP
#define WIGNER3J_HPP

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <vector>

/**
 * @file wigner3j.hpp
 * @brief Full-grid Wigner 3-j symbols in the Woodhouse convention.
 */

/**
 * @brief Computes and stores every Wigner 3-j symbol on the (m1, m2) grid
 *        for a fixed triple (l1, is, l2).
 *
 * The constructor fills the matrix
 *
 *     A(m1, m2) = (-1)^m1 * / l1    is     l2 \
 *                           \ -m1  m1-m2   m2 /
 *
 * for all -l1 <= m1 <= l1 and -l2 <= m2 <= l2, where the bracketed array is
 * the ordinary Wigner 3-j symbol. The three projections sum to zero by
 * construction: (-m1) + (m1 - m2) + m2 = 0. The extra (-1)^m1 phase and the
 * (-m1, m1-m2, m2) layout reproduce J. H. Woodhouse's original `wig2` routine.
 *
 * Storage is a single row-major buffer; callers index by physical quantum
 * numbers via operator() (unchecked) or at() (bounds-checked). Products of
 * integer factors are formed in floating point, so no intermediate integer
 * overflow occurs for large l (the original F77 default-integer arithmetic
 * overflowed around l ~ 213).
 *
 * @note Accuracy caveat. This is a one-directional recurrence. It is accurate
 *       to machine precision for modest l, but for strongly "stretched"
 *       symbols (is comparable to l1 + l2 with large l) it loses accuracy well
 *       before overflow would occur. As a cheap self-diagnostic the
 *       completeness relation sum_{m1,m2} A(m1,m2)^2 == 1 holds for any valid
 *       triple; a large deviation flags loss of precision. See the notes that
 *       accompany this file for genuinely stable algorithms.
 *
 * @tparam T A floating-point type (default double).
 */
template <std::floating_point T = double>
class Wigner3j {
 public:
  Wigner3j(int l1, int is, int l2)
      : l1_(validated(l1, is, l2)),  // validates before any allocation below
        is_(is),
        l2_(l2),
        rows_(2 * l1 + 1),
        cols_(2 * l2 + 1),
        data_(static_cast<std::size_t>(rows_) * static_cast<std::size_t>(cols_),
              T{0}) {
    compute();
  }

  /// Unchecked access. Precondition: |m1| <= l1 and |m2| <= l2.
  [[nodiscard]] T operator()(int m1, int m2) const noexcept {
    return data_[index(m1 + l1_, m2 + l2_)];
  }

  /// Bounds-checked access; throws std::out_of_range if (m1, m2) is invalid.
  [[nodiscard]] T at(int m1, int m2) const {
    if (m1 < -l1_ || m1 > l1_ || m2 < -l2_ || m2 > l2_) {
      throw std::out_of_range("Wigner3j: (m1, m2) out of valid bounds.");
    }
    return (*this)(m1, m2);
  }

  [[nodiscard]] int l1() const noexcept { return l1_; }
  [[nodiscard]] int is() const noexcept { return is_; }
  [[nodiscard]] int l2() const noexcept { return l2_; }
  [[nodiscard]] int rows() const noexcept { return rows_; }
  [[nodiscard]] int cols() const noexcept { return cols_; }

  /// True iff (l1, is, l2) satisfies the triangle inequality (otherwise the
  /// whole grid is identically zero).
  [[nodiscard]] bool satisfies_triangle() const noexcept {
    return is_ >= std::abs(l1_ - l2_) && is_ <= l1_ + l2_;
  }

 private:
  int l1_, is_, l2_;
  int rows_, cols_;
  std::vector<T> data_;

  // Runs first in the initializer list so bad arguments are rejected before
  // any size is computed from them. Returns l1 unchanged on success.
  static int validated(int l1, int is, int l2) {
    if (l1 < 0 || is < 0 || l2 < 0) {
      throw std::invalid_argument(
          "Wigner3j: l1, is, l2 must all be non-negative.");
    }
    return l1;
  }

  [[nodiscard]] std::size_t index(int r, int c) const noexcept {
    return static_cast<std::size_t>(r) * static_cast<std::size_t>(cols_) +
           static_cast<std::size_t>(c);
  }

  /// Promote an integer to T so that products form in floating point rather
  /// than in (potentially overflowing) int arithmetic.
  static constexpr T f(long long n) noexcept { return static_cast<T>(n); }

  void compute() {
    // Triangle inequality: leave the all-zero grid untouched.
    if (!satisfies_triangle()) {
      return;
    }

    auto A = [&](int r, int c) -> T& { return data_[index(r, c)]; };

    // 1. Corner value A(0, 0).
    T r1 = T{1} / std::sqrt(f(2 * std::max(l1_, l2_) + 1));
    const int ldel = std::abs(l1_ - l2_);
    for (int n = 1; n <= is_ - ldel; ++n) {
      const int isc = ldel + n;
      r1 *= std::sqrt(f(l1_ - isc + l2_ + 1) / f(l1_ + l2_ + isc + 1));
    }
    A(0, 0) = r1;

    // 2. First row.
    const int num1 = std::min(is_ - l1_ + l2_, l2_);
    for (int n = 1; n <= num1; ++n) {
      const int m2 = -l2_ + n;
      A(0, n) =
          -A(0, n - 1) * std::sqrt(f(l1_ + is_ + m2) * f(is_ - l1_ - m2 + 1) /
                                   (f(l2_ - m2 + 1) * f(l2_ + m2)));
    }

    // 3. First column.
    const int num2 = is_ - l2_ + l1_;
    for (int n = 1; n <= num2; ++n) {
      const int m1 = -l1_ + n;
      A(n, 0) =
          -A(n - 1, 0) * std::sqrt(f(l2_ + is_ + m1) * f(-l2_ + is_ - m1 + 1) /
                                   (f(l1_ + m1) * f(l1_ - m1 + 1)));
    }

    // 4. Iterate south-east.
    const long long iss =
        1LL * is_ * (is_ + 1) - 1LL * l1_ * (l1_ + 1) - 1LL * l2_ * (l2_ + 1);
    const int numd = std::min(2 * is_ + 1, is_ + l1_);
    for (int nd = 1; nd <= numd; ++nd) {
      const int m1b = std::max(-l1_, is_ - l2_ - nd + 1);
      const int m2b = std::max(-l2_, -is_ - l1_ + nd - 1);
      const int numit = std::min(-m2b, l1_ - m1b);
      for (int nit = 1; nit <= numit; ++nit) {
        const int m1 = m1b + nit;
        const int m2 = m2b + nit;
        const int i = m1 + l1_;
        const int j = m2 + l2_;

        const T diag = (i >= 2 && j >= 2) ? A(i - 2, j - 2) : T{0};

        T val = -diag * std::sqrt(f(l1_ + m1 - 1) * f(l1_ - m1 + 2) *
                                  f(l2_ - m2 + 2) * f(l2_ + m2 - 1));
        val +=
            static_cast<T>(iss + 2LL * (m1 - 1) * (m2 - 1)) * A(i - 1, j - 1);
        val /= std::sqrt(f(l1_ + m1) * f(l1_ - m1 + 1) * f(l2_ - m2 + 1) *
                         f(l2_ + m2));
        A(i, j) = val;
      }
    }

    // 5. Fill the right half from the left half by symmetry.
    const T sgn = ((l1_ + l2_ + is_) % 2 != 0) ? T{-1} : T{1};
    if (l2_ != 0) {
      for (int i = 0; i < rows_; ++i) {
        for (int m2 = 1; m2 <= l2_; ++m2) {
          A(i, l2_ + m2) = sgn * A(2 * l1_ - i, l2_ - m2);
        }
      }
    }

    // 6. Flip signs on alternating rows.
    const int start_i = (l2_ + is_ + 1) % 2;
    for (int i = start_i; i < rows_; i += 2) {
      for (int j = 0; j < cols_; ++j) {
        A(i, j) = -A(i, j);
      }
    }
  }
};

#endif  // WIGNER3J_HPP
