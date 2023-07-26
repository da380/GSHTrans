#ifndef WIGNER_GUARD_H
#define WIGNER_GUARD_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <limits>
#include <numbers>
#include <numeric>
#include <tuple>
#include <vector>

#include "Indexing.h"

namespace GSHT {

// Define the WignerValues class.
template <std::floating_point Float>
class WignerValues {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = typename std::vector<Float>::iterator;
  using difference_type = std::vector<Float>::difference_type;

  // Delete default constructor.
  WignerValues() = delete;

  // Constructor taking in execution policy as argument.
  template <typename ExecutionPolicy>
  WignerValues(int L, int M, int N, Float theta, ExecutionPolicy policy);

  // Constructor using default construction policy.
  WignerValues(int L, int M, int N, Float theta)
      : WignerValues(L, M, N, theta, std::execution::seq) {}

  // Geters for basic data.

  Float Angle() const { return theta; }
  int MaxDegree() const { return L; }
  int MaxOrder() const { return M; }
  int MaxUpperIndex() const { return N; }

  // Functions to return iterators to the data.
  iterator begin() { return data.begin(); }
  iterator end() { return data.end(); }
  iterator begin(int n) { return std::next(begin(), Count(n - 1)); }
  iterator end(int n) { return std::next(begin(), Count(n)); }
  iterator begin(int n, int l) { return std::next(begin(), Count(n, l - 1)); }
  iterator end(int n, int l) { return std::next(begin(), Count(n, l)); }

 private:
  int L;  // Maximum degree
  int M;  // Maximum order
  int N;  // Maximum upper index

  Float theta;              // angle
  std::vector<Float> data;  // vector containing struct

  // Returns the number of elements up to and including upper index n.
  // Note that n = -1 is a valid input with return value 0 which is needed.
  constexpr difference_type Count(int n) const {
    return (n + 1) * ((M + 1) * (M + 1) + (L - M) * (2 * M + 1)) -
           (n * (n + 1) * (2 * n + 1)) / 6;
  }

  // Returns the number of elements up to and including degree l at upper index
  // n. Note that l = -1 is a valid input with this case being needed.
  constexpr difference_type Count(int n, int l) const {
    int i = Count(n - 1);
    if (l >= n) {
      i -= n * n;
      if (l <= M) {
        i += (l + 1) * (l + 1);
      } else {
        i += (M + 1) * (M + 1) + (l - M) * (2 * M + 1);
      }
    }
    return i;
  }

  // Declare/define some utility functions.
  int Sign(int m) { return m % 2 ? -1 : 1; }
  Float WignerValueMaxOrderAtUpperIndex(int, int, Float, Float, bool, bool);
  Float WignerValueMinOrderAtUpperIndex(int, int, Float, Float, bool, bool);
  Float WignerValueMaxUpperIndexAtOrder(int, int, Float, Float, bool, bool);
};

template <std::floating_point Float>
template <typename ExecutionPolicy>
WignerValues<Float>::WignerValues(int L, int M, int N, Float theta,
                                  ExecutionPolicy policy)
    : L{L}, M{M}, N{N}, theta{theta} {
  // Check the maximum degree is non-negative.
  assert(L >= 0);

  // Check the maximum order is in range.
  assert(M >= 0 && M <= L);

  // Check the upper index is in range.
  assert(N >= 0 && N <= L);

  // If N > M swap their values
  if (N > M) std::swap(M, N);

  // Initialise the data vector.
  data = std::vector<Float>(Count(N));

  // Pre-compute and store trigonometric terms
  auto cos = std::cos(theta);
  auto logSinHalf = std::sin(theta / 2);
  auto logCosHalf = std::cos(theta / 2);
  auto atLeft = logSinHalf < std::numeric_limits<Float>::min();
  auto atRight = logCosHalf < std::numeric_limits<Float>::min();
  logSinHalf = atLeft ? static_cast<Float>(0) : std::log(logSinHalf);
  logCosHalf = atRight ? static_cast<Float>(0) : std::log(logCosHalf);

  // Pre-compute and store square roots and their inverses for natural numbers
  // up to L + M
  std::vector<Float> sqInt(L + M + 1);
  std::transform(
      policy, sqInt.begin(), sqInt.end(), sqInt.begin(),
      [&](auto& x) { return std::sqrt(static_cast<Float>(&x - &sqInt[0])); });
  std::vector<Float> sqIntInv(L + M + 1);
  std::transform(
      policy, sqInt.begin(), sqInt.end(), sqIntInv.begin(), [](auto x) {
        return x > static_cast<Float>(0) ? 1 / x : static_cast<Float>(0);
      });

  for (int n = 0; n <= N; n++) {
    // Set the values for l == n
    {
      auto l = n;
      auto mMax = std::min(l, M);
      auto start = begin(l, l);
      auto finish = end(l, l);
      std::transform(policy, start, finish, start, [&](auto& p) {
        auto m = static_cast<int>(&p - &*start) - mMax;
        return WignerValueMaxUpperIndexAtOrder(l, m, logSinHalf, logCosHalf,
                                               atLeft, atRight);
      });
    }

    // Set the values for l == n+1 if needed.
    if (n < L) {
      // Degree and maximum order
      auto l = n + 1;
      auto mMax = std::min(l, M);

      // Set iterators
      auto startMinus1 = begin(n, l - 1);
      auto finishMinus1 = end(n, l - 1);
      auto start = begin(n, l);

      // Add in value at m == -l if needed.
      if (l <= M) {
        *start++ = WignerValueMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                                   atLeft, atRight);
      }

      // Add in interior orders using one-term recursion.
      {
        auto alpha = (2 * l - 1) * l * cos * sqIntInv[l + n];
        auto beta = (2 * l - 1) * sqIntInv[l + n];
        std::transform(
            policy, startMinus1, finishMinus1, start, [&](auto& minus1) {
              auto m = static_cast<int>(&minus1 - &*startMinus1) - mMax + 1;
              auto f1 = (alpha - beta * m) * sqIntInv[l - m] * sqIntInv[l + m];
              return f1 * minus1;
            });
      }

      // Add in value at m == l if needed
      if (l <= M) {
        *std::next(start, 2 * l - 1) = WignerValueMaxOrderAtUpperIndex(
            l, n, logSinHalf, logCosHalf, atLeft, atRight);
      }
    }

    // Now do the remaining degrees.
    for (int l = n + 2; l <= L; l++) {
      // Starting order within two-term recursion
      auto mStart = std::min(l, M);

      // Set iterators
      auto startMinus2 = begin(n, l - 2);
      auto finishMinus2 = end(n, l - 2);
      auto startMinus1 = begin(n, l - 1);
      auto start = begin(n, l);

      // Add in lower boundary terms if still growing.
      if (l <= M) {
        // Add in the m == -l term.
        *start++ = WignerValueMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                                   atLeft, atRight);
        // Now do the m == -l+1 term using one-point recursion.
        {
          auto m = -l + 1;
          auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                    sqIntInv[l - n] * sqIntInv[l + n] * sqIntInv[l - m] *
                    sqIntInv[l + m] / static_cast<Float>(l - 1);
          *start++ = f1 * (*startMinus1++);
        }
        // Update the starting order for two-term recursion.
        mStart -= 2;
      }

      // Add in the lower boundary term at the critical degree
      if (l == M + 1) {
        auto m = -M;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                  sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                  static_cast<Float>(l - 1);
        *start++ = f1 * (*startMinus1++);
        // Update the starting order for two-term recursion.
        mStart -= 1;
      }

      // Apply two-term recusion for the interior orders.
      {
        auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - n] * sqIntInv[l + n];
        auto beta = (2 * l - 1) * n * sqIntInv[l - n] * sqIntInv[l + n] /
                    static_cast<Float>(l - 1);
        auto gamma = l * sqInt[l - 1 - n] * sqInt[l - 1 + n] * sqIntInv[l - n] *
                     sqIntInv[l + n] / static_cast<Float>(l - 1);

        std::transform(
            policy, startMinus2, finishMinus2, startMinus1, start,
            [&](auto& minus2, auto& minus1) {
              auto m = static_cast<int>(&minus2 - &*startMinus2) - mStart;
              auto denom = sqIntInv[l - m] * sqIntInv[l + m];
              auto f1 = (alpha - beta * m) * denom;
              auto f2 = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
              return f1 * minus1 - f2 * minus2;
            });
      }

      // Add in the upper boundary terms if still growing.
      if (l <= M) {
        // Update the iterator
        start = std::next(start, 2 * l - 2);
        // Add in m == -l+1 term using one-point recursion.
        {
          auto m = -l + 1;
          auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                    sqIntInv[l - n] * sqIntInv[l + n] * sqIntInv[l - m] *
                    sqIntInv[l + m] / static_cast<Float>(l - 1);
          *start++ = f1 * (*startMinus1++);
        }
        // Now do m == l.
        *start++ = WignerValueMaxOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                                   atLeft, atRight);
      }

      // Add in the upper boundary term at the crtiical degree.
      if (l == M + 1) {
        // Update the iterators.
        startMinus1 = std::next(startMinus1, M - 2);
        start = std::next(start, M - 2);
	auto m = M;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                  sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                  static_cast<Float>(l - 1);
        *start++ = f1 * (*startMinus1++);
      }
    }
  }
}

template <std::floating_point Float>
Float WignerValues<Float>::WignerValueMaxOrderAtUpperIndex(int l, int n,
                                                           Float logSinHalf,
                                                           Float logCosHalf,
                                                           bool atLeft,
                                                           bool atRight) {
  // Check the inputs.
  assert(l >= 0);
  assert(std::abs(n) <= l);

  // Deal with l == 0 case
  if (l == 0) return static_cast<Float>(1);

  // Deal with special case at the left boundary.
  if (atLeft) {
    return n == l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with special case at the right boundary.
  if (atRight) {
    return n == -l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with the general case.
  auto Fl = static_cast<Float>(l);
  auto Fn = static_cast<Float>(n);
  using std::exp;
  using std::lgamma;
  return Sign(l + n) * exp(0.5 * (lgamma(2 * Fl + 1) - lgamma(Fl + Fn + 1) -
                                  lgamma(Fl - Fn + 1)) +
                           (Fl - Fn) * logSinHalf + (Fl + Fn) * logCosHalf);
}

template <std::floating_point Float>
Float WignerValues<Float>::WignerValueMinOrderAtUpperIndex(int l, int n,
                                                           Float logSinHalf,
                                                           Float logCosHalf,
                                                           bool atLeft,
                                                           bool atRight) {
  // Check the inputs.
  assert(l >= 0);
  assert(std::abs(n) <= l);

  // Deal with l == 0 case
  if (l == 0) return static_cast<Float>(1);

  // Deal with special case at the left boundary.
  if (atLeft) {
    return n == -l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with special case at the right boundary.
  if (atRight) {
    return n == l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with the general case.
  auto Fl = static_cast<Float>(l);
  auto Fn = static_cast<Float>(n);
  using std::exp;
  using std::lgamma;
  return exp(
      0.5 * (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
      (Fl + Fn) * logSinHalf + (Fl - Fn) * logCosHalf);
}

template <std::floating_point Float>
Float WignerValues<Float>::WignerValueMaxUpperIndexAtOrder(int l, int m,
                                                           Float logSinHalf,
                                                           Float logCosHalf,
                                                           bool atLeft,
                                                           bool atRight) {
  // Check the inputs.
  assert(l >= 0);
  assert(std::abs(m) <= l);

  // Deal with l == 0 case
  if (l == 0) return static_cast<Float>(1);

  // Deal with special case at the left boundary.
  if (atLeft) {
    return m == l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with special case at the right boundary.
  if (atRight) {
    return m == -l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with the general case.
  auto Fl = static_cast<Float>(l);
  auto Fm = static_cast<Float>(m);
  using std::exp;
  using std::lgamma;
  return exp(
      0.5 * (lgamma(2 * Fl + 1) - lgamma(Fl - Fm + 1) - lgamma(Fl + Fm + 1)) +
      (Fl - Fm) * logSinHalf + (Fl + Fm) * logCosHalf);
}

}  // namespace GSHT

#endif  // WIGNER_GUARD_H
