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

// Declarations for some utility functions.
int Sign(int m) { return m % 2 ? -1 : 1; }
template <std::floating_point Float>
Float WignerValueMaxOrderAtUpperIndex(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerValueMinOrderAtUpperIndex(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerValueMaxUpperIndexAtOrder(int, int, Float, Float, bool, bool);

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
  std::vector<Float> data;  // vector containing values

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

  template <typename ExecutionPolicy>
  void SetInitialValues(int, Float, Float, bool, bool, ExecutionPolicy);

  template <typename ExecutionPolicy>
  void OneTermRecursion(int, int, Float, std::vector<Float>&,
                        std::vector<Float>&, ExecutionPolicy);

  template <typename ExecutionPolicy>
  void TwoTermRecursion(int, int, Float, std::vector<Float>&,
                        std::vector<Float>&, ExecutionPolicy);
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
  auto sinHalf = std::sin(theta / 2);
  auto cosHalf = std::cos(theta / 2);
  constexpr auto tiny = std::numeric_limits<Float>::min();
  auto logSinHalf = sinHalf > tiny ? std::log(sinHalf) : static_cast<Float>(0);
  auto logCosHalf = cosHalf > tiny ? std::log(cosHalf) : static_cast<Float>(0);
  auto atLeft = sinHalf < tiny;
  auto atRight = cosHalf < tiny;

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

    // Set the values for l == n;
    SetInitialValues(n, logSinHalf, logCosHalf, atLeft, atRight, policy);


    // Set values at the next degree using formulae and one-term recursion.
    if (n + 1 <= L) {
      auto l = n + 1;
      auto mMax = std::min(l, M);

      // Deal with m == -l using formulae if needed.
      if (l <= M) {
        data.push_back(WignerValueMinOrderAtUpperIndex(
            l, n, logSinHalf, logCosHalf, atLeft, atRight));
      }

      // Deal with interior orders using one-term recusion.
      {
        // Set the iterators.
        auto startMinusOne = begin(n, l - 1);
        auto finishMinusOne = end(n, l - 1);
        auto start = begin(n, l);
        // Pre-compute some factors outside the loop.
        auto alpha = (2 * l - 1) * l * cos * sqIntInv[l + n];
        auto beta = (2 * l - 1) * sqIntInv[l + n];
        std::transform(
            policy, startMinusOne, finishMinusOne, start, [&](auto& minusOne) {
              auto m = static_cast<int>(&minusOne - &*startMinusOne) - mMax + 1;
              auto fOne =
                  (alpha - beta * m) * sqIntInv[l - m] * sqIntInv[l + m];
              return fOne * minusOne;
            });
      }

      // Deal with m == l using formulae if needed
      if (l <= M) {
        data.push_back(WignerValueMaxOrderAtUpperIndex(
            l, n, logSinHalf, logCosHalf, atLeft, atRight));
      }
    }

    // Loop over the degrees until the maximum order stops growing.
    for (int l = n + 2; l <= M; l++) {
      // Deal with m == -l using formulae
      data.push_back(WignerValueMinOrderAtUpperIndex(
          l, n, logSinHalf, logCosHalf, atLeft, atRight));

      // Deal with interior orders using two-term recusion.
      {
        // set the iterators.
        auto startMinusTwo = begin(n, l - 2);
        auto finishMinusTwo = end(n, l - 2);
        auto startMinusOne = std::next(begin(n, l - 1));
        auto start = std::next(begin(n, l), 2);
        // Pre-compute some common factors outside the loop.
        auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - n] * sqIntInv[l + n];
        auto beta = (2 * l - 1) * n * sqIntInv[l - n] * sqIntInv[l + n] /
                    static_cast<Float>(l - 1);
        auto gamma = l * sqInt[l - 1 - n] * sqInt[l - 1 + n] * sqIntInv[l - n] *
                     sqIntInv[l + n] / static_cast<Float>(l - 1);

        std::transform(
            policy, startMinusTwo, finishMinusTwo, startMinusOne, start,
            [&](auto& minusTwo, auto& minusOne) {
              auto m = static_cast<int>(&minusTwo - &*startMinusTwo) - l + 2;
              auto denom = sqIntInv[l - m] * sqIntInv[l + m];
              auto fOne = (alpha - beta * m) * denom;
              auto fTwo = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
              return fOne * minusOne - fTwo * minusTwo;
            });
      }

      // Deal with m == l using formulae.
      data.push_back(WignerValueMaxOrderAtUpperIndex(
          l, n, logSinHalf, logCosHalf, atLeft, atRight));
    }  // End of first loop over degrees.

    // If needed, deal with l = M+1;
    if (M + 1 <= L) {
      auto l = M + 1;

      // Deal with m == -M using one-term recursion
      {
        auto m = -M;
        auto fOne = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                    sqIntInv[l - n] * sqIntInv[l + n] * sqIntInv[l - m] *
                    sqIntInv[l + m] / static_cast<Float>(l - 1);
        auto minusOne = *begin(n, l - 1);
        data.push_back(fOne * minusOne);
      }

      // Deal with interior orders using two-term recusion.
      {
        // set the iterators.
        auto startMinusTwo = begin(n, l - 2);
        auto finishMinusTwo = end(n, l - 2);
        auto startMinusOne = std::next(begin(n, l - 1));
        auto start = std::next(begin(n, l));

        // Pre-compute some common factors outside the loop.
        auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - n] * sqIntInv[l + n];
        auto beta = (2 * l - 1) * n * sqIntInv[l - n] * sqIntInv[l + n] /
                    static_cast<Float>(l - 1);
        auto gamma = l * sqInt[l - 1 - n] * sqInt[l - 1 + n] * sqIntInv[l - n] *
                     sqIntInv[l + n] / static_cast<Float>(l - 1);

        std::transform(
            policy, startMinusTwo, finishMinusTwo, startMinusOne, start,
            [&](auto& minusTwo, auto& minusOne) {
              auto m = static_cast<int>(&minusTwo - &*startMinusTwo) - M + 1;
              auto denom = sqIntInv[l - m] * sqIntInv[l + m];
              auto fOne = (alpha - beta * m) * denom;
              auto fTwo = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
              return fOne * minusOne - fTwo * minusTwo;
            });
      }

      // Deal with m == M with one-term recursion.
      {
        auto m = M;
        auto fOne = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                    sqIntInv[l - n] * sqIntInv[l + n] * sqIntInv[l - m] *
                    sqIntInv[l + m] / static_cast<Float>(l - 1);
        auto minusOne = *std::prev(end(n, l - 1));
        data.push_back(fOne * minusOne);
      }
    }

    // If needed, deal with M + 1 < l <= L
    for (int l = M + 2; l <= L; l++) {
      // set the iterators.
      auto startMinusTwo = begin(n, l - 2);
      auto finishMinusTwo = end(n, l - 2);
      auto startMinusOne = begin(n, l - 1);
      auto start = begin(n, l);

      // Pre-compute some common factors outside the loop.
      auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - n] * sqIntInv[l + n];
      auto beta = (2 * l - 1) * n * sqIntInv[l - n] * sqIntInv[l + n] /
                  static_cast<Float>(l - 1);
      auto gamma = l * sqInt[l - 1 - n] * sqInt[l - 1 + n] * sqIntInv[l - n] *
                   sqIntInv[l + n] / static_cast<Float>(l - 1);

      std::transform(
          policy, startMinusTwo, finishMinusTwo, startMinusOne, start,
          [&](auto& minusTwo, auto& minusOne) {
            auto m = static_cast<int>(&minusTwo - &*startMinusTwo) - M;
            auto denom = sqIntInv[l - m] * sqIntInv[l + m];
            auto fOne = (alpha - beta * m) * denom;
            auto fTwo = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
            return fOne * minusOne - fTwo * minusTwo;
          });

    }  // End of second loop over degrees.

  }  // End of loop over upper indices.
}

template <std::floating_point Float>
template <typename ExecutionPolicy>
void WignerValues<Float>::SetInitialValues(int l, Float logSinHalf,
                                           Float logCosHalf, bool atLeft,
                                           bool atRight,
                                           ExecutionPolicy policy) {
  auto mMax = std::min(l, M);
  auto start = begin(l, l);
  auto finish = end(l, l);
  std::transform(policy, start, finish, start, [&](auto& p) {
    auto m = static_cast<int>(&p - &*start) - mMax;
    return WignerValueMaxUpperIndexAtOrder(l, m, logSinHalf, logCosHalf, atLeft,
                                           atRight);
  });
}

template <std::floating_point Float>
template <typename ExecutionPolicy>
void WignerValues<Float>::OneTermRecursion(int n, int l, Float cos,
                                           std::vector<Float>& sqInt,
                                           std::vector<Float>& sqIntInv,
                                           ExecutionPolicy policy) {
  // Work out the offsets for iterators
  int minusOneOffSet = 0;
  int currentOffSet = 0;
  if (l <= M) {
    minusOneOffSet = 1;
    currentOffSet = 2;
  }
  if (l == M + 1) {
    minusOneOffSet = 1;
    currentOffSet = 1;
  }

  // Set the maximum degree
  int mMax = M;
  if (l <= M) {
    mMax = l;
  }
  if (l == M + 1) {
    mMax = M - 1;
  }
}

template <std::floating_point Float>
template <typename ExecutionPolicy>
void WignerValues<Float>::TwoTermRecursion(int n, int l, Float cos,
                                           std::vector<Float>& sqInt,
                                           std::vector<Float>& sqIntInv,
                                           ExecutionPolicy policy) {
  // Work out the offsets for iterators
  int minusOneOffSet = 0;
  int currentOffSet = 0;
  if (l <= M) {
    minusOneOffSet = 1;
    currentOffSet = 2;
  }
  if (l == M + 1) {
    minusOneOffSet = 1;
    currentOffSet = 1;
  }

  // Set the maximum degree
  int mMax = l <= M + 1 ? l - 2 : M;

  // set the iterators.
  auto startMinusTwo = begin(n, l - 2);
  auto finishMinusTwo = end(n, l - 2);
  auto startMinusOne = std::next(begin(n, l - 1), minusOneOffSet);
  auto start = std::next(begin(n, l), currentOffSet);

  // Pre-compute some common factors outside the loop.
  auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - n] * sqIntInv[l + n];
  auto beta = (2 * l - 1) * n * sqIntInv[l - n] * sqIntInv[l + n] /
              static_cast<Float>(l - 1);
  auto gamma = l * sqInt[l - 1 - n] * sqInt[l - 1 + n] * sqIntInv[l - n] *
               sqIntInv[l + n] / static_cast<Float>(l - 1);

  std::transform(
      policy, startMinusTwo, finishMinusTwo, startMinusOne, start,
      [&](auto& minusTwo, auto& minusOne) {
        auto m = static_cast<int>(&minusTwo - &*startMinusTwo) - mMax;
        auto denom = sqIntInv[l - m] * sqIntInv[l + m];
        auto fOne = (alpha - beta * m) * denom;
        auto fTwo = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
        return fOne * minusOne - fTwo * minusTwo;
      });
}

template <std::floating_point Float>
Float WignerValueMaxOrderAtUpperIndex(int l, int n, Float logSinHalf,
                                      Float logCosHalf, bool atLeft,
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
Float WignerValueMinOrderAtUpperIndex(int l, int n, Float logSinHalf,
                                      Float logCosHalf, bool atLeft,
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
Float WignerValueMaxUpperIndexAtOrder(int l, int m, Float logSinHalf,
                                      Float logCosHalf, bool atLeft,
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
