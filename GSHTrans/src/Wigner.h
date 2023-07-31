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

namespace GSHTrans {

// Declare some utility functions
constexpr int Sign(int);
template <std::floating_point Float>
Float WignerMinOrderAtUpperIndex(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerMaxOrderAtUpperIndex(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerMaxUpperIndexAtOrder(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerMinUpperIndexAtOrder(int, int, Float, Float, bool, bool);

// Tags for specifying Wigner class options.
struct AllOrders {};
struct NonNegativeOrders {};
struct FourPiNormalised {};
struct FullyNormalised {};

template <typename Range>
concept OrderRange =
    std::same_as<Range, AllOrders> or std::same_as<Range, NonNegativeOrders>;

template <typename Norm>
concept Normalisation =
    std::same_as<Norm, FourPiNormalised> or std::same_as<Norm, FullyNormalised>;

template <std::floating_point Float, OrderRange Range, Normalisation Norm>
class Wigner {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = typename std::vector<Float>::iterator;
  using const_iterator = typename std::vector<Float>::const_iterator;
  using difference_type = std::vector<Float>::difference_type;

  // Constructor taking in execution policy as argument.
  template <typename ExecutionPolicy>
  Wigner(int L, int M, int n, Float theta, ExecutionPolicy policy);

  Wigner(int L, int M, int n, Float theta)
      : Wigner(L, M, n, theta, std::execution::seq) {}

  // Geters for basic data.
  Float Angle() const { return theta; }
  int MaxDegree() const { return L; }
  int MaxOrder() const { return M; }
  int UpperIndex() const { return n; }

  // Iterators that point to the start of the data.
  iterator begin() { return data.begin(); }
  const_iterator cbegin() const { return data.cbegin(); }

  // Iterators that point to the end of the data.
  iterator end() { return data.end(); }
  const_iterator cend() const { return data.cend(); }

  // Iterators that point to the start of degree l.
  iterator begin(int l) { return std::next(begin(), Count(l - 1)); }
  const_iterator cbegin(int l) const {
    return std::next(cbegin(), Count(l - 1));
  }

  // Iterators that point to the end of degree l.
  iterator end(int l) { return std::next(begin(), Count(l)); }
  const_iterator cend(int l) const { return std::next(cbegin(), Count(l)); }

  // Returns value for given degree and order when all orders are stored.
  Float operator()(int l, int m) const requires std::same_as<Range, AllOrders> {
    assert(0 <= l && l <= L);
    assert(std::abs(m) <= std::min(l, M));
    return *std::next(cbegin(l), l + m);
  }

  // Returns value for given degree and non-negative order when only
  // non-negative orders are stored.
  Float operator()(
      int l, int m) const requires std::same_as<Range, NonNegativeOrders> {
    assert(0 <= l && l <= L);
    assert(0 <= m && m <= std::min(l, M));
    return *std::next(cbegin(l), m);
  }

 private:
  int L;  // Maximum degree
  int M;  // Maximum order
  int n;  // Maximum upper index

  Float theta;              // angle
  std::vector<Float> data;  // vector containing struct

  // Returns the number of values at a given degree when
  // all orders are stored.
  constexpr difference_type Count(
      int l) const requires std::same_as<Range, AllOrders> {
    auto nabs = std::abs(n);
    if (l < nabs) return 0;
    if (l <= M) return (l + 1) * (l + 1) - nabs * nabs;
    return (M + 1) * (M + 1) - nabs * nabs + (l - M) * (2 * M + 1);
  }

  // Returns the number of values at a given degree when
  // only non-negative orders are stored.
  constexpr difference_type Count(
      int l) const requires std::same_as<Range, NonNegativeOrders> {
    auto nabs = std::abs(n);
    if (l < nabs) return 0;
    if (l <= M) return ((l + 1) * (l + 2)) / 2 - (nabs * (nabs + 1)) / 2;
    return ((M + 1) * (M + 2)) / 2 - (nabs * (nabs + 1)) / 2 +
           (l - M) * (M + 1);
  }

  // Returns total number of values.
  constexpr difference_type const Count() { return Count(L); }
};

template <std::floating_point Float, OrderRange Range, Normalisation Norm>
template <typename ExecutionPolicy>
Wigner<Float, Range, Norm>::Wigner(int L, int M, int n, Float theta,
                                   ExecutionPolicy policy)
    : L{L}, M{M}, n{n}, theta{theta} {
  // Check the maximum degree is non-negative.
  assert(L >= 0);

  // Check the maximum order is in range.
  assert(M >= 0 && M <= L);

  // Check the upper index is in range.
  assert(abs(n) <= M);

  // Initialise the data vector.
  data = std::vector<Float>(Count());

  // Pre-compute and store trigonometric terms.
  auto cos = std::cos(theta);
  auto logSinHalf = std::sin(theta / 2);
  auto logCosHalf = std::cos(theta / 2);
  auto atLeft = logSinHalf < std::numeric_limits<Float>::min();
  auto atRight = logCosHalf < std::numeric_limits<Float>::min();
  logSinHalf = atLeft ? static_cast<Float>(0) : std::log(logSinHalf);
  logCosHalf = atRight ? static_cast<Float>(0) : std::log(logCosHalf);

  // Pre-compute and store square roots and their inverses up to L + M.
  std::vector<Float> sqInt(L + M + 1);
  std::transform(
      policy, sqInt.begin(), sqInt.end(), sqInt.begin(),
      [&](auto &x) { return std::sqrt(static_cast<Float>(&x - &sqInt[0])); });
  std::vector<Float> sqIntInv(L + M + 1);
  std::transform(
      policy, sqInt.begin(), sqInt.end(), sqIntInv.begin(), [](auto x) {
        return x > static_cast<Float>(0) ? 1 / x : static_cast<Float>(0);
      });

  // Set the values for l == |n|
  {
    auto l = std::abs(n);
    int mStart = 0;
    if constexpr (std::same_as<Range, AllOrders>) {
      mStart = -std::min(l, M);
    }
    auto start = begin(l);
    auto finish = end(l);
    if (n >= 0) {
      std::transform(policy, start, finish, start, [&](auto &p) {
        auto m = mStart + static_cast<int>(&p - &*start);
        return WignerMaxUpperIndexAtOrder(l, m, logSinHalf, logCosHalf, atLeft,
                                          atRight);
      });
    } else {
      std::transform(policy, start, finish, start, [&](auto &p) {
        auto m = mStart + static_cast<int>(&p - &*start);
        return WignerMinUpperIndexAtOrder(l, m, logSinHalf, logCosHalf, atLeft,
                                          atRight);
      });
    }
  }

  // Set the values for l == n+1 if needed.
  if (std::abs(n) < L) {
    auto l = std::abs(n) + 1;
    auto mStart = -std::min(l, M);
    if constexpr (std::same_as<Range, NonNegativeOrders>) {
      mStart = 0;
    }

    // Set iterators
    auto startMinus1 = begin(l - 1);
    auto finishMinus1 = end(l - 1);
    auto start = begin(l);

    // Add in value at m == -l if needed.
    if constexpr (std::same_as<Range, AllOrders>) {
      if (l <= M) {
        *start++ = WignerMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                              atLeft, atRight);
        // Update the starting order for recursion
        mStart += 1;
      }
    }

    // Add in interior orders using one-term recursion.
    {
      auto alpha = (2 * l - 1) * l * cos * sqIntInv[l + n];
      auto beta = (2 * l - 1) * sqIntInv[l + n];
      std::transform(
          policy, startMinus1, finishMinus1, start, [&](auto &minus1) {
            auto m = mStart + static_cast<int>(&minus1 - &*startMinus1);
            auto f1 = (alpha - beta * m) * sqIntInv[l - m] * sqIntInv[l + m];
            return f1 * minus1;
          });
    }

    // Add in value at m == l if needed
    if (l <= M) {
      auto mStep = 2 * l - 1;
      if constexpr (std::same_as<Range, NonNegativeOrders>) {
        mStep = l;
      }
      *std::next(start, mStep) = WignerMaxOrderAtUpperIndex(
          l, n, logSinHalf, logCosHalf, atLeft, atRight);
    }
  }

  // Now do the remaining degrees.
  for (int l = std::abs(n) + 2; l <= L; l++) {
    // Starting order within two-term recursion
    auto mStart = -std::min(l, M);
    if constexpr (std::same_as<Range, NonNegativeOrders>) {
      mStart = 0;
    }

    // Set iterators
    auto startMinus2 = begin(l - 2);
    auto finishMinus2 = end(l - 2);
    auto startMinus1 = begin(l - 1);
    auto start = begin(l);

    // Add in lower boundary terms if still growing.
    if constexpr (std::same_as<Range, AllOrders>) {
      if (l <= M) {
        // Add in the m == -l term.
        *start++ = WignerMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
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
        mStart += 2;
      }

      // Add in the lower boundary term at the critical degree
      if (l == M + 1) {
        auto m = -M;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                  sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                  static_cast<Float>(l - 1);
        *start++ = f1 * (*startMinus1++);
        // Update the starting order for two-term recursion.
        mStart += 1;
      }
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
          [&](auto &minus2, auto &minus1) {
            auto m = mStart + static_cast<int>(&minus2 - &*startMinus2);
            auto denom = sqIntInv[l - m] * sqIntInv[l + m];
            auto f1 = (alpha - beta * m) * denom;
            auto f2 = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
            return f1 * minus1 - f2 * minus2;
          });
    }

    // Add in the upper boundary terms if still growing.
    if (l <= M) {
      // Update the iterator
      auto mStep = 2 * l - 3;
      if constexpr (std::same_as<Range, NonNegativeOrders>) {
        mStep = l - 1;
      }
      startMinus1 = std::next(startMinus1, mStep);
      start = std::next(start, mStep);
      // Add in m == -l+1 term using one-point recursion.
      {
        auto m = l - 1;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                  sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                  static_cast<Float>(l - 1);
        *start++ = f1 * (*startMinus1++);
      }
      // Now do m == l.
      *start++ = WignerMaxOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                            atLeft, atRight);
    }

    // Add in the upper boundary term at the crtiical degree.
    if (l == M + 1) {
      // Update the iterators.
      auto mStep = 2 * M - 1;
      if constexpr (std::same_as<Range, NonNegativeOrders>) {
        mStep = M;
      }
      startMinus1 = std::next(startMinus1, mStep);
      start = std::next(start, mStep);
      auto m = M;
      auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                static_cast<Float>(l - 1);
      *start++ = f1 * (*startMinus1++);
    }
  }

  if constexpr (std::same_as<Norm, FullyNormalised>) {
    for (int l = 0; l <= L; l++) {
      auto start = begin(l);
      auto finish = end(l);
      std::transform(start, finish, start, [l](auto p) {
        return 0.5 * std::sqrt(static_cast<Float>(2 * l + 1)) *
               std::numbers::inv_sqrtpi_v<Float> * p;
      });
    }
  }
}

constexpr int Sign(int m) { return m % 2 ? -1 : 0; }

template <std::floating_point Float>
Float WignerMinOrderAtUpperIndex(int l, int n, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
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
      static_cast<Float>(0.5) *
          (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
      (Fl + Fn) * logSinHalf + (Fl - Fn) * logCosHalf);
}

template <std::floating_point Float>
Float WignerMaxOrderAtUpperIndex(int l, int n, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
  return Sign(n + l) * WignerMinOrderAtUpperIndex(l, -n, logSinHalf, logCosHalf,
                                                  atLeft, atRight);
}

template <std::floating_point Float>
Float WignerMinUpperIndexAtOrder(int l, int m, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
  return WignerMaxOrderAtUpperIndex(l, -m, logSinHalf, logCosHalf, atLeft,
                                    atRight);
}

template <std::floating_point Float>
Float WignerMaxUpperIndexAtOrder(int l, int m, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
  return WignerMinOrderAtUpperIndex(l, -m, logSinHalf, logCosHalf, atLeft,
                                    atRight);
}

}  // namespace GSHTrans

#endif  // WIGNER_GUARD_H
