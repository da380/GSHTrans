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

  // Constructors.
  WignerValues() = delete;

  template <typename ExecutionPolicy>
  WignerValues(ExecutionPolicy policy, int L, int M, int N, Float theta);

  WignerValues(int L, int M, int N, Float theta)
      : WignerValues(std::execution::seq, L, M, N, theta) {}

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
};

template <std::floating_point Float>
template <typename ExecutionPolicy>
WignerValues<Float>::WignerValues(ExecutionPolicy policy, int L, int M, int N,
                                  Float theta)
    : L{L}, M{M}, N{N}, theta{theta} {
  // Check the maximum degree is non-negative.
  assert(L >= 0);

  // Check the maximum order is in range.
  assert(M >= 0 && M <= L);

  // Check the upper index is in range.
  assert(N >= 0 && N <= L);

  // If N > M swap their values
  if (N > M) std::swap(M, N);

  // Reserve space for the values.
  data.reserve(Count(N));

  // Pre-compute and store trigonometric terms
  auto cos = std::cos(theta);
  auto sinHalf = std::sin(theta / 2);
  auto cosHalf = std::cos(theta / 2);
  constexpr auto tiny = std::numeric_limits<Float>::min();
  auto logSinHalf = sinHalf > tiny ? std::log(sinHalf) : static_cast<Float>(0);
  auto logCosHalf = cosHalf > tiny ? std::log(cosHalf) : static_cast<Float>(0);
  auto atLeft = sinHalf < tiny;
  auto atRight = cosHalf < tiny;

  // Pre-compute and store square roots of integers up to L + M
  std::vector<Float> sqrtInts;
  sqrtInts.reserve(L + M);
  std::generate_n(back_inserter(sqrtInts), L + M + 1, [l = 0]() mutable {
    return std::sqrt(static_cast<Float>(l++));
  });




  // Loop over the non-negative uppper indices.
  for (int n = 0; n <= N; n++) {
    // Set values for degree n
    if (n == 0) {
      // Deal with l == 0 separately
      data.push_back(static_cast<Float>(1));
    } else {
      auto l = n;

      /*
      std::generate_n(policy, begin(n, l), 2 * l + 1, [&, m = -l]() {
        return WignerValueMaxUpperIndexAtOrder(l, m, logSinHalf, logCosHalf,
                                               atLeft, atRight);
      });
      */
      
    }

    // Set values for degree n+1.
    {
      auto l = n + 1;
      // Add term at minimum order directly.
      data.push_back(WignerValueMaxUpperIndexAtOrder(
          l, -l, logSinHalf, logCosHalf, atLeft, atRight));
      // Apply one term recursion for |m| < l.
      auto preFac12 = (2 * l - 1) / sqrtInts[l + n];
      auto preFac11 = preFac12 * l * cos;
      std::transform(policy, begin(n, l - 1), end(n, l - 1), begin(n, l),
                     [=, &sqrtInts, m = -l + 1](auto minus1) mutable {
                       auto fac1 = (preFac11 - m * preFac12) /
                                   (sqrtInts[l - m] * sqrtInts[l + m]);
                       m++;
                       return fac1 * minus1;
                     });

      // Add in term at maximum order directly.
      data.push_back(WignerValueMaxUpperIndexAtOrder(
          l, l, logSinHalf, logCosHalf, atLeft, atRight));
    }

  }  // End loop over upper index

  /*


  // Loop over the upper indices.
  for (int n = 0; n <= N; n++) {
    // Compute the values at l = n for appropriate m directly.
    {
      int l = n;
      for (int m : Orders(l)) {
        data.push_back(ValueMaxUpperIndexAtOrder(l, m));
      }
    }

    // Do the next degree using one-term recusion.
    if (n < L) {
      int l = n + 1;

      // Add in the new term at m = -l
      data.push_back(ValueMinOrderAtUpperIndex(l,n));

      // Add in the new terms for |m| < l using one-term recusion
      auto minusOne = beginForDegreeAtUpperIndex(l - 1, n);
      for (int m : Orders(l-1, M)) {
        Float f1 = Numerator1(l,m,n)/Denominator(l,m,n);
        data.push_back(f1 * (*minusOne++));

      }

      // Add in the new term at m = l
      data.push_back(ValueMaxOrderAtUpperIndex(l,n));
    }

    // Now do the rest using two-term recusion.
    for (int l = n + 2; l <= L; l++) {

      // Add in the new term at m = -l
      data.push_back(ValueMinOrderAtUpperIndex(l,n));

      // Add in the new terms for |m| < l using two-term recusion
      auto minusTwo = beginForDegreeAtUpperIndex(l - 2, n);
      auto minusOne = beginForDegreeAtUpperIndex(l - 1, n);
      for (int m : Orders(l-1, M)) {
        Float d = Denominator(l,m,n);
        Float f1 = Numerator1(l,m,n)/d;
        Float f2 = Numerator2(l,m,n)/d;
        data.push_back( f1 * (*minusOne++) + f2 * (*minusTwo++));
      }

      // Add in the new term at m = l
      data.push_back(ValueMaxOrderAtUpperIndex(l,n));

    }
  }

  */
}

template <std::floating_point Float>
auto WignerValueMaxOrderAtUpperIndex(int l, int n, Float logSinHalf,
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
