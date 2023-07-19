#ifndef WIGNER_GUARD_H
#define WIGNER_GUARD_H

#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>
#include <tuple>
#include <vector>

#include "Indexing.h"

namespace GSHT {

template <std::floating_point Float>
class WignerValues {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = typename std::vector<Float>::iterator;

  // Constructors.
  WignerValues() = delete;
  WignerValues(int L, int M, int N, Float theta);

  // Geters.
  Float Angle() const { return theta; }
  int MaxDegree() const { return L; }
  int MaxOrder() const { return M; }
  int MaxUpperIndex() const { return N; }

  // Functions to return iterators to the data.
  iterator begin() { return data.begin(); }
  iterator end() { return data.end(); }
  iterator beginForUpperIndex(const int n) {
    return std::next(begin(), ElementsIncludingUpperIndex(n - 1));
  }
  iterator endForUpperIndex(const int n) {
    return std::next(begin(), ElementsIncludingUpperIndex(n));
  }
  iterator beginForDegreeAtUpperIndex(const int l, const int n) {
    return std::next(begin(), ElementsIncludingDegreeAtUpperIndex(l - 1, n));
  }
  iterator endForDegreeAtUpperIndex(const int l, const int n) {
    return std::next(begin(), ElementsIncludingDegreeAtUpperIndex(l, n));
  }

 private:
  int L;  // Maximum degree
  int M;  // Maximum order
  int N;  // Maximum upper index

  Float theta;           // angle
  bool atLeft = false;   // true if theta = 0
  bool atRight = false;  // true if theta = pi
  Float cos;             // cos(theta)
  Float logSinHalf;      // log of sin(theta/2)
  Float logCosHalf;      // log of cos(theta/2)

  std::vector<Float> data;

  constexpr int ElementsIncludingUpperIndex(int n) const {
    return (n + 1) * ((M + 1) * (M + 1) + (L - M) * (2 * M + 1)) -
           (n * (n + 1) * (2 * n + 1)) / 6;
  }

  constexpr int ElementsIncludingDegreeAtUpperIndex(int l, int n) const {
    int i = ElementsIncludingUpperIndex(n - 1);
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

  constexpr int ElementsBeforeDegreeAtUpperIndex(int l, int n) const {
    return ElementsIncludingDegreeAtUpperIndex(l - 1, n);
  }

  Float Sign(int m) { return m % 2 ? -1.0 : 1.0; }
  Float Numerator1(int, int, int);
  Float Numerator2(int, int, int);
  Float Denominator(int, int, int);
  auto RecursionCoefficients(int, int, int);

  Float ValueMaxOrderAtUpperIndex(int, int);
  Float ValueMinOrderAtUpperIndex(int, int);
  Float ValueMaxUpperIndexAtOrder(int, int);


  
};

template <std::floating_point Float>
WignerValues<Float>::WignerValues(int L, int M, int N,
                                  Float theta)
    : L{L}, M{M}, N{N}, theta{theta} {
  // Check the maximum degree is non-negative.
  assert(L >= 0);

  // Check the maximum order is in range.
  assert(M >= 0 && M <= L);

  // Check the upper index is in range.
  assert(N >= 0 && N <= L);

  // If N > M swap their values
  if(N > M) std::swap(M,N);

  // Reserve space for the values.
  data.reserve(ElementsIncludingUpperIndex(N));

  // Compute and store trigonometric terms
  Float cos = std::cos(theta);
  Float sinHalf = std::sin(theta / 2);
  Float cosHalf = std::cos(theta / 2);
  using std::log;
  if (sinHalf > std::numeric_limits<Float>::min()) {
    logSinHalf = log(sinHalf);
  } else {
    atLeft = true;
  }
  if (cosHalf > std::numeric_limits<Float>::min()) {
    logCosHalf = log(cosHalf);
  } else {
    atRight = true;
  }

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
}

// Returns the value of the generalised Legendre function, P_{l,m}^{l}(theta),
// following the definitions in Appendix C of Dahlen & Tromp (1998).
template <std::floating_point Float>
Float WignerValues<Float>::ValueMaxUpperIndexAtOrder(int l, int m) {
  // Deal with l == 0 case
  if (l == 0) return static_cast<Float>(1);

  // Deal with special case at the left boundary.
  if (atLeft) {
    return m == l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with special case at the right boundary.
  if (atRight) {
    return -m == l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with the general case.
  Float Fl = static_cast<Float>(l);
  Float Fm = static_cast<Float>(m);
  using std::exp;
  using std::lgamma;
  return std::exp(
      0.5 * (lgamma(2 * Fl + 1) - lgamma(Fl - Fm + 1) - lgamma(Fl + Fm + 1)) +
      (Fl - Fm) * logSinHalf + (Fl + Fm) * logCosHalf);
}


template< std::floating_point Float>
Float WignerValues<Float>::ValueMaxOrderAtUpperIndex(int l, int n) {

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
  Float Fl = static_cast<Float>(l);
  Float Fn = static_cast<Float>(n);
  using std::exp;
  using std::lgamma;
  return Sign(l+n) * std::exp(
      0.5 * (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
      (Fl - Fn) * logSinHalf + (Fl + Fn) * logCosHalf);
}


template< std::floating_point Float>
Float WignerValues<Float>::ValueMinOrderAtUpperIndex(int l, int n) {

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
  Float Fl = static_cast<Float>(l);
  Float Fn = static_cast<Float>(n);
  using std::exp;
  using std::lgamma;
  return std::exp(
      0.5 * (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
      (Fl + Fn) * logSinHalf + (Fl - Fn) * logCosHalf);
}  


template<std::floating_point Float>
auto RecusionCoefficients(int l, int m, int n) {
  Float f1 = 0.0;
  Float f2 = 0.0;
  return std::pair(f1,f2);
}
  

  
template <std::floating_point Float>
Float WignerValues<Float>::Numerator1(int l, int m, int n) {
  Float Fl = static_cast<Float>(l);
  Float Fm = static_cast<Float>(m);
  Float Fn = static_cast<Float>(n);
  return (2 * Fl - 1) * (Fl * (Fl - 1) * cos - m * n);
}

template <std::floating_point Float>
Float WignerValues<Float>::Numerator2(int l, int m, int n) {
  Float Fl = static_cast<Float>(l);
  Float Fm = static_cast<Float>(m);
  Float Fn = static_cast<Float>(n);
  return -Fl * std::sqrt((Fl - Fn - 1) * (Fl + Fn - 1) * (Fl - Fm - 1) *
                         (Fl + Fm - 1));
}

template <std::floating_point Float>
Float WignerValues<Float>::Denominator(int l, int m, int n) {
  Float Fl = static_cast<Float>(l);
  Float Fm = static_cast<Float>(m);
  Float Fn = static_cast<Float>(n);
  return (Fl - 1) * std::sqrt((Fl - Fn) * (Fl + Fn) * (Fl - Fm) * (Fl + Fm));
}

}  // namespace GSHT

#endif  // WIGNER_GUARD_H
