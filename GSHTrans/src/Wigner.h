#ifndef GSH_TRANS_WIGNER_GUARD_H
#define GSH_TRANS_WIGNER_GUARD_H

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <execution>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <vector>

#include "Concepts.h"

namespace GSHTrans {

// Define some tag-classes.
struct All {};
struct NonNegative {};

// Define some useful concepts.
template <typename Orders>
concept OrderRange =
    std::same_as<Orders, All> or std::same_as<Orders, NonNegative>;

// Define enum class for normalisation options.
enum class Normalisation { FourPi, Ortho };

// Declare/Define some utility functions.
constexpr int MinusOneToPower(int m) { return m % 2 ? -1 : 1; }

template <std::floating_point Real>
Real WignerMinOrderAtUpperIndex(int, int, Real, Real, bool, bool);

template <std::floating_point Real>
Real WignerMaxOrderAtUpperIndex(int, int, Real, Real, bool, bool);

template <std::floating_point Real>
Real WignerMaxUpperIndexAtOrder(int, int, Real, Real, bool, bool);

template <std::floating_point Real>
Real WignerMinUpperIndexAtOrder(int, int, Real, Real, bool, bool);

template <OrderRange Orders = All>
std::size_t WignerNStorage(int lMax, int mMax, int n);

template <OrderRange Orders>
requires std::same_as<Orders, All> std::size_t WignerNStorage(int lMax,
                                                              int mMax, int n) {
  auto nabs = std::abs(n);
  assert(lMax >= 0);
  assert(nabs >= 0 && nabs <= lMax);
  assert(mMax >= 0 && mMax <= lMax);
  return (mMax + 1) * (mMax + 1) - nabs * nabs + (lMax - mMax) * (2 * mMax + 1);
}

template <OrderRange Orders>
requires std::same_as<Orders, NonNegative> std::size_t WignerNStorage(int lMax,
                                                                      int mMax,
                                                                      int n) {
  auto nabs = std::abs(n);
  assert(lMax >= 0);
  assert(nabs >= 0 && nabs <= lMax);
  assert(mMax >= 0 && mMax <= lMax);
  return ((mMax + 1) * (mMax + 2)) / 2 - (nabs * (nabs + 1)) / 2 +
         (lMax - mMax) * (mMax + 1);
}

template <std::floating_point Real>
Real Wigner(int, int, int, Real, Normalisation);

/////////////////////////////////////////////////////////////////////////
//                         WignerN class                               //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Real, OrderRange Orders = All>
class WignerN {
 public:
  // Set member types.
  using value_type = Real;
  using iterator = Real *;
  using const_iterator = Real const *;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

  // Default constructor.
  WignerN() = default;

  // Copy constructor.
  WignerN(WignerN const &) = default;

  // Move constructor.
  WignerN(WignerN &&) = default;

  // Constructor with storage passed as iterators.
  template <RealFloatingPointIterator Iterator>
  WignerN(Iterator start, Iterator finish, int lMax, int mMax, int n,
          Real theta, Normalisation norm = Normalisation::Ortho)
      : _start{&*start}, _finish{&*finish}, _lMax{lMax}, _mMax{mMax}, _n{n} {
    ComputeValues(theta, norm);
  }

  // Constructor with storage passed as range.
  template <RealFloatingPointRange Range>
  WignerN(Range storage, int lMax, int mMax, int n, Real theta,
          Normalisation norm = Normalisation::Ortho)
      : _start{&*std::begin(storage)},
        _finish{&*std::end(storage)},
        _lMax{lMax},
        _mMax{mMax},
        _n{n} {
    ComputeValues(theta, norm);
  }

  // Copy assigment.
  WignerN &operator=(WignerN const &) = default;

  // Move assignment.
  WignerN &operator=(WignerN &&) = default;

  // Geters for basic data.
  int MaxDegree() const { return _lMax; }
  int MaxOrder() const { return _mMax; }
  int UpperIndex() const { return _n; }

  // Returns lowest order at given degree.
  int StartingOrder(int l) requires std::same_as<Orders, All> {
    assert(l <= _lMax && l >= std::abs(_n));
    return -std::min(l, _mMax);
  }
  int StartingOrder(int l) requires std::same_as<Orders, NonNegative> {
    assert(l <= _lMax && l >= std::abs(_n));
    return 0;
  }

  // Returns the number of values at a given degree when
  // all orders are stored.
  constexpr size_type Count(int l) const requires std::same_as<Orders, All> {
    auto nabs = std::abs(_n);
    if (l < nabs) return 0;
    if (l <= _mMax) return (l + 1) * (l + 1) - nabs * nabs;
    return (_mMax + 1) * (_mMax + 1) - nabs * nabs +
           (l - _mMax) * (2 * _mMax + 1);
  }

  // Returns the number of values at a given degree when
  // only non-negative orders are stored.
  constexpr size_type Count(
      int l) const requires std::same_as<Orders, NonNegative> {
    auto nabs = std::abs(_n);
    if (l < nabs) return 0;
    if (l <= _mMax) return ((l + 1) * (l + 2)) / 2 - (nabs * (nabs + 1)) / 2;
    return ((_mMax + 1) * (_mMax + 2)) / 2 - (nabs * (nabs + 1)) / 2 +
           (l - _mMax) * (_mMax + 1);
  }

  // Returns total number of values.
  constexpr size_type const Count() { return Count(_lMax); }

  // Iterators that point to the start of the data.
  iterator begin() { return _start; }
  const_iterator cbegin() const { return _start; }

  // Iterators that point to the end of the data.
  iterator end() { return _finish; }
  const_iterator cend() const { return _finish; }

  // Iterators that point to the start of degree l.
  iterator begin(int l) { return std::next(begin(), Count(l - 1)); }
  const_iterator cbegin(int l) const {
    return std::next(cbegin(), Count(l - 1));
  }

  // Iterators that point to the end of degree l.
  iterator end(int l) { return std::next(begin(), Count(l)); }
  const_iterator cend(int l) const { return std::next(cbegin(), Count(l)); }

  // Returns value for given degree and order when all orders are stored.
  auto operator()(int l, int m) const requires std::same_as<Orders, All> {
    assert(l >= std::abs(_n) && l <= _lMax);
    auto mMaxAbs = std::min(l, _mMax);
    assert(std::abs(m) <= mMaxAbs);
    //    std::cout << _start << std::endl;
    return *std::next(cbegin(l), mMaxAbs + m);
  }

  // Returns value for given degree and order, m >= 0, when only non-negative
  // orders are stored.
  auto operator()(int l,
                  int m) const requires std::same_as<Orders, NonNegative> {
    assert(l >= std::abs(_n) && l <= _lMax);
    assert(0 <= m && m <= std::min(l, _mMax));
    //    std::cout << _start << std::endl;
    return *std::next(cbegin(l), m);
  }

 private:
  // Set the execution policy.
  static constexpr auto _policy = std::execution::seq;

  int _lMax;  // Maximum degree.
  int _mMax;  // Maximum order.
  int _n;     // Upper index.

  // Store iterators to the data.
  iterator _start;
  iterator _finish;

  // Storage passed as iterators.
  void ComputeValues(Real theta, Normalisation norm);
};

template <std::floating_point Real, OrderRange Orders>
void WignerN<Real, Orders>::ComputeValues(Real theta, Normalisation norm) {
  // Check the maximum degree is non-negative.
  assert(_lMax >= 0);

  // Check the maximum order is in range.
  assert(_mMax >= 0 && _mMax <= _lMax);

  // Check storage space is appropriate
  assert(std::distance(_start, _finish) >= Count());

  // Pre-compute and store trigonometric terms.
  auto cos = std::cos(theta);
  auto logSinHalf = std::sin(0.5 * theta);
  auto logCosHalf = std::cos(0.5 * theta);
  auto atLeft = logSinHalf < std::numeric_limits<value_type>::min();
  auto atRight = logCosHalf < std::numeric_limits<value_type>::min();
  logSinHalf = atLeft ? static_cast<value_type>(0) : std::log(logSinHalf);
  logCosHalf = atRight ? static_cast<value_type>(0) : std::log(logCosHalf);

  // Pre-compute and store square roots and their inverses up to lMax + mMax.
  std::vector<value_type> sqInt(_lMax + _mMax + 1);
  std::transform(_policy, sqInt.begin(), sqInt.end(), sqInt.begin(),
                 [&](auto &x) {
                   return std::sqrt(static_cast<value_type>(&x - &sqInt[0]));
                 });
  std::vector<value_type> sqIntInv(_lMax + _mMax + 1);
  std::transform(
      _policy, sqInt.begin(), sqInt.end(), sqIntInv.begin(), [](auto x) {
        return x > static_cast<value_type>(0) ? 1 / x
                                              : static_cast<value_type>(0);
      });

  // Set the values for l == |n|
  const int nabs = std::abs(_n);
  {
    auto l = nabs;
    auto mStart = StartingOrder(l);
    auto start = begin(l);
    auto finish = end(l);
    if (_n >= 0) {
      std::transform(_policy, start, finish, start, [&](auto &p) {
        int m = mStart + std::distance(&*start, &p);
        return WignerMaxUpperIndexAtOrder(l, m, logSinHalf, logCosHalf, atLeft,
                                          atRight);
      });
    } else {
      std::transform(_policy, start, finish, start, [&](auto &p) {
        int m = mStart + std::distance(&*start, &p);
        return WignerMinUpperIndexAtOrder(l, m, logSinHalf, logCosHalf, atLeft,
                                          atRight);
      });
    }
  }

  // Set the values for l == n+1 if needed.
  if (nabs < _lMax) {
    auto l = nabs + 1;
    auto mStart = StartingOrder(l);

    // Set iterators
    auto startMinusOne = begin(l - 1);
    auto finishMinusOne = end(l - 1);
    auto start = begin(l);

    // Add in value at m == -l if needed.
    if constexpr (std::same_as<Orders, All>) {
      if (l <= _mMax) {
        *start++ = WignerMinOrderAtUpperIndex(l, _n, logSinHalf, logCosHalf,
                                              atLeft, atRight);
        // Update the starting order for recursion
        mStart += 1;
      }
    }

    // Add in interior orders using one-term recursion.
    {
      auto alpha = (2 * l - 1) * l * cos * sqIntInv[l + nabs];
      auto beta = (2 * l - 1) * sqIntInv[l + nabs];
      if (_n < 0) beta *= -1;
      std::transform(
          _policy, startMinusOne, finishMinusOne, start, [&](auto &minusOne) {
            int m = mStart + std::distance(&*startMinusOne, &minusOne);
            auto f1 = (alpha - beta * m) * sqIntInv[l - m] * sqIntInv[l + m];
            return f1 * minusOne;
          });
    }

    // Add in value at m == l if needed
    if (l <= _mMax) {
      auto mStep = 2 * l - 1;
      if constexpr (std::same_as<Orders, NonNegative>) {
        mStep = l;
      }
      *std::next(start, mStep) = WignerMaxOrderAtUpperIndex(
          l, _n, logSinHalf, logCosHalf, atLeft, atRight);
    }
  }

  // Now do the remaining degrees.
  for (int l = nabs + 2; l <= _lMax; l++) {
    // Starting order within two-term recursion
    auto mStart = StartingOrder(l);

    // Set iterators
    auto startMinusTwo = begin(l - 2);
    auto finishMinusTwo = end(l - 2);
    auto startMinusOne = begin(l - 1);
    auto start = begin(l);

    // Add in lower boundary terms if still growing.
    if constexpr (std::same_as<Orders, All>) {
      if (l <= _mMax) {
        // Add in the m == -l term.
        *start++ = WignerMinOrderAtUpperIndex(l, _n, logSinHalf, logCosHalf,
                                              atLeft, atRight);
        // Now do the m == -l+1 term using one-point recursion.
        {
          auto m = -l + 1;
          auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * _n) *
                    sqIntInv[l - _n] * sqIntInv[l + _n] * sqIntInv[l - m] *
                    sqIntInv[l + m] / static_cast<value_type>(l - 1);
          *start++ = f1 * (*startMinusOne++);
        }
        // Update the starting order for two-term recursion.
        mStart += 2;
      }

      // Add in the lower boundary term at the critical degree
      if (l == _mMax + 1) {
        auto m = -_mMax;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * _n) *
                  sqIntInv[l - _n] * sqIntInv[l + _n] * sqIntInv[l - m] *
                  sqIntInv[l + m] / static_cast<value_type>(l - 1);
        *start++ = f1 * (*startMinusOne++);
        // Update the starting order for two-term recursion.
        mStart += 1;
      }
    }

    // Apply two-term recusion for the interior orders.
    {
      auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - _n] * sqIntInv[l + _n];
      auto beta = (2 * l - 1) * _n * sqIntInv[l - _n] * sqIntInv[l + _n] /
                  static_cast<value_type>(l - 1);
      auto gamma = l * sqInt[l - 1 - _n] * sqInt[l - 1 + _n] *
                   sqIntInv[l - _n] * sqIntInv[l + _n] /
                   static_cast<value_type>(l - 1);
      std::transform(
          _policy, startMinusTwo, finishMinusTwo, startMinusOne, start,
          [&](auto &minusTwo, auto &minusOne) {
            int m = mStart + std::distance(&*startMinusTwo, &minusTwo);
            auto denom = sqIntInv[l - m] * sqIntInv[l + m];
            auto f1 = (alpha - beta * m) * denom;
            auto f2 = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
            return f1 * minusOne - f2 * minusTwo;
          });
    }

    // Add in the upper boundary terms if still growing.
    if (l <= _mMax) {
      // Update the iterator
      auto mStep = 2 * l - 3;
      if constexpr (std::same_as<Orders, NonNegative>) {
        mStep = l - 1;
      }
      startMinusOne = std::next(startMinusOne, mStep);
      start = std::next(start, mStep);
      // Add in m == l - 1 term using one-point recursion.
      {
        auto m = l - 1;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * _n) *
                  sqIntInv[l - _n] * sqIntInv[l + _n] * sqIntInv[l - m] *
                  sqIntInv[l + m] / static_cast<value_type>(l - 1);
        *start++ = f1 * (*startMinusOne++);
      }
      // Now do m == l.
      *start++ = WignerMaxOrderAtUpperIndex(l, _n, logSinHalf, logCosHalf,
                                            atLeft, atRight);
    }

    // Add in the upper boundary term at the crtiical degree.
    if (l == _mMax + 1) {
      // Update the iterators.
      auto mStep = 2 * _mMax - 1;
      if constexpr (std::same_as<Orders, NonNegative>) {
        mStep = _mMax;
      }
      startMinusOne = std::next(startMinusOne, mStep);
      start = std::next(start, mStep);
      auto m = _mMax;
      auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * _n) * sqIntInv[l - _n] *
                sqIntInv[l + _n] * sqIntInv[l - m] * sqIntInv[l + m] /
                static_cast<value_type>(l - 1);
      *start++ = f1 * (*startMinusOne++);
    }
  }

  if (norm == Normalisation::Ortho) {
    for (int l = 0; l <= _lMax; l++) {
      auto start = begin(l);
      auto finish = end(l);
      std::transform(start, finish, start, [l](auto p) {
        return 0.5 * std::sqrt(static_cast<value_type>(2 * l + 1)) *
               std::numbers::inv_sqrtpi_v<value_type> * p;
      });
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//      Simple function to compute Wigner values using +/- recursion in m //
////////////////////////////////////////////////////////////////////////////////

template <std::floating_point Real>
Real Wigner(int l, int m, int n, Real theta, Normalisation norm) {
  // Check the inputs.
  assert(l >= 0);
  assert(std::abs(m) <= l);
  assert(std::abs(n) <= l);

  // Deal with l = 0 separately.
  if (l == 0) {
    if (norm == Normalisation::Ortho) {
      return 0.5 * std::numbers::inv_sqrtpi_v<Real>;
    } else {
      return 1;
    }
  }

  // Pre-compute some trigonometric terms.
  auto logSinHalf = std::sin(0.5 * theta);
  auto logCosHalf = std::cos(0.5 * theta);
  auto atLeft = logSinHalf < std::numeric_limits<Real>::min();
  auto atRight = logCosHalf < std::numeric_limits<Real>::min();

  // Deal with values at the end points
  if (atLeft) {
    if (n == -l) {
      if (norm == Normalisation::Ortho) {
        return 0.5 * std::sqrt(static_cast<Real>(2 * l + 1)) *
               std::numbers::inv_sqrtpi_v<Real>;
      } else {
        return 1;
      }
    } else {
      return 0;
    }
  }

  if (atRight) {
    if (n == l) {
      if (norm == Normalisation::Ortho) {
        return MinusOneToPower(l + n) * 0.5 *
               std::sqrt(static_cast<Real>(2 * l + 1)) *
               std::numbers::inv_sqrtpi_v<Real>;
      } else {
        return MinusOneToPower(l + n);
      }
    } else {
      return 0;
    }
  }

  // Compute remaining trigonometric terms
  logSinHalf = std::log(logSinHalf);
  logCosHalf = std::log(logCosHalf);
  auto cosec = static_cast<Real>(1) / std::sin(theta);
  auto cot = std::cos(theta) * cosec;

  // Pre-compute and store square roots and their inverses up to 2*l+1.
  std::vector<Real> sqInt(2 * l + 1);
  std::transform(sqInt.begin(), sqInt.end(), sqInt.begin(), [&](auto &x) {
    return std::sqrt(static_cast<Real>(&x - &sqInt[0]));
  });
  std::vector<Real> sqIntInv(2 * l + 1);
  std::transform(sqInt.begin(), sqInt.end(), sqIntInv.begin(), [](auto x) {
    return x > static_cast<Real>(0) ? 1 / x : static_cast<Real>(0);
  });

  // Compute the optimal meeting point for the recursion.
  int mOpt = n * std::cos(theta);

  // Apply upward recursion from m = -l.
  if (m <= mOpt) {
    Real minusOne = 0;
    Real current =
        WignerMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf, false, false);
    for (int mp = -l; mp < m; mp++) {
      current = 2 * (n * cosec - mp * cot) * sqIntInv[l - mp] *
                    sqIntInv[l + mp + 1] * current -
                sqInt[l + mp] * sqInt[l - mp + 1] * sqIntInv[l - mp] *
                    sqIntInv[l + mp + 1] * minusOne;
      minusOne = current;
    }
    if (norm = Normalisation::Ortho) {
      return 0.5 * std::sqrt(static_cast<Real>(2 * l + 1)) *
             std::numbers::inv_sqrtpi_v<Real> * current;

    } else {
      return current;
    }
  }

  // Apply downward recursion from m = l.
  if (m > mOpt) {
    Real plusOne = 0;
    Real current =
        WignerMaxOrderAtUpperIndex(l, n, logSinHalf, logCosHalf, false, false);
    for (int mp = l; mp > m; mp--) {
      current = 2 * (n * cosec - mp * cot) * sqIntInv[l + mp] *
                    sqIntInv[l - mp + 1] * current -
                sqInt[l - mp] * sqInt[l + mp + 1] * sqIntInv[l + mp] *
                    sqIntInv[l - mp + 1] * plusOne;
      plusOne = current;
    }
    if (norm = Normalisation::Ortho) {
      return 0.5 * std::sqrt(static_cast<Real>(2 * l + 1)) *
             std::numbers::inv_sqrtpi_v<Real> * current;

    } else {
      return current;
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//                           Utility functions                         //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Real>
Real WignerMinOrderAtUpperIndex(int l, int n, Real logSinHalf, Real logCosHalf,
                                bool atLeft, bool atRight) {
  // Check the inputs.
  assert(l >= 0);
  assert(std::abs(n) <= l);

  // Deal with l == 0 case
  if (l == 0) return static_cast<Real>(1);

  // Deal with special case at the left boundary.
  if (atLeft) {
    return n == -l ? static_cast<Real>(1) : static_cast<Real>(0);
  }

  // Deal with special case at the right boundary.
  if (atRight) {
    return n == l ? static_cast<Real>(1) : static_cast<Real>(0);
  }

  // Deal with the general case.
  auto Fl = static_cast<Real>(l);
  auto Fn = static_cast<Real>(n);
  using std::exp;
  using std::lgamma;
  return exp(
      static_cast<Real>(0.5) *
          (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
      (Fl + Fn) * logSinHalf + (Fl - Fn) * logCosHalf);
}

template <std::floating_point Real>
Real WignerMaxOrderAtUpperIndex(int l, int n, Real logSinHalf, Real logCosHalf,
                                bool atLeft, bool atRight) {
  return MinusOneToPower(n + l) * WignerMinOrderAtUpperIndex(l, -n, logSinHalf,
                                                             logCosHalf, atLeft,
                                                             atRight);
}

template <std::floating_point Real>
Real WignerMinUpperIndexAtOrder(int l, int m, Real logSinHalf, Real logCosHalf,
                                bool atLeft, bool atRight) {
  return WignerMaxOrderAtUpperIndex(l, -m, logSinHalf, logCosHalf, atLeft,
                                    atRight);
}

template <std::floating_point Real>
Real WignerMaxUpperIndexAtOrder(int l, int m, Real logSinHalf, Real logCosHalf,
                                bool atLeft, bool atRight) {
  return WignerMinOrderAtUpperIndex(l, -m, logSinHalf, logCosHalf, atLeft,
                                    atRight);
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_WIGNER_GUARD_H
