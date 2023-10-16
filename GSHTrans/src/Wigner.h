#ifndef GSH_TRANS_WIGNER_GUARD_H
#define GSH_TRANS_WIGNER_GUARD_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <vector>

namespace GSHTrans {

// Define some tag-classes.
struct All {};
struct NonNegative {};

// Define some useful concepts.
template <typename Range>
concept IndexRange =
    std::same_as<Range, All> or std::same_as<Range, NonNegative>;

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

/////////////////////////////////////////////////////////////////////////
//                      WignerArrayN class                             //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Real, IndexRange Range = All>
class WignerArrayN {
 public:
  // Define member types.
  using value_type = Real;
  using iterator = std::vector<Real>::iterator;
  using const_iterator = std::vector<Real>::const_iterator;
  using difference_type = std::vector<Real>::difference_type;

  // Constructor
  WignerArrayN(int, int, int, Real, Normalisation norm = Normalisation::Ortho);

  // Geters for basic data.
  int MaxDegree() const { return _lMax; }
  int MaxOrder() const { return _mMax; }
  int UpperIndex() const { return _n; }

  // Returns lowest order at given degree.
  int StartingOrder(int l) requires std::same_as<Range, All> {
    assert(l <= _lMax && l >= std::abs(_n));
    return -std::min(l, _mMax);
  }
  int StartingOrder(int l) requires std::same_as<Range, NonNegative> {
    assert(l <= _lMax && l >= std::abs(_n));
    return 0;
  }

  // Iterators that point to the start of the data.
  iterator begin() { return _data.begin(); }
  const_iterator cbegin() const { return _data.cbegin(); }

  // Iterators that point to the end of the _data.
  iterator end() { return _data.end(); }
  const_iterator cend() const { return _data.cend(); }

  // Iterators that point to the start of degree l.
  iterator begin(int l) { return std::next(begin(), Count(l - 1)); }
  const_iterator cbegin(int l) const {
    return std::next(cbegin(), Count(l - 1));
  }

  // Iterators that point to the end of degree l.
  iterator end(int l) { return std::next(begin(), Count(l)); }
  const_iterator cend(int l) const { return std::next(cbegin(), Count(l)); }

  // Returns value for given degree and order when all orders are stored.
  Real operator()(int l, int m) const requires std::same_as<Range, All> {
    assert(l <= _lMax);
    auto mMaxAbs = std::min(l, _mMax);
    assert(std::abs(m) <= mMaxAbs);
    return *std::next(cbegin(l), mMaxAbs + m);
  }

  // Returns value for given degree and order, m >= 0, when only non-negative
  // orders are stored.
  Real operator()(int l,
                  int m) const requires std::same_as<Range, NonNegative> {
    assert(l >= std::abs(_n) && l <= _lMax);
    assert(0 <= m && m <= std::min(l, _mMax));
    return *std::next(cbegin(l), m);
  }

 private:
  // Set the execution policy.
  static constexpr auto _policy = std::execution::seq;

  int _lMax;  // Maximum degree.
  int _mMax;  // Maximum order.
  int _n;     // Upper index.

  std::vector<Real> _data;  // Vector containing values.

  // Returns the number of values at a given degree when
  // all orders are stored.
  constexpr difference_type Count(
      int l) const requires std::same_as<Range, All> {
    auto nabs = std::abs(_n);
    if (l < nabs) return 0;
    if (l <= _mMax) return (l + 1) * (l + 1) - nabs * nabs;
    return (_mMax + 1) * (_mMax + 1) - nabs * nabs +
           (l - _mMax) * (2 * _mMax + 1);
  }

  // Returns the number of values at a given degree when
  // only non-negative orders are stored.
  constexpr difference_type Count(
      int l) const requires std::same_as<Range, NonNegative> {
    auto nabs = std::abs(_n);
    if (l < nabs) return 0;
    if (l <= _mMax) return ((l + 1) * (l + 2)) / 2 - (nabs * (nabs + 1)) / 2;
    return ((_mMax + 1) * (_mMax + 2)) / 2 - (nabs * (nabs + 1)) / 2 +
           (l - _mMax) * (_mMax + 1);
  }

  // Returns total number of values.
  constexpr difference_type const Count() { return Count(_lMax); }
};

template <std::floating_point Real, IndexRange Range>
WignerArrayN<Real, Range>::WignerArrayN(int lMax, int mMax, int n, Real theta,
                                        Normalisation norm)
    : _lMax{lMax}, _mMax{mMax}, _n{n} {
  // Check the maximum degree is non-negative.
  assert(_lMax >= 0);

  // Check the maximum order is in range.
  assert(_mMax >= 0 && _mMax <= _lMax);

  // Initialise the data vector.
  _data = std::vector<Real>(Count());

  // Pre-compute and store trigonometric terms.
  auto cos = std::cos(theta);
  auto logSinHalf = std::sin(0.5 * theta);
  auto logCosHalf = std::cos(0.5 * theta);
  auto atLeft = logSinHalf < std::numeric_limits<Real>::min();
  auto atRight = logCosHalf < std::numeric_limits<Real>::min();
  logSinHalf = atLeft ? static_cast<Real>(0) : std::log(logSinHalf);
  logCosHalf = atRight ? static_cast<Real>(0) : std::log(logCosHalf);

  // Pre-compute and store square roots and their inverses up to lMax + mMax.
  std::vector<Real> sqInt(_lMax + _mMax + 1);
  std::transform(
      _policy, sqInt.begin(), sqInt.end(), sqInt.begin(),
      [&](auto &x) { return std::sqrt(static_cast<Real>(&x - &sqInt[0])); });
  std::vector<Real> sqIntInv(_lMax + _mMax + 1);
  std::transform(
      _policy, sqInt.begin(), sqInt.end(), sqIntInv.begin(), [](auto x) {
        return x > static_cast<Real>(0) ? 1 / x : static_cast<Real>(0);
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
    if constexpr (std::same_as<Range, All>) {
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
      if constexpr (std::same_as<Range, NonNegative>) {
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
    if constexpr (std::same_as<Range, All>) {
      if (l <= _mMax) {
        // Add in the m == -l term.
        *start++ = WignerMinOrderAtUpperIndex(l, _n, logSinHalf, logCosHalf,
                                              atLeft, atRight);
        // Now do the m == -l+1 term using one-point recursion.
        {
          auto m = -l + 1;
          auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * _n) *
                    sqIntInv[l - _n] * sqIntInv[l + _n] * sqIntInv[l - m] *
                    sqIntInv[l + m] / static_cast<Real>(l - 1);
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
                  sqIntInv[l + m] / static_cast<Real>(l - 1);
        *start++ = f1 * (*startMinusOne++);
        // Update the starting order for two-term recursion.
        mStart += 1;
      }
    }

    // Apply two-term recusion for the interior orders.
    {
      auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - _n] * sqIntInv[l + _n];
      auto beta = (2 * l - 1) * _n * sqIntInv[l - _n] * sqIntInv[l + _n] /
                  static_cast<Real>(l - 1);
      auto gamma = l * sqInt[l - 1 - _n] * sqInt[l - 1 + _n] *
                   sqIntInv[l - _n] * sqIntInv[l + _n] /
                   static_cast<Real>(l - 1);
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
      if constexpr (std::same_as<Range, NonNegative>) {
        mStep = l - 1;
      }
      startMinusOne = std::next(startMinusOne, mStep);
      start = std::next(start, mStep);
      // Add in m == l - 1 term using one-point recursion.
      {
        auto m = l - 1;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * _n) *
                  sqIntInv[l - _n] * sqIntInv[l + _n] * sqIntInv[l - m] *
                  sqIntInv[l + m] / static_cast<Real>(l - 1);
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
      if constexpr (std::same_as<Range, NonNegative>) {
        mStep = _mMax;
      }
      startMinusOne = std::next(startMinusOne, mStep);
      start = std::next(start, mStep);
      auto m = _mMax;
      auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * _n) * sqIntInv[l - _n] *
                sqIntInv[l + _n] * sqIntInv[l - m] * sqIntInv[l + m] /
                static_cast<Real>(l - 1);
      *start++ = f1 * (*startMinusOne++);
    }
  }

  if (norm == Normalisation::Ortho) {
    for (int l = 0; l <= _lMax; l++) {
      auto start = begin(l);
      auto finish = end(l);
      std::transform(start, finish, start, [l](auto p) {
        return 0.5 * std::sqrt(static_cast<Real>(2 * l + 1)) *
               std::numbers::inv_sqrtpi_v<Real> * p;
      });
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//                   Definition of WignerArray class                   //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Real, IndexRange MRange = All,
          IndexRange NRange = All>
class WignerArray {
 public:
  using value_type = Real;
  using iterator = WignerArrayN<Real, MRange>::iterator;
  using const_iterator = WignerArrayN<Real, MRange>::const_iterator;
  using difference_type = WignerArrayN<Real, MRange>::difference_type;
  template <int Storage>

  WignerArray(int lMax, int mMax, int nMax, Real theta,
              Normalisation norm = Normalisation::Ortho)
      : _lMax{lMax}, _mMax{mMax}, _nMax{nMax} {
    assert(0 <= _lMax);
    assert(0 <= _mMax && _mMax <= _lMax);
    assert(0 <= _nMax && _nMax <= _lMax);
    for (int n = nMin(); n <= _nMax; n++) {
      _data.push_back(std::make_unique<Array>(_lMax, _mMax, n, theta));
    }
  }

  int MaxDegree() const { return _lMax; }
  int MaxOrder() const { return _mMax; }
  int MaxUpperIndex() const { return _nMax; }

  // Return value for given degree, order and
  // upper index.
  Real operator()(int l, int m, int n) const {
    assert(0 <= l && l <= _lMax);
    assert(mMin() <= m && m <= _mMax);
    assert(nMin() <= n && n <= _nMax);
    return _data[nIndex(n)]->operator()(l, m);
  }

  // Iterators to the data for given upper index.
  iterator begin(int n) { return _data[nIndex(n)]->begin(); }
  iterator cbegin(int n) const { return _data[nIndex(n)]->cbegin(); }
  iterator end(int n) { return _data[nIndex(n)]->end(); }
  iterator cend(int n) const { return _data[nIndex(n)]->cend(); }

  // Iterators to the data for given degree and upper index.
  iterator begin(int l, int n) { return _data[nIndex(n)]->begin(l); }
  iterator cbegin(int l, int n) const { return _data[nIndex(n)]->cbegin(l); }
  iterator end(int l, int n) { return _data[nIndex(n)]->end(l); }
  iterator cend(int l, int n) const { return _data[nIndex(n)]->cend(l); }

 private:
  int _lMax;
  int _mMax;
  int _nMax;
  using Array = WignerArrayN<Real, MRange>;
  using ArrayPointer = std::unique_ptr<Array>;
  std::vector<ArrayPointer> _data;

  constexpr int mMin() const {
    if constexpr (std::same_as<MRange, All>) {
      return -_mMax;
    }
    if constexpr (std::same_as<MRange, NonNegative>) {
      return 0;
    }
  }

  constexpr int nMin() const {
    if constexpr (std::same_as<NRange, All>) {
      return -_nMax;
    }
    if constexpr (std::same_as<NRange, NonNegative>) {
      return 0;
    }
  }

  constexpr int nIndex(int n) const {
    if constexpr (std::same_as<NRange, All>) {
      return _nMax + n;
    }
    if constexpr (std::same_as<NRange, NonNegative>) {
      return n;
    }
  }

  constexpr int mIndex(int m) const {
    if constexpr (std::same_as<MRange, All>) {
      return _mMax + m;
    }
    if constexpr (std::same_as<MRange, NonNegative>) {
      return m;
    }
  }
};

/////////////////////////////////////////////////////////////////////////
//                          WignerArrayLN class                        //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Real>
class WignerArrayLN {
 public:
  // Define member types.
  using value_type = Real;
  using iterator = std::vector<Real>::iterator;
  using const_iterator = std::vector<Real>::const_iterator;
  using difference_type = std::vector<Real>::difference_type;

  // Constructor.
  WignerArrayLN(int, int, Real, Normalisation norm = Normalisation::Ortho);

  // Basic data functions.
  int Degree() const { return l; }
  int UpperIndex() const { return n; }
  Real Angle() const { return theta; }

  // Iterators to the data.
  iterator begin() { return _data.begin(); }
  const_iterator cbegin() const { return _data.cbegin(); }
  iterator end() { return _data.end(); }
  const_iterator cend() { return _data.end(); }

  // Returns value at given order.
  Real operator()(int m) const { return _data[l + m]; }

 private:
  int l;                    // Degree.
  int n;                    // Upper index.
  Real theta;               // Angle.
  std::vector<Real> _data;  // Stored values.
};

template <std::floating_point Real>
WignerArrayLN<Real>::WignerArrayLN(int l, int n, Real theta, Normalisation norm)
    : l{l}, n{n}, theta{theta} {
  // Check the inputs
  assert(0 <= l);
  assert(std::abs(n) <= l);

  // Deal with l = 0 separately.
  if (l == 0) {
    _data.push_back(1.0);
    if (norm == Normalisation::Ortho) {
      _data[0] *= 0.5 * std::numbers::inv_sqrtpi_v<Real>;
    }
    return;
  }

  // Allocate the data vector.
  _data = std::vector<Real>(2 * l + 1);

  // Pre-compute some trigonometric terms.
  auto logSinHalf = std::sin(0.5 * theta);
  auto logCosHalf = std::cos(0.5 * theta);
  auto atLeft = logSinHalf < std::numeric_limits<Real>::min();
  auto atRight = logCosHalf < std::numeric_limits<Real>::min();

  // Deal with values at the end points
  if (atLeft) {
    _data[l + n] = 1;
    if (norm == Normalisation::Ortho) {
      _data[l + n] *= 0.5 * std::sqrt(static_cast<Real>(2 * l + 1)) *
                      std::numbers::inv_sqrtpi_v<Real>;
    }
    return;
  }

  if (atRight) {
    _data[l - n] = MinusOneToPower(l + n);
    if (norm == Normalisation::Ortho) {
      _data[l - n] *= 0.5 * std::sqrt(static_cast<Real>(2 * l + 1)) *
                      std::numbers::inv_sqrtpi_v<Real>;
    }
    return;
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

  // Set value at minimum order directly.
  _data[0] =
      WignerMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf, false, false);

  //  Upwards recusion.
  Real minusOne = 0;
  for (int m = -l; m < mOpt; m++) {
    Real current = _data[l + m];
    _data[l + m + 1] = 2 * (n * cosec - m * cot) * sqIntInv[l - m] *
                           sqIntInv[l + m + 1] * current -
                       sqInt[l + m] * sqInt[l - m + 1] * sqIntInv[l - m] *
                           sqIntInv[l + m + 1] * minusOne;
    minusOne = current;
  }

  // Set value at maximum order directly.
  _data[2 * l] =
      WignerMaxOrderAtUpperIndex(l, n, logSinHalf, logCosHalf, false, false);

  // Downwards recusion.
  Real plusOne = 0;
  for (int m = l; m > mOpt + 1; m--) {
    Real current = _data[l + m];
    _data[l + m - 1] = 2 * (n * cosec - m * cot) * sqIntInv[l + m] *
                           sqIntInv[l - m + 1] * current -
                       sqInt[l - m] * sqInt[l + m + 1] * sqIntInv[l + m] *
                           sqIntInv[l - m + 1] * plusOne;
    plusOne = current;
  }

  if (norm == Normalisation::Ortho) {
    std::transform(_data.begin(), _data.end(), _data.begin(), [l](auto p) {
      return 0.5 * std::sqrt(static_cast<Real>(2 * l + 1)) *
             std::numbers::inv_sqrtpi_v<Real> * p;
    });
  }
}

/////////////////////////////////////////////////////////////////////////
//                Definition of some utility functions                 //
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

template <std::floating_point Real>
Real Winger(int l, int m, int n, Real theta, Normalisation norm) {
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

}  // namespace GSHTrans

#endif  // GSH_TRANS_WIGNER_GUARD_H
