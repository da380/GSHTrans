#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

#include <cassert>
#include <ranges>

#include "Concepts.h"

namespace GSHTrans {

// Return view to the degrees for a given maximum value.
auto Degrees(int lMax) { return std::ranges::views::iota(0, lMax + 1); }

// View to the orders for a given degree
template <OrderRange ORange>
auto Orders(int l) {
  if constexpr (std::same_as<ORange, All>) {
    return std::ranges::views::iota(-l, l + 1);
  } else {
    return std::ranges::views::iota(0, l + 1);
  }
}

// Return view to the upper indices for given nMax.
template <OrderRange NRange>
auto UpperIndices(int nMax) {
  if constexpr (std::same_as<NRange, All>) {
    return std::ranges::views::iota(-nMax, nMax + 1);
  } else {
    return std::ranges::views::iota(0, nMax + 1);
  }
}

// Class for indexing spherical harmonic arrays based on
// the "natural" ordering.
template <OrderRange ORange = All>
class SphericalHarmonicIndices {
  using Integer = std::ptrdiff_t;

 public:
  SphericalHarmonicIndices() = delete;

  // General constructor taking in lMax, mMax, and n.
  SphericalHarmonicIndices(Integer lMax, Integer mMax, Integer n)
      : _lMax{lMax}, _mMax{mMax}, _nAbs{std::abs(n)} {
    assert(_lMax >= 0);
    assert(_mMax <= _lMax);
    assert(_nAbs <= _lMax);
    _l = _nAbs;
    if constexpr (std::same_as<ORange, All>) {
      _m = -std::min(_l, _mMax);
    } else {
      _m = 0;
    }
  }

  // Reduced constructor that assumes that mMax = lMax;
  SphericalHarmonicIndices(Integer lMax, Integer n)
      : SphericalHarmonicIndices(lMax, lMax, n) {}

  // Return iterators to the begining and end.
  const auto& begin() const { return *this; }
  const auto& end() const { return *this; }

  // Not equals required for terminating loops.
  auto operator!=(const SphericalHarmonicIndices&) const { return _l <= _lMax; }

  // Increment operator.
  void operator++() {
    if (_m == std::min(_l, _mMax)) {
      ++_l;
      if constexpr (std::same_as<ORange, All>) {
        _m = -std::min(_l, _mMax);
      } else {
        _m = 0;
      }
    } else {
      ++_m;
    }
  }

  // Dereference by returning (l,m) as a pair.
  auto operator*() const { return std::pair(_l, _m); }

  // Application operator returns the index of the (l,m)th value.
  auto operator()(Integer l,
                  Integer m) const requires std::same_as<ORange, All> {
    assert(l >= _nAbs && l <= _lMax);
    assert(m >= -l && m <= l);
    if (_mMax >= _nAbs) {
      return l <= _mMax ? l * (l + 1) + m - _nAbs * _nAbs
                        : (_mMax + 1) * (_mMax + 1) - _nAbs * _nAbs +
                              (l - 1 - _mMax) * (2 * _mMax + 1) + _mMax + m;
    } else {
      return (l - _nAbs) * (2 * _mMax + 1) + _mMax + m;
    }
  }
  auto operator()(Integer l,
                  Integer m) const requires std::same_as<ORange, NonNegative> {
    assert(l >= _nAbs && l <= _lMax);
    assert(m >= 0 && m <= l);
    if (_mMax >= _nAbs) {
      return l <= _mMax
                 ? (l * (l + 1)) / 2 + m - (_nAbs * (_nAbs + 1)) / 2
                 : ((_mMax + 1) * (_mMax + 2)) / 2 - (_nAbs * (_nAbs + 1)) / 2 +
                       (l - 1 - _mMax) * (_mMax + 1) + m;
    } else {
      return (l - _nAbs) * (_mMax + 1) + m;
    }
  }

  // Return the total number of (l,m) values.
  auto size() const { return operator()(_lMax, _mMax) + 1; }

 private:
  Integer _lMax;
  Integer _mMax;
  Integer _nAbs;
  Integer _l;
  Integer _m;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H
