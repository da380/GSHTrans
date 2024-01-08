#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

#include <algorithm>
#include <cassert>
#include <numeric>
#include <ranges>

#include "Concepts.h"

namespace GSHTrans {

// Return view to the degrees for a given maximum value.
auto Degrees(int lMax) { return std::ranges::views::iota(0, lMax + 1); }

// View to the orders for a given degree
template <IndexRange MRange>
auto Orders(int mMax) {
  if constexpr (std::same_as<MRange, All>) {
    return std::ranges::views::iota(-mMax, mMax + 1);
  } else {
    return std::ranges::views::iota(0, mMax + 1);
  }
}

// Return view to the upper indices for given nMax.
template <IndexRange NRange>
auto UpperIndices(int nMax) {
  if constexpr (std::same_as<NRange, All>) {
    return std::ranges::views::iota(-nMax, nMax + 1);
  } else {
    return std::ranges::views::iota(0, nMax + 1);
  }
}

template <IndexRange MRange = All>
class SHIndices {
  using Integer = std::ptrdiff_t;

 public:
  SHIndices() = delete;

  // General constructor taking in lMax, mMax, and n.
  SHIndices(Integer lMax, Integer mMax, Integer n)
      : _lMax{lMax}, _mMax{mMax}, _nAbs{std::abs(n)} {
    assert(_lMax >= 0);
    assert(_mMax <= _lMax);
    assert(_nAbs <= _lMax);
    _l = _nAbs;
    if constexpr (std::same_as<MRange, All>) {
      _m = -std::min(_l, _mMax);
    } else {
      _m = 0;
    }
  }

  // Reduced constructor that assumes that mMax = lMax;
  SHIndices(Integer lMax, Integer n) : SHIndices(lMax, lMax, n) {}

  // Reduced constructor that assumes mMax = lMax and nMax = 0;
  SHIndices(Integer lMax) : SHIndices(lMax, lMax, 0) {}

  // Return iterators to the begining and end.
  const auto& begin() const { return *this; }
  const auto& end() const { return *this; }

  // Not equals required for terminating loops.
  auto operator!=(const SHIndices&) const { return _l <= _lMax; }

  // Increment operator.
  void operator++() {
    if (_m == std::min(_l, _mMax)) {
      ++_l;
      if constexpr (std::same_as<MRange, All>) {
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
                  Integer m) const requires std::same_as<MRange, All> {
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
                  Integer m) const requires std::same_as<MRange, NonNegative> {
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

  // Return minimum order for given degree.
  auto MinOrder(Integer l) const {
    if constexpr (std::same_as<MRange, All>) {
      return -std::min(l, _mMax);
    } else {
      return 0;
    }
  }

  // Return maximum order for given degree.
  auto MaxOrder(Integer l) const { return std::min(l, _mMax); }

  // Return offset for values at degree l.
  auto OffsetForDegree(Integer l) const {
    return l == 0 ? 0 : operator()(l - 1, MinOrder(l));
  }

  // Return size for value at degree l.
  auto sizeForDegree(Integer l) const {
    if constexpr (std::same_as<MRange, All>) {
      if (l < _mMax) {
        return 2 * l + 1;
      } else {
        return 2 * _mMax + 1;
      }
    } else {
      if (l < _mMax) {
        return l + 1;
      } else {
        return _mMax + 1;
      }
    }
  }

 private:
  Integer _lMax;
  Integer _mMax;
  Integer _nAbs;
  Integer _l;
  Integer _m;
};

template <IndexRange MRange, RealOrComplexFloatingPointIterator Iterator>
class SHViewFixedNL {
  using Integer = std::iterator_traits<Iterator>::difference_type;

 public:
  SHViewFixedNL() = default;

  SHViewFixedNL(Integer l, Integer mMax, Iterator start)
      : _l{l}, _mMax{mMax}, _start{start} {
    assert(_l >= 0);
    assert(_mMax >= 0);
  }

  auto size() const {
    if constexpr (std::same_as<MRange, All>) {
      return _l > _mMax ? 2 * _l + 1 : 2 * _mMax + 1;
    } else {
      return _l > _mMax ? _l + 1 : _mMax + 1;
    }
  }

  auto begin() { return _start; }
  auto end() { return std::next(_start, size()); }

  auto operator()(Integer m) const {
    if constexpr (std::same_as<MRange, All>) {
      return _start[_mMax + m];
    } else {
      return _start[m];
    }
  }

  auto Orders() const {
    auto m = std::min(_l, _mMax);
    if constexpr (std::same_as<MRange, All>) {
      return std::ranges::views::iota(-m, m + 1);
    } else {
      return std::ranges::views::iota(0, m + 1);
    }
  }

 private:
  Integer _l;
  Integer _mMax;
  Iterator _start;
};

template <IndexRange MRange, RealOrComplexFloatingPointIterator Iterator>
class SHViewFixedN {
  using Integer = std::iterator_traits<Iterator>::difference_type;

 public:
  SHViewFixedN() = default;

  SHViewFixedN(Integer lMax, Integer mMax, Integer n, Iterator start)
      : _indices{SHIndices<MRange>(lMax, mMax, n)}, _start{start} {}

  auto size() const { return _indices.size(); }
  auto begin() { return _start; }
  auto end() { return std::next(_start, size()); }

  auto operator()(Integer l) {
    auto offset = _indices(l, _indices.MinOrder(l));
    auto start = std::next(_start, Offset(l));
    return SHViewFixedNL<MRange, Iterator>(l, _indices.MaxOrder(l), start);
  }

  auto Indices() const { return _indices; }

 private:
  SHIndices<MRange> _indices;
  Iterator _start;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H
