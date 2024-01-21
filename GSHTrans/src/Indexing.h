#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <ranges>

#include "Concepts.h"

namespace GSHTrans {

template <IndexRange MRange>
class GSHIndices {
  using Integer = std::ptrdiff_t;

 public:
  GSHIndices() = delete;

  // General constructor taking in lMax, mMax, and n.
  GSHIndices(Integer lMax, Integer mMax, Integer n)
      : _lMax{lMax}, _mMax{mMax}, _n{n}, _nAbs{std::abs(_n)} {
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
  GSHIndices(Integer lMax, Integer n) : GSHIndices(lMax, lMax, n) {}

  // Reduced constructor that assumes mMax = lMax and nMax = 0;
  GSHIndices(Integer lMax) : GSHIndices(lMax, lMax, 0) {}

  // Return iterators to the begining and end.
  const auto& begin() const { return *this; }
  const auto& end() const { return *this; }

  // Not equals required for terminating loops.
  auto operator!=(const GSHIndices&) const { return _l <= _lMax; }

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

  // Return the upper index and its absolute value.
  auto UpperIndex() const { return _n; }
  auto AbsUpperIndex() const { return _nAbs; }

  // Return the minimum and maximum degrees.
  auto MinDegree() const { return _nAbs; }
  auto MaxDegree() const { return _lMax; }

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
  Integer _n;
  Integer _nAbs;
  Integer _l;
  Integer _m;
};

template <RealOrComplexFloatingPoint Scalar, IndexRange MRange>
class GSHViewDegree {
 public:
  using value_type = Scalar;
  using iterator = Scalar*;
  using const_iterator = Scalar const*;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

 private:
  using Int = difference_type;

 public:
  GSHViewDegree() = delete;

  template <RealOrComplexFloatingPointIterator Iterator>
  GSHViewDegree(Int l, Int mMax, Iterator start)
      : _l{l}, _mMax{mMax}, _start{&*start} {
    assert(_l >= 0);
    assert(_mMax >= 0);
  }

  auto Degree() const { return _l; }

  auto MinOrder() const {
    if constexpr (std::same_as<MRange, All>) {
      return -std::min(_l, _mMax);
    } else {
      return Int{0};
    }
  }

  auto MaxOrder() const { return std::min(_l, _mMax); }

  auto Orders() const {
    return std::ranges::views::iota(MinOrder(), MaxOrder() + 1);
  }

  auto NonNegativeOrders() const {
    return std::ranges::views::iota(0, MaxOrder() + 1);
  }

  auto NegativeOrders() const {
    if constexpr (std::same_as<MRange, All>) {
      return std::ranges::views::iota(MinOrder(), 0);

    } else {
      return std::ranges::views::empty<Int>;
    }
  }

  auto size() const { return Orders().size(); }

  auto begin() { return _start; }
  auto end() { return std::next(_start, size()); }

  auto cbegin() const { return _start; }
  auto cend() const { return std::next(_start, size()); }

  auto operator()(Int m) const {
    if constexpr (std::same_as<MRange, All>) {
      return _start[_mMax + m];
    } else {
      return _start[m];
    }
  }

  auto& operator()(Int m) {
    if constexpr (std::same_as<MRange, All>) {
      return _start[_mMax + m];
    } else {
      return _start[m];
    }
  }

 private:
  Int _l;
  Int _mMax;
  iterator _start;
};

template <RealOrComplexFloatingPoint Scalar, IndexRange MRange>
class GSHView {
 public:
  using value_type = Scalar;
  using iterator = Scalar*;
  using const_iterator = Scalar const*;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

 private:
  using Int = difference_type;

 public:
  GSHView() = delete;

  template <RealOrComplexFloatingPointIterator Iterator>
  GSHView(Int lMax, Int mMax, Int n, Iterator start)
      : _indices{GSHIndices<MRange>(lMax, mMax, n)}, _start{&*start} {}

  auto size() const { return _indices.size(); }

  auto begin() { return _start; }
  auto end() { return std::next(_start, size()); }

  auto cbegin() const { return _start; }
  auto cend() const { return std::next(_start, size()); }

  auto MinDegree() const { return _indices.MinDegree(); }
  auto MaxDegree() const { return _indices.MaxDegree(); }
  auto Degrees() const {
    return std::ranges::views::iota(MinDegree(), MaxDegree() + 1);
  }

  auto operator()(Int l) {
    auto offset = _indices(l, _indices.MinOrder(l));
    auto start = std::next(_start, offset);
    return GSHViewDegree<Scalar, MRange>(l, _indices.MaxOrder(l), start);
  }

  auto Indices() const { return _indices; }

 private:
  GSHIndices<MRange> _indices;
  iterator _start;
};

template <RealOrComplexFloatingPoint Scalar, IndexRange MRange>
class GSHViewAngleRange {
 public:
  using value_type = Scalar;
  using iterator = Scalar*;
  using const_iterator = Scalar const*;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

 private:
  using Int = difference_type;

 public:
  GSHViewAngleRange() = delete;

  template <RealOrComplexFloatingPointIterator Iterator>
  GSHViewAngleRange(Int lMax, Int mMax, Int n, Int nTheta, Iterator start)
      : _lMax{lMax}, _mMax{mMax}, _n{n}, _nTheta{nTheta}, _start{&*start} {}

  auto size() const { GSHIndices<MRange>(_lMax, _mMax, _n).size() * _nTheta; }

  auto begin() { return _start; }
  auto end() { return std::next(_start, size()); }

  auto cbegin() const { return _start; }
  auto cend() const { return std::next(_start, size()); }

  auto NumberOfAngles() const { return _nTheta; }
  auto AngleIndices() const { return std::ranges::views::iota(0, _nTheta); }

  auto operator()(Int iTheta = 0) {
    auto offset = GSHIndices<MRange>(_lMax, _mMax, _n).size() * iTheta;
    auto start = std::next(_start, offset);
    return GSHView<Scalar, MRange>(_lMax, _mMax, _n, start);
  }

 private:
  Int _lMax;
  Int _mMax;
  Int _n;
  Int _nTheta;
  iterator _start;
};

class GLGIndices {
  using Int = std::ptrdiff_t;

 public:
  GLGIndices() = default;

  GLGIndices(Int lMax)
      : _nTheta{lMax + 1}, _nPhi{2 * lMax}, _iTheta{0}, _iPhi{0} {}

  // Return iterators to the begining and end.
  const auto& begin() const { return *this; }
  const auto& end() const { return *this; }

  // Increment operator.
  void operator++() {
    if (++_iPhi == _nPhi) {
      _iPhi = 0;
      ++_iTheta;
    }
  }

  // Dereference by returning (l,m) as a pair.
  auto operator*() const { return std::pair(_iTheta, _iPhi); }

  // Not equals required for terminating loops.
  auto operator!=(const GLGIndices&) const { return _iTheta < _nTheta; }

 private:
  Int _iPhi;
  Int _iTheta;
  Int _nTheta;
  Int _nPhi;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H
