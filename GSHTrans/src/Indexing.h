#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>
#include <utility>

#include "Concepts.h"

namespace GSHTrans {

//-----------------------------------------------//
//           Spherical harmonic indices          //
//-----------------------------------------------//

template <OrderIndexRange MRange>
class GSHSubIndices {
  using Int = std::ptrdiff_t;

 public:
  constexpr GSHSubIndices(Int l, Int mMax) : _l{l}, _mMax{std::min(l, mMax)} {
    assert(_l >= 0);
    assert(_mMax >= 0);
  }

  constexpr auto Degree() const { return _l; }

  constexpr auto MinOrder() const {
    if constexpr (std::same_as<MRange, All>) {
      return -_mMax;
    } else {
      return Int{0};
    }
  }

  constexpr auto MaxOrder() const { return _mMax; }

  constexpr auto Orders() const {
    return std::ranges::views::iota(MinOrder(), MaxOrder() + 1);
  }

  constexpr auto NegativeOrders() const {
    return std::ranges::views::iota(MinOrder(), 0);
  }

  constexpr auto NonNegativeOrders() const {
    return std::ranges::views::iota(0, MaxOrder() + 1);
  }

  constexpr auto size() const { return MaxOrder() - MinOrder() + 1; }

  constexpr auto Index(Int m) const {
    if constexpr (std::same_as<MRange, All>) {
      assert(m >= -_l && m <= _l);
      return m + _mMax;
    } else {
      assert(m >= 0 && m <= _l);
      return m;
    }
  }

 private:
  Int _l;
  Int _mMax;
};

template <OrderIndexRange MRange>
class GSHIndices {
  using Int = std::ptrdiff_t;

 public:
  GSHIndices() = default;
  constexpr GSHIndices(Int lMax, Int mMax, Int n)
      : _lMax{lMax}, _mMax{std::min(lMax, mMax)}, _n{n} {
    assert(_lMax >= 0);
    assert(_mMax >= 0);
    assert(std::abs(n) <= _lMax);
  }

  constexpr auto UpperIndex() const { return _n; }

  constexpr auto MaxOrder() const { return _mMax; }

  constexpr auto MinDegree() const { return std::abs(_n); }
  constexpr auto MaxDegree() const { return _lMax; }
  constexpr auto Degrees() const {
    return std::ranges::views::iota(MinDegree(), MaxDegree() + 1);
  }

  constexpr auto Indices() const {
    return Degrees() | std::ranges::views::transform([this](auto l) {
             return std::ranges::views::cartesian_product(
                 std::ranges::views::single(l),
                 GSHSubIndices<MRange>(l, _mMax).Orders());
           }) |
           std::ranges::views::join;
  }

  constexpr auto OffsetForDegree(Int l) const
  requires std::same_as<MRange, All>
  {
    assert(l >= MinDegree() && l <= MaxDegree());
    auto nAbs = std::abs(_n);
    if (_mMax >= nAbs) {
      return l <= _mMax ? l * l - nAbs * nAbs
                        : (_mMax + 1) * (_mMax + 1) - nAbs * nAbs +
                              (l - 1 - _mMax) * (2 * _mMax + 1);
    } else {
      return (l - nAbs) * (2 * _mMax + 1);
    }
  }

  constexpr auto OffsetForDegree(Int l) const
  requires std::same_as<MRange, NonNegative>
  {
    assert(l >= MinDegree() && l <= MaxDegree());
    auto nAbs = std::abs(_n);
    if (_mMax >= nAbs) {
      return l <= _mMax
                 ? (l * (l + 1)) / 2 - (nAbs * (nAbs + 1)) / 2
                 : ((_mMax + 1) * (_mMax + 2)) / 2 - (nAbs * (nAbs + 1)) / 2 +
                       (l - 1 - _mMax) * (_mMax + 1);
    } else {
      return (l - nAbs) * (_mMax + 1);
    }
  }

  constexpr auto SizeForDegree(Int l) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return GSHSubIndices<MRange>(l, _mMax).size();
  }

  constexpr auto Index(Int l) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return std::pair(OffsetForDegree(l), GSHSubIndices<MRange>(l, _mMax));
  }

  constexpr auto Index(Int l, Int m) const {
    auto [offset, indices] = Index(l);
    return offset + indices.Index(m);
  }

  constexpr auto size() const {
    return OffsetForDegree(_lMax) + SizeForDegree(_lMax);
  }

 private:
  Int _lMax;
  Int _mMax;
  Int _n;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H