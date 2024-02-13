#ifndef GSH_TRANS_INDEXING_NEW_GUARD_H
#define GSH_TRANS_INDEXING_NEW_GUARD_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>

#include "Concepts.h"

namespace GSHTrans {

namespace Testing {

template <OrderIndexRange MRange>
class GSHIndices {
  using Int = std::ptrdiff_t;

  class SubIndices {
   public:
    SubIndices(Int l, Int mMax) : _l{l}, _mMax{std::min(l, mMax)} {
      assert(_l >= 0);
      assert(_mMax >= 0);
    }

    auto Degree() const { return _l; }

    auto MinOrder() const {
      if constexpr (std::same_as<MRange, All>) {
        return -_mMax;
      } else {
        return Int{0};
      }
    }

    auto MaxOrder() const { return _mMax; }

    auto Orders() const {
      return std::ranges::views::iota(MinOrder(), MaxOrder() + 1);
    }

    auto NegativeOrders() const {
      return std::ranges::views::iota(MinOrder(), 0);
    }

    auto NonNegativeOrders() const {
      return std::ranges::views::iota(0, MaxOrder() + 1);
    }

    auto size() const { return MaxOrder() - MinOrder() + 1; }

    auto Index(Int m) const {
      if constexpr (std::same_as<MRange, All>) {
        return m + _mMax;
      } else {
        return m;
      }
    }

   private:
    Int _l;
    Int _mMax;
  };

 public:
  GSHIndices() = default;

  GSHIndices(Int lMax, Int mMax, Int n)
      : _lMax{lMax}, _mMax{std::min(lMax, mMax)}, _n{n} {
    assert(_lMax >= 0);
    assert(_mMax >= 0);
    assert(std::abs(n) <= _lMax);
  }

  auto UpperIndex() const { return _n; }

  auto MinDegree() const { return std::abs(_n); }
  auto MaxDegree() const { return _lMax; }
  auto Degrees() const {
    return std::ranges::views::iota(MinDegree(), MaxDegree() + 1);
  }

  auto OffsetForDegree(Int l) const
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

  auto OffsetForDegree(Int l) const
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

  auto SizeForDegree(Int l) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return SubIndices(l, _mMax).size();
  }

  auto Index(Int l) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return std::pair(OffsetForDegree(l), SubIndices(l, _mMax));
  }

 private:
  Int _lMax;
  Int _mMax;
  Int _n;
};

}  // namespace Testing

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_NEW_GUARD_H