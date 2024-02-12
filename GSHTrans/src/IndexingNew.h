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

template <OrderIndexRange MRange>
class CoefficientIndices {
  using Int = std::ptrdiff_t;

 public:
  CoefficientIndices(Int lMax, Int mMax, Int n)
      : _lMax{lMax}, _mMax{std::min(lMax, mMax)}, _n{n} {
    assert(_lMax >= 0);
    assert(_mMax >= 0);
    assert(std::abs(n) <= _lMax);
  }

  Int operator()(Int l, Int m) const
  requires std::same_as<MRange, All>
  {
    auto nAbs = std::abs(_n);
    if (_mMax >= nAbs) {
      return l <= _mMax ? l * (l + 1) + m - nAbs * nAbs
                        : (_mMax + 1) * (_mMax + 1) - nAbs * nAbs +
                              (l - 1 - _mMax) * (2 * _mMax + 1) + _mMax + m;
    } else {
      return (l - nAbs) * (2 * _mMax + 1) + _mMax + m;
    }
  }

  Int operator()(Int l, Int m) const
  requires std::same_as<MRange, NonNegative>
  {
    auto nAbs = std::abs(_n);
    if (_mMax >= nAbs) {
      return l <= _mMax
                 ? (l * (l + 1)) / 2 + m - (nAbs * (nAbs + 1)) / 2
                 : ((_mMax + 1) * (_mMax + 2)) / 2 - (nAbs * (nAbs + 1)) / 2 +
                       (l - 1 - _mMax) * (_mMax + 1) + m;
    } else {
      return (l - nAbs) * (_mMax + 1) + m;
    }
  }

  auto size() const { return operator()(_lMax, _mMax) + 1; }

 private:
  Int _lMax;
  Int _mMax;
  Int _n;
};

//------------------------------------------------------------//
//                View to data at a single degree             //
//------------------------------------------------------------//

template <std::ranges::view View, OrderIndexRange MRange>
class CoefficientSubView
    : public std::ranges::view_interface<CoefficientSubView<View, MRange>> {
  using Int = std::ptrdiff_t;

 public:
  CoefficientSubView(Int l, Int mMax, View view, MRange)
      : _l{l}, _mMax{std::min(l, mMax)}, _view{view} {
    assert(_l >= 0);
    assert(_mMax >= 0);
    assert(this->size() == MaxOrder() - MinOrder() + 1);
  }

  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto Degree() const { return _l; }
  auto MaxOrder() const { return std::min(_l, _mMax); }
  auto MinOrder() const {
    if constexpr (std::same_as<MRange, All>) {
      return -std::min(_l, _mMax);
    } else {
      return 0;
    }
  }
  auto Orders() const {
    return std::ranges::views::iota(MinOrder(), MaxOrder() + 1);
  }
  auto NegativeOrders() const {
    return std::ranges::views::iota(MinOrder(), 0);
  }
  auto NonNegativeOrders() const {
    return std::ranges::views::iota(0, MaxOrder() + 1);
  }

  auto Index(Int m) {
    if constexpr (std::same_as<MRange, All>) {
      return m + MinOrder();
    } else {
      return m;
    }
  }

  // Return value at given order.
  auto operator()(Int m) const { return this->operator[](Index(m)); }
  auto& operator()(Int m) { return this->operator[](Index(m)); }

 private:
  Int _l;
  Int _mMax;
  View _view;
};

// Deduction guide to construct from a range.
template <std::ranges::viewable_range Range, OrderIndexRange MRange>
CoefficientSubView(std::ptrdiff_t, std::ptrdiff_t, Range&&, MRange)
    -> CoefficientSubView<std::ranges::views::all_t<Range>, MRange>;

//------------------------------------------------------------//
//            View to data at a single upper index            //
//------------------------------------------------------------//

template <std::ranges::view View, OrderIndexRange MRange>
class CoefficientView
    : public std::ranges::view_interface<CoefficientView<View, MRange>> {
  using Int = std::ptrdiff_t;

 public:
  CoefficientView(Int lMax, Int mMax, Int n, View view, MRange)
      : _lMax{lMax}, _mMax{std::min(lMax, mMax)}, _n{n}, _view{view} {
    assert(_lMax >= 0);
    assert(mMax >= 0);
    assert(std::abs(_n) <= _lMax);
    assert(_view.size() == CoefficientIndices<MRange>(_lMax, _mMax, n).size());
  }

  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto UpperIndex() const { return _n; }

  auto MinDegree() const { return std::abs(_n); }
  auto MaxDegree() const { return _lMax; }
  auto Degrees() const {
    return std::ranges::views::iota(MinDegree(), MaxDegree() + 1);
  }

  auto Index(Int l, Int m) const {
    GSHIndices<MRange>(_lMax, _mMax, _n)(l, m);
  };

  auto OffsetForDegree(Int l) const {
    return GSHIndices<MRange>(l, _mMax, _n).size() - 1;
  }

  auto SizeForDegree(Int l) const {
    if constexpr (std::same_as<MRange, All>) {
      return 2 * std::min(l, _mMax) + 1;
    } else {
      return std::min(l, _mMax) + 1;
    }
  }

  auto operator()(Int l) const {
    auto start = std::next(begin(), OffsetForDegree(l));
    auto finish = std::next(start, SizeForDegree(l));
    auto view =
        std::ranges::subrange(start, finish) | std::ranges::views::as_const;
    return CoefficientSubView(l, _mMax, view, MRange{});
  }

  auto operator()(Int l) {
    auto start = std::next(begin(), OffsetForDegree(l));
    auto finish = std::next(start, SizeForDegree(l));
    auto view = std::ranges::subrange(start, finish);
    return CoefficientSubView(l, _mMax, view, MRange{});
  }

 private:
  Int _lMax;
  Int _mMax;
  Int _n;
  View _view;
};

// Deduction guide to construct from a range.
template <std::ranges::viewable_range Range, OrderIndexRange MRange>
CoefficientView(std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t, Range&&, MRange)
    -> CoefficientView<std::ranges::views::all_t<Range>, MRange>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_NEW_GUARD_H