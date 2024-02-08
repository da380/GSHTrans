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

template <std::ranges::view View, OrderIndexRange MRange>
class OrdersView
    : public std::ranges::view_interface<OrdersView<View, MRange>> {
  using Int = std::ptrdiff_t;

 public:
  OrdersView(Int l, Int mMax, View view, MRange)
      : _l{l}, _mMax{std::min(l, mMax)}, _view{view} {
    assert(_l >= 0);
    assert(_mMax >= 0);
    if constexpr (std::same_as<MRange, All>) {
      assert(this->size() == 2 * _mMax + 1);
    } else {
      assert(this->size() == _mMax + 1);
    }
  }

  // Set begin and end from which other functions are
  // deduced via CRTP.
  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  // Return information.
  auto Degree() const { return _l; }
  auto MaxOrder() const { return _mMax; }
  auto MinOrder() const {
    if constexpr (std::same_as<MRange, All>) {
      return -_mMax;
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

  // Return value at given order.
  auto operator()(Int m) const { return this->operator[](m + MinOrder()); }
  auto& operator()(Int m) { return this->operator[](m + MinOrder()); }

  // Return "enumerated" values.
  auto Enumerate() const {
    return std::ranges::views::zip(Orders(),
                                   _view | std::ranges::views::as_const);
  }
  auto Enumerate() { return std::ranges::views::zip(Orders(), _view); }

  auto EnumerateNegative() const {
    return std::ranges::views::zip(
        NegativeOrders(), _view | std::ranges::views::take(MinOrder()) |
                              std::ranges::views::as_const);
  }

  auto EnumerateNegative() {
    return std::ranges::views::zip(
        NegativeOrders(), _view | std::ranges::views::take(MinOrder()));
  }

 private:
  Int _l;
  Int _mMax;
  View _view;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_NEW_GUARD_H