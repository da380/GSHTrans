#ifndef GSH_TRANS_TENSOR_INDICES_GUARD_H
#define GSH_TRANS_TENSOR_INDICES_GUARD_H

#include <array>
#include <map>
#include <ranges>
#include <tuple>

namespace GSHTrans {

template <std::size_t N>
class CartesianPower
    : public std::ranges::range_adaptor_closure<CartesianPower<N>> {
 public:
  CartesianPower() = default;

  template <std::ranges::view View>
  auto operator()(View v) const {
    auto repeat = std::array<View, N>();
    repeat.fill(v);
    return std::apply(std::ranges::views::cartesian_product, repeat);
  }
};

template <std::size_t Rank>
class TensorIndex : public std::ranges::view_interface<TensorIndex<Rank>> {
  using Int = std::ptrdiff_t;

 public:
  TensorIndex() = default;
  TensorIndex(std::array<Int, Rank>& index) : _index{index} {}
  TensorIndex(std::array<Int, Rank>&& index) : TensorIndex(index) {}

  template <typename... I>
  requires(std::convertible_to<I, Int> && ...) && (sizeof...(I) == Rank)
  TensorIndex(I... alpha) : TensorIndex(std::array<Int, Rank>{alpha...}) {}

  auto& operator-() {
    for (auto& val : _index) val *= -1;
    return *this;
  }

  TensorIndex(const TensorIndex&) = default;
  TensorIndex(TensorIndex&&) = default;

  TensorIndex& operator=(const TensorIndex&) = default;
  TensorIndex& operator=(TensorIndex&&) = default;

  auto operator<=>(const TensorIndex&) const = default;

  auto begin() { return _index.begin(); }
  auto end() { return _index.end(); }

 private:
  std::array<Int, Rank> _index;
};

template <std::size_t Rank>
auto abs(TensorIndex<Rank>& index) {
  return std::ranges::fold_left_first(index, std::plus<>()).value_or(0);
}

template <std::size_t Rank>
auto abs(TensorIndex<Rank>&& index) {
  return abs(index);
}

template <std::size_t Rank>
auto TensorIndices() {
  return std::ranges::views::iota(-1, 2) | CartesianPower<Rank>() |
         std::ranges::views::transform([](auto tuple) {
           return std::apply(
               [](auto... alpha) { return TensorIndex<Rank>(alpha...); },
               tuple);
         });
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_INDICES_GUARD_H