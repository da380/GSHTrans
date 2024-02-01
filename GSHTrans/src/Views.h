#ifndef GSH_TRANS_VIEWS_GUARD_H
#define GSH_TRANS_VIEWS_GUARD_H

#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <ranges>
#include <type_traits>

#include "Concepts.h"

namespace GSHTrans {

namespace Views {

template <std::ranges::view View>
requires Field<std::ranges::range_value_t<View>>
auto operator-(View v) {
  return v | std::ranges::views::transform([](auto x) { return -x; });
}

template <std::ranges::view View, Field S, typename Function>
requires requires(View v, S s, Function f) {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
  { f(v[0], v[0]) } -> std::convertible_to<std::ranges::range_value_t<View>>;
}
auto ScalarOperationOnView(View v, S s, Function f) {
  using T = std::ranges::range_value_t<View>;
  return v | std::ranges::views::transform([s, f](auto x) { return f(x, s); });
}

template <std::ranges::view View, Field S>
requires requires() {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
}
auto operator*(View v, S s) {
  return ScalarOperationOnView(v, s, std::multiplies<>());
}

template <std::ranges::view View, Field S>
requires requires() {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
}
auto operator*(S s, View v) {
  return v * s;
}

template <std::ranges::view View, Field S>
requires requires() {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
}
auto operator+(View v, S s) {
  return ScalarOperationOnView(v, s, std::plus<>());
}

template <std::ranges::view View, Field S>
requires requires() {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
}
auto operator+(S s, View v) {
  return v + s;
}

template <std::ranges::view View, Field S>
requires requires() {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
}
auto operator-(View v, S s) {
  return ScalarOperationOnView(v, s, std::minus<>());
}

template <std::ranges::view View, Field S>
requires requires() {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
}
auto operator-(S s, View v) {
  return s + (-v);
}

template <std::ranges::view View, Field S>
requires requires() {
  requires Field<std::ranges::range_value_t<View>>;
  requires std::constructible_from<std::ranges::range_value_t<View>, S>;
}
auto operator/(View v, S s) {
  return ScalarOperationOnView(v, s, std::divides<>());
}

template <std::ranges::view View1, std::ranges::view View2, typename Function>
requires requires(View1 v1, View2 v2, Function f) {
  requires Field<std::ranges::range_value_t<View1>>;
  requires Field<std::ranges::range_value_t<View2>>;
  typename std::common_type<std::ranges::range_value_t<View1>,
                            std::ranges::range_value_t<View2>>::type;
  {
    f(v1[0], v2[0])
  } -> std::convertible_to<
      typename std::common_type<std::ranges::range_value_t<View1>,
                                std::ranges::range_value_t<View2>>::type>;
}
auto BinaryOperation(View1 v1, View2 v2, Function f) {
  return std::ranges::views::zip_transform(f, v1, v2);
}

template <std::ranges::view View1, std::ranges::view View2>
requires requires() {
  requires Field<std::ranges::range_value_t<View1>>;
  requires Field<std::ranges::range_value_t<View2>>;
  typename std::common_type<std::ranges::range_value_t<View1>,
                            std::ranges::range_value_t<View2>>::type;
}
auto operator+(View1 v1, View2 v2) {
  return BinaryOperation(v1, v2, std::plus<>());
}

template <std::ranges::view View1, std::ranges::view View2>
requires requires() {
  requires Field<std::ranges::range_value_t<View1>>;
  requires Field<std::ranges::range_value_t<View2>>;
  typename std::common_type<std::ranges::range_value_t<View1>,
                            std::ranges::range_value_t<View2>>::type;
}
auto operator-(View1 v1, View2 v2) {
  return BinaryOperation(v1, v2, std::minus<>());
}

template <std::ranges::view View1, std::ranges::view View2>
requires requires() {
  requires Field<std::ranges::range_value_t<View1>>;
  requires Field<std::ranges::range_value_t<View2>>;
  typename std::common_type<std::ranges::range_value_t<View1>,
                            std::ranges::range_value_t<View2>>::type;
}
auto operator*(View1 v1, View2 v2) {
  return BinaryOperation(v1, v2, std::multiplies<>());
}

}  // namespace Views

}  // namespace GSHTrans

#endif  // GSH_TRANS_VIEWS_GUARD_H