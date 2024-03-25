#ifndef GSH_TRANS_COMPONENTS_GUARD_HVector
#define GSH_TRANS_COMPONENTS_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <ranges>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

//----------------------------------------------------------//
//                   Component field view                   //
//----------------------------------------------------------//
template <typename _Grid, std::ranges::view _View>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class ComponentField
    : public std::ranges::view_interface<ComponentField<_Grid, _View>> {
  using Int = std::ptrdiff_t;

 public:
  // Public type aliases.
  using Grid = _Grid;
  using View = _View;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar = std::ranges::range_value_t<_View>;

  // Constructors.
  ComponentField() = default;

  ComponentField(_Grid grid, Int n, _View view)
      : _grid{grid}, _n{n}, _view{view} {
    assert(std::ranges::contains(_grid.UpperIndices(), UpperIndex()));
    assert(_view.size() == _grid.ComponentSize());
  }

  ComponentField(const ComponentField&) = default;
  ComponentField(ComponentField&&) = default;

  // Assignment.
  ComponentField& operator=(const ComponentField&) = default;
  ComponentField& operator=(ComponentField&&) = default;

  ComponentField& operator=(Scalar s)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    for (auto& val : _view) val = s;
    return *this;
  }

  template <typename Function>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::regular_invocable<Function, Real, Real>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Real>,
                                 Scalar>;
  }
  ComponentField& operator=(Function f) {
    auto iter = begin();
    for (auto [theta, phi] : _grid.Points()) {
      *iter++ = f(theta, phi);
    }
    return *this;
  }

  template <typename OtherGrid, std::ranges::view OtherView>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::convertible_to<Scalar, std::ranges::range_value_t<OtherView>>;
  }
  ComponentField& operator=(ComponentField<OtherGrid, OtherView> other) {
    assert(std::ranges::equal(other.CoLatitudes(), CoLatitudes()));
    assert(std::ranges::equal(other.Longitudes(), Longitudes()));
    std::ranges::copy(other, begin());
    _n = other.UpperIndex();
    return *this;
  }

  // Methods related to the grid.
  auto GetGrid() const { return _grid; }
  auto NumberOfLongitudes() const { return _grid.NumberOfLongitudes(); }
  auto NumberOfCoLatitudes() const { return _grid.NumberOfCoLatitudes(); }
  auto Longitudes() const { return _grid.Longitudes(); }
  auto CoLatitudes() const { return _grid.CoLatitudes(); }

  // Data access methods.
  auto UpperIndex() const { return _n; }
  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto Index(Int iTheta, Int iPhi) const {
    return iTheta * NumberOfLongitudes() + iPhi;
  }

  auto operator()(Int iTheta, Int iPhi) const {
    return this->operator[](Index(iTheta, iPhi));
  }

 private:
  _Grid _grid;
  _View _view;
  Int _n;
};

template <typename Grid, std::ranges::range R>
ComponentField(Grid, std::ptrdiff_t, R&&)
    -> ComponentField<Grid, std::ranges::views::all_t<R>>;

//-------------------------------------------------------//
//         Range adaptor to form Component fields        //
//-------------------------------------------------------//
template <typename _Grid>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormComponentField
    : public std::ranges::range_adaptor_closure<FormComponentField<_Grid>> {
  using Int = std::ptrdiff_t;

 public:
  FormComponentField(_Grid grid, Int n) : _grid{grid}, _n{n} {}

  template <std::ranges::view View>
  auto operator()(View view) {
    return ComponentField(_grid, _n, view);
  }

 private:
  _Grid _grid;
  Int _n;
};

//-----------------------------------------------------//
//          Operations on component fields             //
//-----------------------------------------------------//

template <typename Grid, std::ranges::view View, typename Function>
auto ComponentFieldUnary(ComponentField<Grid, View> u, Function f) {
  return u | std::ranges::views::transform(f) |
         FormComponentField(u.GetGrid(), u.UpperIndex());
}

template <typename Grid, std::ranges::view View>
auto operator-(ComponentField<Grid, View> u) {
  return ComponentFieldUnary(u, [](auto x) { return -x; });
}

template <typename Grid, std::ranges::view View>
auto real(ComponentField<Grid, View> u) {
  return ComponentFieldUnary(u, [](auto x) { return std::real(x); });
}

template <typename Grid, std::ranges::view View>
auto imag(ComponentField<Grid, View> u) {
  return ComponentFieldUnary(u, [](auto x) { return std::imag(x); });
}

template <typename Grid, std::ranges::view View>
auto abs(ComponentField<Grid, View> u) {
  return ComponentFieldUnary(u, [](auto x) { return std::abs(x); });
}

template <typename Grid, std::ranges::view View>
auto conj(ComponentField<Grid, View> u) {
  if constexpr (RealFloatingPoint<std::ranges::range_value_t<View>>) {
    return u | FormComponentField(u.GetGrid(), -u.UpperIndex());
  } else {
    return u |
           std::ranges::views::transform([](auto x) { return std::conj(x); }) |
           FormComponentField(u.GetGrid(), -u.UpperIndex());
  }
}

template <typename Grid, std::ranges::view View, typename Function, typename S>
requires requires() {
  requires std::integral<S> || RealOrComplexFloatingPoint<S>;
}
auto ComponentFieldUnaryWithScalar(ComponentField<Grid, View> u, Function f,
                                   S s) {
  using U = std::ranges::range_value_t<View>;
  using T =
      std::conditional_t<std::integral<S>, U,
                         std::conditional_t<ComplexFloatingPoint<U>, U, S>>;
  auto t = T(s);
  return u | std::ranges::views::transform([t, f](auto x) { return f(x, t); }) |
         FormComponentField(u.GetGrid(), u.UpperIndex());
}

template <typename Grid, std::ranges::view View, typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
auto operator*(ComponentField<Grid, View> u, S s) {
  return ComponentFieldUnaryWithScalar(u, std::multiplies<>(), s);
}

template <typename Grid, std::ranges::view View, typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
auto operator*(S s, ComponentField<Grid, View> u) {
  return u * s;
}

template <typename Grid, std::ranges::view View, typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
auto operator+(ComponentField<Grid, View> u, S s) {
  return ComponentFieldUnaryWithScalar(u, std::plus<>(), s);
}

template <typename Grid, std::ranges::view View, typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
auto operator+(S s, ComponentField<Grid, View> u) {
  return u + s;
}

template <typename Grid, std::ranges::view View, typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
auto operator-(ComponentField<Grid, View> u, S s) {
  return ComponentFieldUnaryWithScalar(u, std::minus<>(), s);
}

template <typename Grid, std::ranges::view View, typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
auto operator/(ComponentField<Grid, View> u, S s) {
  return ComponentFieldUnaryWithScalar(u, std::divides<>(), s);
}

template <typename Grid, std::ranges::view View, typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
auto pow(ComponentField<Grid, View> u, S s) {
  return ComponentFieldUnaryWithScalar(
      u, [](auto x, auto y) { return std::pow(x, y); }, s);
}

template <typename Grid, std::ranges::view View1, std::ranges::view View2>
auto operator+(ComponentField<Grid, View1> u1, ComponentField<Grid, View2> u2) {
  assert(u1.UpperIndex() == u2.UpperIndex());
  return std::ranges::views::zip_transform(std::plus<>(), u1, u2) |
         FormComponentField(u1.GetGrid(), u1.UpperIndex());
}

template <typename Grid, std::ranges::view View1, std::ranges::view View2>
auto operator-(ComponentField<Grid, View1> u1, ComponentField<Grid, View2> u2) {
  assert(u1.UpperIndex() == u2.UpperIndex());
  return std::ranges::views::zip_transform(std::minus<>(), u1, u2) |
         FormComponentField(u1.GetGrid(), u1.UpperIndex());
}

template <typename Grid, std::ranges::view View1, std::ranges::view View2>
auto operator*(ComponentField<Grid, View1> u1, ComponentField<Grid, View2> u2) {
  auto indices = std::vector<std::ptrdiff_t>{};
  return std::ranges::views::zip_transform(std::multiplies<>(), u1, u2) |
         FormComponentField(u1.GetGrid(), u1.UpperIndex() + u2.UpperIndex());
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_COMPONENTS_GUARD_H
