#ifndef GSH_TRANS_FIELD_GUARD_H
#define GSH_TRANS_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <ranges>
#include <utility>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

template <std::ranges::view _View, typename _Grid, RealOrComplexValued _Value>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
  requires(std::same_as<_Value, ComplexValued> &&
           ComplexFloatingPoint<std::ranges::range_value_t<_View>>) ||
              (std::same_as<_Value, RealValued> &&
               RealFloatingPoint<std::ranges::range_value_t<_View>>);
  requires(std::same_as<_Value, ComplexValued> &&
           std::same_as<typename _Grid::MRange, All>) ||
              std::same_as<_Value, RealValued>;
}
class Field : public std::ranges::view_interface<Field<_View, _Grid, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  // Public type aliases.
  using Grid = _Grid;
  using View = _View;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Constructors.
  Field() = default;

  Field(_View view, _Grid grid) : _view{view}, _grid{grid} {
    assert(_view.size() == _grid.FieldSize());
  }

  Field(const Field&) = default;
  Field(Field&&) = default;

  // Assignment.
  Field& operator=(const Field&) = default;
  Field& operator=(Field&&) = default;

  Field& operator=(Scalar s)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()), begin());
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
  Field& operator=(Function f) {
    std::ranges::copy(_grid.InterpolateFunction(f), begin());
    return *this;
  }

  template <typename OtherGrid, std::ranges::view OtherView,
            RealOrComplexValued OtherValue>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::convertible_to<std::ranges::range_value_t<OtherView>, Scalar>;
  }
  Field& operator=(Field<OtherGrid, OtherView, OtherValue> other) {
    assert(this->size() == other.size());
    std::ranges::copy(other, begin());
    return *this;
  }

  // Methods related to the grid.
  auto GetGrid() const { return _grid; }
  auto MaxDegree() const { return _grid.MaxDegree(); }
  auto NumberOfLongitudes() const { return _grid.NumberOfLongitudes(); }
  auto NumberOfCoLatitudes() const { return _grid.NumberOfCoLatitudes(); }
  auto Longitudes() const { return _grid.Longitudes(); }
  auto CoLatitudes() const { return _grid.CoLatitudes(); }

  // Data access methods.
  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto Index(Int iTheta, Int iPhi) const {
    return iTheta * NumberOfLongitudes() + iPhi;
  }

  auto operator()(Int iTheta, Int iPhi) const {
    return this->operator[](Index(iTheta, iPhi));
  }

  auto& operator()(Int iTheta, Int iPhi)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    return this->operator[](Index(iTheta, iPhi));
  }

 private:
  _View _view;
  _Grid _grid;
};

// Deduction guides for construction from ranges.
template <std::ranges::viewable_range R, typename Grid,
          RealOrComplexValued Value>
Field(R&&, Grid) -> Field<std::ranges::views::all_t<R>, Grid, Value>;

// Type aliases for real and complex fields.
template <std::ranges::view View, typename Grid>
using RealField = Field<View, Grid, RealValued>;

template <std::ranges::view View, typename Grid>
using ComplexField = Field<View, Grid, ComplexValued>;

//-------------------------------------------------------//
//                Range adaptors for Fields              //
//-------------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormField
    : public std::ranges::range_adaptor_closure<FormField<_Grid, _Value>> {
 public:
  FormField(_Grid grid) : _grid{grid} {}

  template <std::ranges::view _View>
  auto operator()(_View view) {
    return Field<_View, _Grid, _Value>(view, _grid);
  }

 private:
  _Grid _grid;
};

//---------------------------------------------------//
//                  Field -> Field                   //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(Field<View, Grid, Value> u) {
  return u | std::ranges::views::transform([](auto x) { return -x; }) |
         FormField<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto real(Field<View, Grid, Value> u) {
  return u |
         std::ranges::views::transform([](auto x) { return std::real(x); }) |
         FormField<Grid, RealValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto imag(Field<View, Grid, Value> u) {
  return u |
         std::ranges::views::transform([](auto x) { return std::imag(x); }) |
         FormField<Grid, RealValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto abs(Field<View, Grid, Value> u) {
  return u | std::ranges::views::transform([](auto x) { return std::abs(x); }) |
         FormField<Grid, RealValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto conj(Field<View, Grid, Value> u) {
  if constexpr (std::same_as<Value, RealValued>) {
    return u;
  } else {
    return u |
           std::ranges::views::transform([](auto x) { return std::conj(x); }) |
           FormField<Grid, Value>(u.GetGrid());
  }
}

//---------------------------------------------------//
//               Field x Scalar -> Field             //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(Field<View, Grid, Value> u, std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         FormField<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(std::ranges::range_value_t<View> s, Field<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(Field<View, Grid, Value> u, std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         FormField<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(std::ranges::range_value_t<View> s, Field<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(Field<View, Grid, Value> u, std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         FormField<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator/(Field<View, Grid, Value> u, std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         FormField<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto pow(Field<View, Grid, Value> u, std::ranges::range_value_t<View> s) {
  return u |
         std::ranges::views::transform([s](auto x) { return std::pow(x, s); }) |
         FormField<Grid, Value>(u.GetGrid());
}

//---------------------------------------------------//
//               Field x Complex -> Field            //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator*(Field<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         FormField<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator*(typename Grid::Complex s, Field<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator+(Field<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         FormField<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator+(typename Grid::Complex s, Field<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator-(Field<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         FormField<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator/(Field<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         FormField<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto pow(Field<View, Grid, Value> u, typename Grid::Complex s) {
  return u |
         std::ranges::views::transform([s](auto x) { return std::pow(x, s); }) |
         FormField<Grid, ComplexValued>(u.GetGrid());
}

//---------------------------------------------------//
//                Field x Field -> Field             //
//---------------------------------------------------//
template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator+(Field<View1, Grid, Value1> u1, Field<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x + y; },
                                           u1, u2) |
         FormField<Grid, Value>(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator-(Field<View1, Grid, Value1> u1, Field<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x - y; },
                                           u1, u2) |
         FormField<Grid, Value>(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator*(Field<View1, Grid, Value1> u1, Field<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x * y; },
                                           u1, u2) |
         FormField<Grid, Value>(u1.GetGrid());
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELD_GUARD_H