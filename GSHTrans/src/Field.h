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
class FieldView
    : public std::ranges::view_interface<FieldView<_View, _Grid, _Value>> {
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
  FieldView() = default;

  FieldView(_View view, _Grid grid) : _view{view}, _grid{grid} {
    assert(_view.size() == _grid.FieldSize());
  }

  FieldView(const FieldView&) = default;
  FieldView(FieldView&&) = default;

  // Assignment.
  FieldView& operator=(const FieldView&) = default;
  FieldView& operator=(FieldView&&) = default;

  FieldView& operator=(Scalar s)
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
  FieldView& operator=(Function f) {
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
  FieldView& operator=(FieldView<OtherGrid, OtherView, OtherValue> other) {
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

//-------------------------------------------------------//
//                Range adaptors for Fields              //
//-------------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class Field : public std::ranges::range_adaptor_closure<Field<_Grid, _Value>> {
 public:
  Field(_Grid grid) : _grid{grid} {}

  template <std::ranges::view _View>
  auto operator()(_View view) {
    return FieldView<_View, _Grid, _Value>(view, _grid);
  }

 private:
  _Grid _grid;
};

// Define type aliases for real and complex adaptors.
template <typename _Grid>
using RealField = Field<_Grid, RealValued>;

template <typename _Grid>
using ComplexField = Field<_Grid, ComplexValued>;

//---------------------------------------------------//
//                  Field -> Field                   //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(FieldView<View, Grid, Value> u) {
  return u | std::ranges::views::transform([](auto x) { return -x; }) |
         Field<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto real(FieldView<View, Grid, Value> u) {
  return u |
         std::ranges::views::transform([](auto x) { return std::real(x); }) |
         Field<Grid, RealValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto imag(FieldView<View, Grid, Value> u) {
  return u |
         std::ranges::views::transform([](auto x) { return std::imag(x); }) |
         Field<Grid, RealValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto abs(FieldView<View, Grid, Value> u) {
  return u | std::ranges::views::transform([](auto x) { return std::abs(x); }) |
         Field<Grid, RealValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto conj(FieldView<View, Grid, Value> u) {
  if constexpr (std::same_as<Value, RealValued>) {
    return u;
  } else {
    return u |
           std::ranges::views::transform([](auto x) { return std::conj(x); }) |
           Field<Grid, Value>(u.GetGrid());
  }
}

//---------------------------------------------------//
//               Field x Scalar -> Field             //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(FieldView<View, Grid, Value> u,
               std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         Field<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(std::ranges::range_value_t<View> s,
               FieldView<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(FieldView<View, Grid, Value> u,
               std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         Field<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(std::ranges::range_value_t<View> s,
               FieldView<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(FieldView<View, Grid, Value> u,
               std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         Field<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator/(FieldView<View, Grid, Value> u,
               std::ranges::range_value_t<View> s) {
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         Field<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto pow(FieldView<View, Grid, Value> u, std::ranges::range_value_t<View> s) {
  return u |
         std::ranges::views::transform([s](auto x) { return std::pow(x, s); }) |
         Field<Grid, Value>(u.GetGrid());
}

//---------------------------------------------------//
//               Field x Complex -> Field            //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator*(FieldView<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         Field<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator*(typename Grid::Complex s, FieldView<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator+(FieldView<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         Field<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator+(typename Grid::Complex s, FieldView<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator-(FieldView<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         Field<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto operator/(FieldView<View, Grid, Value> u, typename Grid::Complex s) {
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         Field<Grid, ComplexValued>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<typename Grid::MRange, All>
auto pow(FieldView<View, Grid, Value> u, typename Grid::Complex s) {
  return u |
         std::ranges::views::transform([s](auto x) { return std::pow(x, s); }) |
         Field<Grid, ComplexValued>(u.GetGrid());
}

//---------------------------------------------------//
//                Field x Field -> Field             //
//---------------------------------------------------//
template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator+(FieldView<View1, Grid, Value1> u1,
               FieldView<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x + y; },
                                           u1, u2) |
         Field<Grid, Value>(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator-(FieldView<View1, Grid, Value1> u1,
               FieldView<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x - y; },
                                           u1, u2) |
         Field<Grid, Value>(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator*(FieldView<View1, Grid, Value1> u1,
               FieldView<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x * y; },
                                           u1, u2) |
         Field<Grid, Value>(u1.GetGrid());
}

//------------------------------------------------//
//                Helper functions                //
//------------------------------------------------//
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
auto AllocateRealField(Grid grid) {
  auto data =
      std::make_unique<FFTWpp::vector<typename Grid::Real>>(grid.FieldSize());
  return std::pair(std::ranges::views::all(*data) | RealField(grid),
                   std::move(data));
}

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
auto AllocateComplexField(Grid grid) {
  auto data = std::make_unique<FFTWpp::vector<typename Grid::Complex>>(
      grid.FieldSize());
  return std::pair(std::ranges::views::all(*data) | ComplexField(grid),
                   std::move(data));
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELD_GUARD_H