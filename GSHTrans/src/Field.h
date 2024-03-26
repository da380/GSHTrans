#ifndef GSH_TRANS_FIELD_GUARD_H
#define GSH_TRANS_FIELD_GUARD_H

#include <algorithm>
#include <complex>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

//--------------------------------------------------------------//
//                       Class declarations                     //
//--------------------------------------------------------------//
template <std::ranges::view _View, typename _Grid>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
  requires ComplexFloatingPoint<std::ranges::range_value_t<_View>>;
  requires std::same_as<typename _Grid::MRange, All>;
}
class Field;

template <std::ranges::view _View, typename _Grid>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
  requires RealFloatingPoint<std::ranges::range_value_t<_View>>;
}
class RealField;

//----------------------------------------------------------//
//                      Field definition                    //
//----------------------------------------------------------//
template <std::ranges::view _View, typename _Grid>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
  requires ComplexFloatingPoint<std::ranges::range_value_t<_View>>;
  requires std::same_as<typename _Grid::MRange, All>;
}
class Field : public std::ranges::view_interface<Field<_View, _Grid>> {
  using Int = std::ptrdiff_t;

 public:
  // Public type aliases.
  using Grid = _Grid;
  using View = _View;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // Constructors.
  Field() = default;

  Field(_View view, _Grid grid) : _view{view}, _grid{grid} {
    assert(_view.size() == _grid.ComponentSize());
  }

  Field(const Field&) = default;
  Field(Field&&) = default;

  // Assignment.
  Field& operator=(const Field&) = default;
  Field& operator=(Field&&) = default;

  Field& operator=(Complex s)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    std::ranges::copy(std::ranges::views::repeat(s, _grid.ComponentSize()),
                      begin());
    return *this;
  }

  template <typename Function>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::regular_invocable<Function, Real, Real>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Real>,
                                 Complex>;
  }
  Field& operator=(Function f) {
    std::ranges::copy(_grid.InterpolateFunction(f), begin());
    return *this;
  }

  template <typename OtherGrid, std::ranges::view OtherView>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::convertible_to<std::ranges::range_value_t<OtherView>,
                                 Complex>;
  }
  Field& operator=(Field<OtherGrid, OtherView> other) {
    assert(std::ranges::equal(other.CoLatitudes(), CoLatitudes()));
    assert(std::ranges::equal(other.Longitudes(), Longitudes()));
    std::ranges::copy(other, begin());
    return *this;
  }

  template <typename OtherGrid, std::ranges::view OtherView>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::convertible_to<std::ranges::range_value_t<OtherView>,
                                 Complex>;
  }
  Field& operator=(RealField<OtherGrid, OtherView> other) {
    assert(std::ranges::equal(other.CoLatitudes(), CoLatitudes()));
    assert(std::ranges::equal(other.Longitudes(), Longitudes()));
    std::ranges::copy(other | std::ranges::views::transform(
                                  [](auto x) -> Complex { return x; }),
                      begin());
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

 private:
  _View _view;
  _Grid _grid;
};

// Deduction guides to construct Fields from ranges.
template <std::ranges::viewable_range R, typename Grid>
Field(R&&, Grid) -> Field<std::ranges::views::all_t<R>, Grid>;

//----------------------------------------------------------//
//                      RealField definition                //
//----------------------------------------------------------//
template <std::ranges::view _View, typename _Grid>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
  requires RealFloatingPoint<std::ranges::range_value_t<_View>>;
}
class RealField : public std::ranges::view_interface<RealField<_View, _Grid>> {
  using Int = std::ptrdiff_t;

 public:
  // Public type aliases.
  using Grid = _Grid;
  using View = _View;
  using Real = typename _Grid::Real;

  // Constructors.
  RealField() = default;

  RealField(_View view, _Grid grid) : _view{view}, _grid{grid} {
    assert(_view.size() == _grid.ComponentSize());
  }

  RealField(const RealField&) = default;
  RealField(RealField&&) = default;

  // Assignment.
  RealField& operator=(const RealField&) = default;
  RealField& operator=(RealField&&) = default;

  RealField& operator=(Real s)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    std::ranges::copy(std::ranges::views::repeat(s, _grid.ComponentSize()),
                      begin());
    return *this;
  }

  template <typename Function>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::regular_invocable<Function, Real, Real>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Real>,
                                 Real>;
  }
  RealField& operator=(Function f) {
    std::ranges::copy(_grid.InterpolateFunction(f), begin());
    return *this;
  }

  template <typename OtherGrid, std::ranges::view OtherView>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::convertible_to<std::ranges::range_value_t<OtherView>, Real>;
  }
  RealField& operator=(RealField<OtherGrid, OtherView> other) {
    assert(std::ranges::equal(other.CoLatitudes(), CoLatitudes()));
    assert(std::ranges::equal(other.Longitudes(), Longitudes()));
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

 private:
  _View _view;
  _Grid _grid;
};

// Deduction guides to construct RealFields from ranges.
template <std::ranges::viewable_range R, typename Grid>
RealField(R&&, Grid) -> RealField<std::ranges::views::all_t<R>, Grid>;

//-------------------------------------------------------//
//        Range adaptors for Fields and RealFields       //
//-------------------------------------------------------//
template <typename _Grid>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormField : public std::ranges::range_adaptor_closure<FormField<_Grid>> {
 public:
  FormField(_Grid grid) : _grid{grid} {}

  template <std::ranges::view View>
  auto operator()(View view) {
    return Field(view, _grid);
  }

 private:
  _Grid _grid;
};

template <typename _Grid>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormRealField
    : public std::ranges::range_adaptor_closure<FormRealField<_Grid>> {
 public:
  FormRealField(_Grid grid) : _grid{grid} {}

  template <std::ranges::view View>
  auto operator()(View view) {
    return RealField(view, _grid);
  }

 private:
  _Grid _grid;
};

//---------------------------------------------------//
//                  Field -> Field                   //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, typename Function>
auto FieldUnary(Field<View, Grid> u, Function f) {
  return u | std::ranges::views::transform(f) | FormField(u.GetGrid());
}

template <std::ranges::view View, typename Grid>
auto operator-(Field<View, Grid> u) {
  return FieldUnary(u, [](auto x) { return -x; });
}

template <std::ranges::view View, typename Grid>
auto conj(Field<View, Grid> u) {
  return FieldUnary(u, [](auto x) { return std::conj(x); });
}

//---------------------------------------------------//
//                RealField -> RealField             //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, typename Function>
auto RealFieldUnary(RealField<View, Grid> u, Function f) {
  return u | std::ranges::views::transform(f) | FormRealField(u.GetGrid());
}

template <std::ranges::view View, typename Grid>
auto operator-(RealField<View, Grid> u) {
  return RealFieldUnary(u, [](auto x) { return -x; });
}

template <std::ranges::view View, typename Grid>
auto conj(RealField<View, Grid> u) {
  return u;
}

//---------------------------------------------------//
//                 Field -> RealField                //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid>
auto real(Field<View, Grid> u) {
  return u |
         std::ranges::views::transform([](auto x) { return std::real(x); }) |
         FormRealField(u.GetGrid());
}

template <std::ranges::view View, typename Grid>
auto imag(Field<View, Grid> u) {
  return u |
         std::ranges::views::transform([](auto x) { return std::imag(x); }) |
         FormRealField(u.GetGrid());
}

template <std::ranges::view View, typename Grid>
auto abs(Field<View, Grid> u) {
  return u | std::ranges::views::transform([](auto x) { return std::abs(x); }) |
         FormRealField(u.GetGrid());
}

//---------------------------------------------------//
//               Field x Complex -> Field            //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, typename Function>
auto FieldUnaryWithScalar(Field<View, Grid> u, Function f,
                          typename Grid::Complex s) {
  return u | std::ranges::views::transform([s, f](auto x) { return f(x, s); }) |
         FormField(u.GetGrid());
}

template <std::ranges::view View, typename Grid>
auto operator*(Field<View, Grid> u, typename Grid::Complex s) {
  return FieldUnaryWithScalar(u, std::multiplies<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator*(typename Grid::Complex s, Field<View, Grid> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid>
auto operator+(Field<View, Grid> u, typename Grid::Complex s) {
  return FieldUnaryWithScalar(u, std::plus<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator+(typename Grid::Complex s, Field<View, Grid> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid>
auto operator-(Field<View, Grid> u, typename Grid::Complex s) {
  return FieldUnaryWithScalar(u, std::minus<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator/(Field<View, Grid> u, typename Grid::Complex s) {
  return FieldUnaryWithScalar(u, std::divides<>(), s);
}

template <std::ranges::view View, typename Grid>
auto pow(Field<View, Grid> u, typename Grid::Complex s) {
  return FieldUnaryWithScalar(
      u, [](auto x, auto y) { return std::pow(x, y); }, s);
}

//---------------------------------------------------//
//             RealField x Real -> RealField         //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, typename Function>
auto RealFieldUnaryWithScalar(RealField<View, Grid> u, Function f,
                              typename Grid::Real s) {
  return u | std::ranges::views::transform([s, f](auto x) { return f(x, s); }) |
         FormRealField(u.GetGrid());
}

template <std::ranges::view View, typename Grid>
auto operator*(RealField<View, Grid> u, typename Grid::Real s) {
  return RealFieldUnaryWithScalar(u, std::multiplies<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator*(typename Grid::Real s, RealField<View, Grid> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid>
auto operator+(RealField<View, Grid> u, typename Grid::Real s) {
  return RealFieldUnaryWithScalar(u, std::plus<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator+(typename Grid::Real s, RealField<View, Grid> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid>
auto operator-(RealField<View, Grid> u, typename Grid::Real s) {
  return RealFieldUnaryWithScalar(u, std::minus<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator/(RealField<View, Grid> u, typename Grid::Real s) {
  return RealFieldUnaryWithScalar(u, std::divides<>(), s);
}

template <std::ranges::view View, typename Grid>
auto pow(RealField<View, Grid> u, typename Grid::Real s) {
  return RealFieldUnaryWithScalar(
      u, [](auto x, auto y) { return std::pow(x, y); }, s);
}

//---------------------------------------------------//
//             RealField x Complex -> Field          //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, typename Function>
auto RealFieldUnaryWithScalar(RealField<View, Grid> u, Function f,
                              typename Grid::Complex s) {
  return u | std::ranges::views::transform([s, f](auto x) { return f(x, s); }) |
         FormField(u.GetGrid());
}

template <std::ranges::view View, typename Grid>
auto operator*(RealField<View, Grid> u, typename Grid::Complex s) {
  return RealFieldUnaryWithScalar(u, std::multiplies<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator*(typename Grid::Complex s, RealField<View, Grid> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid>
auto operator+(RealField<View, Grid> u, typename Grid::Complex s) {
  return RealFieldUnaryWithScalar(u, std::plus<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator+(typename Grid::Complex s, RealField<View, Grid> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid>
auto operator-(RealField<View, Grid> u, typename Grid::Complex s) {
  return RealFieldUnaryWithScalar(u, std::minus<>(), s);
}

template <std::ranges::view View, typename Grid>
auto operator/(RealField<View, Grid> u, typename Grid::Complex s) {
  return RealFieldUnaryWithScalar(u, std::divides<>(), s);
}

template <std::ranges::view View, typename Grid>
auto pow(RealField<View, Grid> u, typename Grid::Complex s) {
  return RealFieldUnaryWithScalar(
      u, [](auto x, auto y) { return std::pow(x, y); }, s);
}

//---------------------------------------------------//
//                Field x Field -> Field             //
//---------------------------------------------------//
template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator+(Field<View1, Grid> u1, Field<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(std::plus<>(), u1, u2) |
         FormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator-(Field<View1, Grid> u1, Field<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(std::minus<>(), u1, u2) |
         FormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator*(Field<View1, Grid> u1, Field<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(std::multiplies<>(), u1, u2) |
         FormField(u1.GetGrid());
}

//---------------------------------------------------//
//          RealField x RealField -> RealField       //
//---------------------------------------------------//
template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator+(RealField<View1, Grid> u1, RealField<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(std::plus<>(), u1, u2) |
         RealFormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator-(RealField<View1, Grid> u1, RealField<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(std::minus<>(), u1, u2) |
         RealFormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator*(RealField<View1, Grid> u1, RealField<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(std::multiplies<>(), u1, u2) |
         RealFormField(u1.GetGrid());
}

//---------------------------------------------------//
//              RealField x Field -> Field           //
//---------------------------------------------------//
template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator+(RealField<View1, Grid> u1, Field<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(
             std::plus<>(),
             u1 | std::ranges::views::transform(
                      [](auto x) -> typename Grid::Complex { return x; }),
             u2) |
         FormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator-(RealField<View1, Grid> u1, Field<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(
             std::minus<>(),
             u1 | std::ranges::views::transform(
                      [](auto x) -> typename Grid::Complex { return x; }),
             u2) |
         FormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator*(RealField<View1, Grid> u1, Field<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(
             std::multiplies<>(),
             u1 | std::ranges::views::transform(
                      [](auto x) -> typename Grid::Complex { return x; }),
             u2) |
         FormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator+(Field<View1, Grid> u1, RealField<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(
             std::plus<>(), u1,
             u2 | std::ranges::views::transform(
                      [](auto x) -> typename Grid::Complex { return x; })) |
         FormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator-(Field<View1, Grid> u1, RealField<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(
             std::minus<>(), u1,
             u2 | std::ranges::views::transform(
                      [](auto x) -> typename Grid::Complex { return x; })) |
         FormField(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator*(Field<View1, Grid> u1, RealField<View2, Grid> u2) {
  assert(u1.MaxDegree() == u2.MaxDegree());
  return std::ranges::views::zip_transform(
             std::multiplies<>(), u1,
             u2 | std::ranges::views::transform(
                      [](auto x) -> typename Grid::Complex { return x; })) |
         FormField(u1.GetGrid());
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELD_GUARD_H