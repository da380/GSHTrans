#ifndef GSH_TRANS_COEFFICIENT_GUARD_H
#define GSH_TRANS_COEFFICIENT_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <ranges>
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
  requires ComplexFloatingPoint<std::ranges::range_value_t<_View>>;
}
class CoefficientView : public std::ranges::view_interface<
                            CoefficientView<_View, _Grid, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  // Public type aliases.
  using Grid = _Grid;
  using View = _View;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // Constructors.
  CoefficientView() = default;

  CoefficientView(_View view, _Grid grid)
      : _view{view},
        _grid{grid},
        _indices{
            GSHIndices<typename _Grid::MRange>(MaxDegree(), MaxDegree(), 0)} {
    if constexpr (std::same_as<_Value, RealValued>) {
      assert(_view.size() == _grid.RealCoefficientSize(0));
    } else {
      assert(_view.size() == _grid.ComplexCoefficientSize(0));
    }
  }

  CoefficientView(const CoefficientView&) = default;
  CoefficientView(CoefficientView&&) = default;

  // Assignment.
  CoefficientView& operator=(Complex s)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()), begin());
    return *this;
  }

  CoefficientView& operator=(const CoefficientView&) = default;
  CoefficientView& operator=(CoefficientView&&) = default;

  template <typename Function>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::regular_invocable<Function, Int, Int>;
    requires std::convertible_to<std::invoke_result_t<Function, Int, Int>,
                                 Complex>;
  }
  CoefficientView& operator=(Function f) {
    std::ranges::copy(Indices() | std::ranges::views::transform(f), begin());
    return *this;
  }

  template <typename OtherGrid, std::ranges::view OtherView,
            RealOrComplexValued OtherValue>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::convertible_to<std::ranges::range_value_t<OtherView>,
                                 Complex>;
  }
  CoefficientView& operator=(
      CoefficientView<OtherGrid, OtherView, OtherValue> other) {
    assert(this->size() == other.size());
    std::ranges::copy(other, begin());
    return *this;
  }

  // Methods related to the grid.
  auto GetGrid() const { return _grid; }
  auto MaxDegree() const { return _grid.MaxDegree(); }
  auto Degrees() const { return _indices.Degrees(); }
  auto Indices() const { return _indices.Indices(); }

  // Data access methods.
  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto Index(Int l) const { return _indices.Index(l); }
  auto Index(Int l, Int m) const { return _indices.Index(l, m); }

  auto operator()(Int l, Int m) const { return this->operator[](Index(l, m)); }

  auto& operator()(Int l, Int m)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    return this->operator[](Index(l, m));
  }

 private:
  _View _view;
  _Grid _grid;
  GSHIndices<typename _Grid::MRange> _indices;
};

//-------------------------------------------------------//
//             Range adaptors for Coefficient            //
//-------------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class Coefficient
    : public std::ranges::range_adaptor_closure<Coefficient<_Grid, _Value>> {
 public:
  Coefficient(_Grid grid) : _grid{grid} {}

  template <std::ranges::view _View>
  auto operator()(_View view) {
    return CoefficientView<_View, _Grid, _Value>(view, _grid);
  }

 private:
  _Grid _grid;
};

// Define type aliases for real and complex adaptors.
template <typename _Grid>
using RealCoefficient = Coefficient<_Grid, RealValued>;

template <typename _Grid>
using ComplexCoefficient = Coefficient<_Grid, ComplexValued>;

//---------------------------------------------------//
//             Coefficient -> Coefficient            //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(CoefficientView<View, Grid, Value> u) {
  return u | std::ranges::views::transform([](auto x) { return -x; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

//---------------------------------------------------//
//        Coefficient x Complex -> Coefficient       //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(CoefficientView<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<Value, ComplexValued>
auto operator*(typename Grid::Complex s, CoefficientView<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(CoefficientView<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<Value, ComplexValued>
auto operator+(typename Grid::Complex s, CoefficientView<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(CoefficientView<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator/(CoefficientView<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

//---------------------------------------------------//
//          Coefficient x Real -> Coefficient        //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(CoefficientView<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(typename Grid::Real s, CoefficientView<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(CoefficientView<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(typename Grid::Real s, CoefficientView<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(CoefficientView<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator/(CoefficientView<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         Coefficient<Grid, Value>(u.GetGrid());
}

//---------------------------------------------------//
//      Coefficient x Coefficient -> Coefficient     //
//---------------------------------------------------//
template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator+(CoefficientView<View1, Grid, Value1> u1,
               CoefficientView<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x + y; },
                                           u1, u2) |
         Coefficient<Grid, Value>(u1.GetGrid());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value1, RealOrComplexValued Value2>
auto operator-(CoefficientView<View1, Grid, Value1> u1,
               CoefficientView<View2, Grid, Value2> u2) {
  using Value = std::conditional_t<std::same_as<Value1, ComplexValued> ||
                                       std::same_as<Value2, ComplexValued>,
                                   ComplexValued, RealValued>;
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x - y; },
                                           u1, u2) |
         Coefficient<Grid, Value>(u1.GetGrid());
}

//------------------------------------------------//
//                Helper functions                //
//------------------------------------------------//
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
auto AllocateRealCoefficient(Grid grid) {
  auto data = std::make_unique<FFTWpp::vector<typename Grid::Real>>(
      grid.RealCoefficientSize());
  return std::pair(std::ranges::views::all(*data) | RealCoefficient(grid),
                   std::move(data));
}

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
auto AllocateComplexCoefficient(Grid grid) {
  auto data = std::make_unique<FFTWpp::vector<typename Grid::Complex>>(
      grid.ComplexCoefficientSize());
  return std::pair(std::ranges::views::all(*data) | ComplexCoefficient(grid),
                   std::move(data));
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_COEFFICIENT_GUARD_H