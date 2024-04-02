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

template <std::ranges::view _View, RealOrComplexValued _Value>
requires requires() {
  requires ComplexFloatingPoint<std::ranges::range_value_t<_View>>;
}
class SubCoefficient {
  using Int = std::ptrdiff_t;

 public:
  using MRange =
      std::conditional_t<std::same_as<_Value, RealValued>, NonNegative, All>;
  using _GSHSubIndices = GSHSubIndices<MRange>;

  SubCoefficient() = default;

  SubCoefficient(_View view, Int l)
      : _view{view}, _indices{_GSHSubIndices(l, l)} {
    assert(_view.size() == _indices.size());
  }

  SubCoefficient(_View view, _GSHSubIndices indices)
      : _view{view}, _indices{indices} {
    assert(_view.size() == _indices.size());
  }

  auto Degree() const { return _indices.Degree(); }
  auto MinOrder() const { return _indices.MinOrder(); }
  auto MaxOrder() const { return _indices.MaxOrder(); }

  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto Index(Int m) const { return _indices.Index(m); }

  auto operator()(Int m) const { return this->operator[](Index(m)); }
  auto& operator()(Int m)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    return this->operator[](Index(m));
  }

 private:
  _View _view;
  _GSHSubIndices _indices;
};

template <std::ranges::view _View, typename _Grid, RealOrComplexValued _Value>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
  requires ComplexFloatingPoint<std::ranges::range_value_t<_View>>;
  requires(std::same_as<_Value, ComplexValued> &&
           std::same_as<typename _Grid::MRange, All> &&
           std::same_as<typename _Grid::NRange, All>) ||
              std::same_as<_Value, RealValued>;
}
class Coefficient
    : public std::ranges::view_interface<Coefficient<_View, _Grid, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  // Public type aliases.
  using Grid = _Grid;
  using View = _View;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using MRange =
      std::conditional_t<std::same_as<_Value, RealValued>, NonNegative, All>;
  using NRange = MRange;
  using _GSHIndices = GSHIndices<MRange>;

  // Constructors.
  Coefficient() = default;

  Coefficient(_View view, _Grid grid, Int n)
      : _view{view},
        _grid{grid},
        _indices{
            _GSHIndices(_grid.MaxDegree(), _grid.MaxDegree(), std::abs(n))} {
    assert(_view.size() == _indices.size());
  }

  Coefficient(_View view, _Grid grid) : Coefficient(view, grid, 0) {}

  Coefficient(const Coefficient&) = default;
  Coefficient(Coefficient&&) = default;

  // Assignment.
  Coefficient& operator=(Complex s)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()), begin());
    return *this;
  }

  Coefficient& operator=(const Coefficient&) = default;
  Coefficient& operator=(Coefficient&&) = default;

  template <typename Function>
  requires requires() {
    requires std::ranges::output_range<_View,
                                       std::ranges::range_value_t<_View>>;
    requires std::regular_invocable<Function, Int, Int>;
    requires std::convertible_to<std::invoke_result_t<Function, Int, Int>,
                                 Complex>;
  }
  Coefficient& operator=(Function f) {
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
  Coefficient& operator=(Coefficient<OtherGrid, OtherView, OtherValue> other) {
    assert(this->size() == other.size());
    std::ranges::copy(other, begin());
    return *this;
  }

  // Methods related to the grid.
  auto GetGrid() const { return _grid; }
  auto Indices() const { return _indices.Indices(); }
  auto AbsUpperIndex() const { return _indices.UpperIndex(); }
  auto MaxDegree() const { return _indices.MaxDegree(); }
  auto Degrees() const { return _indices.Degrees(); }

  // Data access methods.
  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto Index(Int l) const { return _indices.Index(l); }
  auto Index(Int l, Int m) const { return _indices.Index(l, m); }

  auto operator()(Int l) const {
    auto [offset, indices] = Index(l);
    auto start = std::next(begin(), offset);
    auto finish = std::next(start, indices.size());
    auto view = std::ranges::subrange(start, finish);
    return SubCoefficient<decltype(view), _Value>(view, indices);
  }

  auto operator()(Int l)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    auto [offset, indices] = Index(l);
    auto start = std::next(begin(), offset);
    auto finish = std::next(start, indices.size());
    auto view = std::ranges::subrange(start, finish);
    return SubCoefficient<decltype(view), _Value>(view, indices);
  }

  auto operator()(Int l, Int m) const { return this->operator[](Index(l, m)); }

  auto& operator()(Int l, Int m)
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>
  {
    return this->operator[](Index(l, m));
  }

 private:
  _View _view;
  _Grid _grid;
  _GSHIndices _indices;
};

// Deduction guides for construction from ranges.
template <std::ranges::viewable_range R, typename Grid,
          RealOrComplexValued Value, std::integral I>
Coefficient(R&&, Grid, I)
    -> Coefficient<std::ranges::views::all_t<R>, Grid, Value>;

template <std::ranges::viewable_range R, typename Grid,
          RealOrComplexValued Value>
Coefficient(R&&, Grid)
    -> Coefficient<std::ranges::views::all_t<R>, Grid, Value>;

// Type aliases for real and complex coefficients.
template <std::ranges::view View, typename Grid>
using RealCoefficient = Coefficient<View, Grid, RealValued>;

template <std::ranges::view View, typename Grid>
using ComplexCoefficient = Coefficient<View, Grid, ComplexValued>;

//-------------------------------------------------------//
//             Range adaptors for Coefficient            //
//-------------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormCoefficient : public std::ranges::range_adaptor_closure<
                            FormCoefficient<_Grid, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  FormCoefficient(_Grid grid, Int n) : _grid{grid}, _n{n} {}
  FormCoefficient(_Grid grid) : FormCoefficient(grid, 0) {}

  template <std::ranges::view _View>
  auto operator()(_View view) {
    return Coefficient<_View, _Grid, _Value>(view, _grid);
  }

 private:
  _Grid _grid;
  Int _n;
};

//---------------------------------------------------//
//             Coefficient -> Coefficient            //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(Coefficient<View, Grid, Value> u) {
  return u | std::ranges::views::transform([](auto x) { return -x; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

//---------------------------------------------------//
//        Coefficient x Complex -> Coefficient       //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(Coefficient<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<Value, ComplexValued>
auto operator*(typename Grid::Complex s, Coefficient<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(Coefficient<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
requires std::same_as<Value, ComplexValued>
auto operator+(typename Grid::Complex s, Coefficient<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(Coefficient<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator/(Coefficient<View, Grid, Value> u, typename Grid::Complex s)
requires std::same_as<Value, ComplexValued>
{
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

//---------------------------------------------------//
//          Coefficient x Real -> Coefficient        //
//---------------------------------------------------//
template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(Coefficient<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x * s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator*(typename Grid::Real s, Coefficient<View, Grid, Value> u) {
  return u * s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(Coefficient<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x + s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator+(typename Grid::Real s, Coefficient<View, Grid, Value> u) {
  return u + s;
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator-(Coefficient<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x - s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
auto operator/(Coefficient<View, Grid, Value> u, typename Grid::Real s) {
  return u | std::ranges::views::transform([s](auto x) { return x / s; }) |
         FormCoefficient<Grid, Value>(u.GetGrid(), u.AbsUpperIndex());
}

//---------------------------------------------------//
//      Coefficient x Coefficient -> Coefficient     //
//---------------------------------------------------//
template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value>
auto operator+(Coefficient<View1, Grid, Value> u1,
               Coefficient<View2, Grid, Value> u2) {
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x + y; },
                                           u1, u2) |
         FormCoefficient<Grid, Value>(u1.GetGrid(), u1.AbsUpperIndex());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value>
auto operator-(Coefficient<View1, Grid, Value> u1,
               Coefficient<View2, Grid, Value> u2) {
  assert(u1.size() == u2.size());
  return std::ranges::views::zip_transform([](auto x, auto y) { return x - y; },
                                           u1, u2) |
         FormCoefficient<Grid, Value>(u1.GetGrid(), u1.AbsUpperIndex());
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_COEFFICIENT_GUARD_H