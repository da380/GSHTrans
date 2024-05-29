#ifndef GSH_TRANS_SCALAR_FIELD_VIEW_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_VIEW_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

template <typename _Grid, std::ranges::view _View>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarFieldView : public ScalarFieldBase<ScalarFieldView<_Grid, _View>> {
 public:
  using Int = typename ScalarFieldBase<ScalarFieldView<_Grid, _View>>::Int;
  using Grid = _Grid;
  using Scalar = std::ranges::range_value_t<_View>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // Methods needed to inherit from ScalarFieldView Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Constructors.
  ScalarFieldView() = default;

  // Construct from grid initialising values to zero.
  ScalarFieldView(_Grid grid, _View data) : _grid{grid}, _data{data} {}

  // Construction from grid initialising values with a function.
  template <typename Function>
  requires std::invocable<Function, Real, Real>
  ScalarFieldView(_Grid grid, Function&& f, _View data)
      : ScalarFieldView(grid, data) {
    std::ranges::copy(_grid.InterpolateFunction(f), _data.begin());
  }

  // Default copy and move constructors.
  ScalarFieldView(const ScalarFieldView&) = default;
  ScalarFieldView(ScalarFieldView&&) = default;

  // Default copy and move assigment.
  ScalarFieldView& operator=(const ScalarFieldView&) = default;
  ScalarFieldView& operator=(ScalarFieldView&&) = default;

  // Assign values from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(const ScalarFieldBase<Derived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(ScalarFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Compound assigments.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator+=(const ScalarFieldBase<Derived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      _data[Index(iTheta, iPhi)] += other(iTheta, iPhi);
    }
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator+=(ScalarFieldBase<Derived>&& other) {
    *this += other;
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator-=(const ScalarFieldBase<Derived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      _data[Index(iTheta, iPhi)] -= other(iTheta, iPhi);
    }
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator-=(ScalarFieldBase<Derived>&& other) {
    *this -= other;
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator*=(const ScalarFieldBase<Derived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      _data[Index(iTheta, iPhi)] *= other(iTheta, iPhi);
    }
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator*=(ScalarFieldBase<Derived>&& other) {
    *this *= other;
    return *this;
  }

  // Return view to the data.
  auto Data() { return _data; }

  // Index level assignement and increment.
  void Set(Int iTheta, Int iPhi, Scalar s) {
    this->CheckPointIndices(iTheta, iPhi);
    _data[Index(iTheta, iPhi)] = s;
  }

  void Add(Int iTheta, Int iPhi, Scalar s) {
    this->CheckPointIndices(iTheta, iPhi);
    _data[Index(iTheta, iPhi)] += s;
  }

 private:
  _Grid _grid;
  _View _data;

  template <typename Derived>
  void CopyValues(const ScalarFieldBase<Derived>& other) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      Set(iTheta, iPhi, other(iTheta, iPhi));
    }
  }

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

// Deduction guide to allow construction from a range.
template <typename Grid, std::ranges::viewable_range R>
ScalarFieldView(Grid,
                R&&) -> ScalarFieldView<Grid, std::ranges::views::all_t<R>>;

template <typename Grid, typename Function, std::ranges::viewable_range R>
ScalarFieldView(Grid, Function&&,
                R&&) -> ScalarFieldView<Grid, std::ranges::views::all_t<R>>;

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealScalarFieldView = ScalarFieldView<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexScalarFieldView = ScalarFieldView<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // #ifndef GSH_TRANS_SCALAR_FIELD_GUARD_H
