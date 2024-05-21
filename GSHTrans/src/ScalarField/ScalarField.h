#ifndef GSH_TRANS_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_GUARD_H

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

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarField : public ScalarFieldBase<ScalarField<_Grid, _Value>> {
 public:
  using Int = typename ScalarFieldBase<ScalarField<_Grid, _Value>>::Int;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Constructors.
  ScalarField() = default;

  // Construct from grid initialising values to zero.
  ScalarField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->FieldSize())} {}

  // Construction from grid initialising values with a function.
  template <typename Function>
  requires std::invocable<Function, Real, Real>
  ScalarField(_Grid grid, Function&& f) : ScalarField(grid) {
    std::ranges::copy(_grid.InterpolateFunction(f), _data.begin());
  }

  // Construct from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(const ScalarFieldBase<Derived>& other)
      : ScalarField(other.GetGrid()) {
    assert(this->FieldSize() == other.FieldSize());
    CopyValues(other);
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(ScalarFieldBase<Derived>&& other) : ScalarField(other) {}

  // Default copy and move constructors.
  ScalarField(const ScalarField&) = default;
  ScalarField(ScalarField&&) = default;

  // Default copy and move assigment.
  ScalarField& operator=(const ScalarField&) = default;
  ScalarField& operator=(ScalarField&&) = default;

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
  FFTWpp::vector<Scalar> _data;

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

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealScalarField = ScalarField<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexScalarField = ScalarField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // namespace GSHTrans
