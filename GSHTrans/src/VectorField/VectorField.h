#ifndef GSH_TRANS_VECTOR_FIELD_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "VectorFieldBase.h"

namespace GSHTrans {

//-------------------------------------------------//
//       Vector field that stores its data         //
//-------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class VectorField : public VectorFieldBase<VectorField<_Grid, _Value>> {
 public:
  using Int = typename VectorFieldBase<VectorField<_Grid, _Value>>::Int;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorField() = default;

  VectorField(_Grid grid)
      : _grid{grid},
        _data{FFTWpp::vector<Scalar>(3 * (this->ComponentSize()))} {}

  VectorField(_Grid grid, std::array<Scalar, 3>&& u) : VectorField(grid) {
    for (auto [i, alpha] :
         std::ranges::views::enumerate(this->CanonicalIndices())) {
      this->operator()(alpha) = u[i];
    }
  }

  // Assignment.
  VectorField& operator=(const VectorField&) = default;
  VectorField& operator=(VectorField&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const VectorFieldBase<Derived>& other) {
    assert(this->ComponentSize() == other.ComponentSize());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(VectorFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Methods to make it a range.
  auto size() const { return this->Size(); }
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignment.
  auto& operator()(Int alpha, Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) {
    this->CheckCanonicalIndices(alpha);
    auto start = std::next(begin(), Offset(alpha));
    auto finish = std::next(start, this->ComponentSize());
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  template <typename Derived>
  void CopyValues(const VectorFieldBase<Derived>& other) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator()(alpha, iTheta, iPhi) = other(alpha, iTheta, iPhi);
      }
    }
  }

  auto Offset(Int alpha) const { return (alpha + 1) * (this->ComponentSize()); }

  auto Index(Int alpha, Int iTheta, int iPhi) const {
    return Offset(alpha) + iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

// Type aliases for real and complex vector fields
template <typename Grid>
using RealVectorField = VectorField<Grid, RealValued>;

template <typename Grid>
using ComplexVectorField = VectorField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_GUARD_H