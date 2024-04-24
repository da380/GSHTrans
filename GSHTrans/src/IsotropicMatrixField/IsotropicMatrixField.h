#ifndef GSH_TRANS_ISOTROPIC_MATRIX_FIELD_GUARD_H
#define GSH_TRANS_ISOTROPIC_MATRIX_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "../MatrixField/MatrixFieldBase.h"
#include "IsotropicMatrixFieldBase.h"

namespace GSHTrans {

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class IsotropicMatrixField
    : public IsotropicMatrixFieldBase<IsotropicMatrixField<_Grid, _Value>> {
 public:
  using Int = typename IsotropicMatrixFieldBase<
      IsotropicMatrixField<_Grid, _Value>>::Int;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from IsotropicMatrixField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    return alpha + beta == 0 ? static_cast<Scalar>(MinusOneToPower(alpha)) *
                                   _data[Index(iTheta, iPhi)]
                             : Scalar{0};
  }

  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  IsotropicMatrixField() = default;

  IsotropicMatrixField(_Grid grid, Scalar a)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>((this->FieldSize()), a)} {}

  IsotropicMatrixField(_Grid grid) : IsotropicMatrixField(grid, 0) {}

  // Assignment.
  IsotropicMatrixField& operator=(const IsotropicMatrixField&) = default;
  IsotropicMatrixField& operator=(IsotropicMatrixField&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const IsotropicMatrixFieldBase<Derived>& other) {
    assert(this->FieldSize() == other.FieldSize());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(IsotropicMatrixFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Value assignment.
  auto& operator()(Int alpha, Int beta, Int iTheta, Int iPhi) {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }
  /*
  auto operator()(Int alpha, Int beta) {
    this->CheckCanonicalIndices(alpha, beta);
    auto start = std::next(_data.begin(), Offset(alpha, beta));
    auto finish = std::next(start, this->FieldSize());
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }
  */

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  template <typename Derived>
  void CopyValues(const IsotropicMatrixFieldBase<Derived>& other) {}

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

// Type aliases for real and complex Matrix fields
template <typename Grid>
using RealIsotropicMatrixField = IsotropicMatrixField<Grid, RealValued>;

template <typename Grid>
using ComplexIsotropicMatrixField = IsotropicMatrixField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_ISOTROPIC_MATRIX_FIELD_GUARD_H