#ifndef GSH_TRANS_ISOTROPIC_MATRIX_FIELD_GUARD_H
#define GSH_TRANS_ISOTROPIC_MATRIX_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "MatrixFieldBase.h"

namespace GSHTrans {

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class IsotropicMatrixField
    : public MatrixFieldBase<IsotropicMatrixField<_Grid, _Value>> {
 public:
  using Int =
      typename MatrixFieldBase<IsotropicMatrixField<_Grid, _Value>>::Int;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from MatrixField Base.
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
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->FieldSize(), a)} {}

  IsotropicMatrixField(_Grid grid) : IsotropicMatrixField(grid, 0) {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  IsotropicMatrixField(ScalarFieldBase<Derived>& u) : _grid{u.GetGrid()} {}

  // Assignment.
  IsotropicMatrixField& operator=(const IsotropicMatrixField&) = default;
  IsotropicMatrixField& operator=(IsotropicMatrixField&&) = default;

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_ISOTROPIC_MATRIX_FIELD_GUARD_H