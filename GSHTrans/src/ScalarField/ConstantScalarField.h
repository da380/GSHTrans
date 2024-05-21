#ifndef GSH_TRANS_CONSTANT_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_CONSTANT_SCALAR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

//-------------------------------------------------//
//              Constant scalar field              //
//-------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ConstantScalarField
    : public ScalarFieldBase<ConstantScalarField<_Grid, _Value>> {
 public:
  using Int = typename ScalarFieldBase<ConstantScalarField<_Grid, _Value>>::Int;
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
    return _s;
  }

  // Constructors.
  ConstantScalarField() = default;

  ConstantScalarField(_Grid grid, Scalar s) : _grid{grid}, _s{s} {}

  // Assignment.
  ConstantScalarField& operator=(const ConstantScalarField&) = default;
  ConstantScalarField& operator=(ConstantScalarField&&) = default;

  // Compound assignment.
  auto& operator==(Scalar s) {
    _s = s;
    return *this;
  }

  auto& operator+=(Scalar s) {
    _s += s;
    return *this;
  }

  auto& operator-=(Scalar s) {
    _s -= s;
    return *this;
  }

  auto& operator*=(Scalar s) {
    _s *= s;
    return *this;
  }

 private:
  _Grid _grid;
  Scalar _s;
};

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealConstantScalarField = ConstantScalarField<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexConstantScalarField = ConstantScalarField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_CONSTANT_SCALAR_FIELD_GUARD_H