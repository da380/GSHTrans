#ifndef GSH_TRANS_CONSTANT_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_CONSTANT_SCALAR_FIELD_GUARD_H

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

// Forward declare class.
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ConstantScalarField;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value>
struct Traits<ConstantScalarField<_Grid, _Value>> {
  using Int = std::ptrdiff_t;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<Value, RealValued>, Real, Complex>;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ConstantScalarField
    : public ScalarFieldBase<ConstantScalarField<_Grid, _Value>> {
 public:
  using Int =
      typename Internal::Traits<ConstantScalarField<_Grid, _Value>>::Int;
  using Grid =
      typename Internal::Traits<ConstantScalarField<_Grid, _Value>>::Grid;
  using Value =
      typename Internal::Traits<ConstantScalarField<_Grid, _Value>>::Value;
  using Real =
      typename Internal::Traits<ConstantScalarField<_Grid, _Value>>::Real;
  using Complex =
      typename Internal::Traits<ConstantScalarField<_Grid, _Value>>::Complex;
  using Scalar =
      typename Internal::Traits<ConstantScalarField<_Grid, _Value>>::Scalar;
  using Writeable =
      typename Internal::Traits<ConstantScalarField<_Grid, _Value>>::Writeable;

  // Return the grid.
  auto GetGrid() const { return _grid; }

  // Read access to data.
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _s;
  }

  // Default constructor.
  ConstantScalarField() = default;

  // Construct from grid initialising to constant value.
  ConstantScalarField(_Grid grid, Scalar s) : _grid{grid}, _s{s} {}

  // Default copy and move constructors.
  ConstantScalarField(const ConstantScalarField&) = default;
  ConstantScalarField(ConstantScalarField&&) = default;

  // Default copy and move assigment.
  ConstantScalarField& operator=(const ConstantScalarField&) = default;
  ConstantScalarField& operator=(ConstantScalarField&&) = default;

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
