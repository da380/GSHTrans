#ifndef GSH_TRANS_ABSTRACT_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_ABSTRACT_SCALAR_FIELD_GUARD_H

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
template <typename _Grid, RealOrComplexValued _Value, typename _Function>
requires std::derived_from<_Grid, GridBase<_Grid>> &&
         ScalarFunctionS2<
             _Function, typename _Grid::Real,
             std::conditional_t<std::same_as<_Value, RealValued>,
                                typename _Grid::Real, typename _Grid::Complex>>
class AbstractScalarField;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value, typename _Function>
struct Traits<AbstractScalarField<_Grid, _Value, _Function>> {
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
template <typename _Grid, RealOrComplexValued _Value, typename _Function>
requires std::derived_from<_Grid, GridBase<_Grid>> &&
         ScalarFunctionS2<
             _Function, typename _Grid::Real,
             std::conditional_t<std::same_as<_Value, RealValued>,
                                typename _Grid::Real, typename _Grid::Complex>>
class AbstractScalarField
    : public ScalarFieldBase<AbstractScalarField<_Grid, _Value, _Function>> {
 public:
  using Grid = typename Internal::Traits<
      AbstractScalarField<_Grid, _Value, _Function>>::Grid;
  using Value = typename Internal::Traits<
      AbstractScalarField<_Grid, _Value, _Function>>::Value;
  using Int = typename Internal::Traits<
      AbstractScalarField<_Grid, _Value, _Function>>::Int;
  using Real = typename Internal::Traits<
      AbstractScalarField<_Grid, _Value, _Function>>::Real;
  using Complex = typename Internal::Traits<
      AbstractScalarField<_Grid, _Value, _Function>>::Complex;
  using Scalar = typename Internal::Traits<
      AbstractScalarField<_Grid, _Value, _Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      AbstractScalarField<_Grid, _Value, _Function>>::Writeable;

  // Return the grid.
  auto& GetGrid() const { return _grid; }

  // Read access to data.
  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    auto theta = this->CoLatitudes()[iTheta];
    auto phi = this->Longitudes()[iPhi];
    return _f(theta, phi);
  }

  // Default constructor.
  AbstractScalarField() = default;

  // Construct from grid initialising values to zero.
  AbstractScalarField(_Grid& grid, _Function&& f) : _grid{grid}, _f{f} {}

  // Default copy and move constructors.
  AbstractScalarField(const AbstractScalarField&) = default;
  AbstractScalarField(AbstractScalarField&&) = default;

  // Default copy and move assigment.
  AbstractScalarField& operator=(const AbstractScalarField&) = default;
  AbstractScalarField& operator=(AbstractScalarField&&) = default;

 private:
  _Grid& _grid;
  _Function& _f;
};  // namespace GSHTrans

// Type aliases for real and complex fields.
template <typename Grid, typename Function>
using RealAbstractScalarField = AbstractScalarField<Grid, RealValued, Function>;

template <typename Grid, typename Function>
using ComplexAbstractScalarField =
    AbstractScalarField<Grid, ComplexValued, Function>;

}  // namespace GSHTrans

#endif
