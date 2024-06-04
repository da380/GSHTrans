#ifndef GSH_TRANS_CONSTANT_VECTOR_FIELD_GUARD_H
#define GSH_TRANS_CONSTANT_VECTOR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"
#include "VectorFieldBase.h"
#include "VectorFieldConstComponent.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ConstantVectorField;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value>
struct Traits<ConstantVectorField<_Grid, _Value>> {
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
class ConstantVectorField
    : public VectorFieldBase<ConstantVectorField<_Grid, _Value>> {
 public:
  using Grid =
      typename Internal::Traits<ConstantVectorField<_Grid, _Value>>::Grid;
  using Value =
      typename Internal::Traits<ConstantVectorField<_Grid, _Value>>::Value;
  using Int =
      typename Internal::Traits<ConstantVectorField<_Grid, _Value>>::Int;
  using Real =
      typename Internal::Traits<ConstantVectorField<_Grid, _Value>>::Real;
  using Complex =
      typename Internal::Traits<ConstantVectorField<_Grid, _Value>>::Complex;
  using Scalar =
      typename Internal::Traits<ConstantVectorField<_Grid, _Value>>::Scalar;
  using Writeable =
      typename Internal::Traits<ConstantVectorField<_Grid, _Value>>::Writeable;

  // Return the grid.
  auto GetGrid() const { return _grid; }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const {
    this->CheckCanonicalIndices(alpha);
    this->CheckPointIndices(iTheta, iPhi);
    return _v[alpha];
  }

  // Return read access component
  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return ConstantVectorFieldConstComponent(*this, alpha);
  }

  // Default constructor.
  ConstantVectorField() = default;

  // Construct from grid initialising with Canonical vector.
  ConstantVectorField(_Grid grid, CanonicalVector<Scalar>&& v)
      : _grid{grid}, _v{v} {}

  // Construct from grid initialising with Canonical components.
  ConstantVectorField(_Grid grid, Scalar m, Scalar z, Scalar p)
      : _grid{grid}, _v{CanonicalVector{m, z, p}} {}

  // Default copy and move constructors.
  ConstantVectorField(const ConstantVectorField&) = default;
  ConstantVectorField(ConstantVectorField&&) = default;

  // Default copy and move assigment.
  ConstantVectorField& operator=(const ConstantVectorField&) = default;
  ConstantVectorField& operator=(ConstantVectorField&&) = default;

 private:
  _Grid _grid;
  CanonicalVector<Scalar> _v;
};

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealConstantVectorField = ConstantVectorField<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexConstantVectorField = ConstantVectorField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_GUARD_H
