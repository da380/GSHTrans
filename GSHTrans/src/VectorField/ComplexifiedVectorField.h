#ifndef GSH_TRANS_COMPLEXIFIED_VECTOR_FIELD_GUARD_H
#define GSH_TRANS_COMPLEXIFIED_VECTOR_FIELD_GUARD_H

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
template <typename Derived>
requires std::same_as<RealValued, typename Derived::Value>
class ComplexifiedVectorField;

// Set traits.
namespace Internal {

template <typename Derived>
struct Traits<ComplexifiedVectorField<Derived>> {
  using Int = typename Derived::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = ComplexValued;
  using Scalar = Complex;
  using Writeable = std::false_type;
};

}  // namespace Internal

}  // namespace GSHTrans

#endif  // GSH_TRANS_COMPLEXIFIED_VECTOR_FIELD_GUARD_H