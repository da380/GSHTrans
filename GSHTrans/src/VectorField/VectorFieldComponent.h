#ifndef GSH_TRANS_VECTOR_FIELD_COMPONENT_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_COMPONENT_GUARD_H

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

namespace GSHTrans {

// Forward declare the component class.
template <typename Derived>
requires Derived::Writeable::value
class VectorFieldComponent;

// Set traits for the component class.
namespace Internal {

template <typename Derived>
struct Traits<VectorFieldComponent<Derived>> {
  using Int = typename Derived::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::true_type;
};

}  // namespace Internal

// Define the component class.
template <typename Derived>
requires Derived::Writeable::value
class VectorFieldComponent
    : public ScalarFieldBase<VectorFieldComponent<Derived>> {
 public:
  using Grid = typename Internal::Traits<VectorFieldComponent<Derived>>::Grid;
  using Value = typename Internal::Traits<VectorFieldComponent<Derived>>::Value;
  using Int = typename Internal::Traits<VectorFieldComponent<Derived>>::Int;
  using Real = typename Internal::Traits<VectorFieldComponent<Derived>>::Real;
  using Complex =
      typename Internal::Traits<VectorFieldComponent<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<VectorFieldComponent<Derived>>::Scalar;
  using Writeable =
      typename Internal::Traits<VectorFieldComponent<Derived>>::Writeable;

  // Return the grid.
  auto GetGrid() const { return _v.GetGrid(); }

  // Read access to the data.
  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _v[_alpha, iTheta, iPhi];
  }

  // Write access to the data.
  auto& operator[](Int iTheta, Int iPhi)
  requires Writeable::value
  {
    this->CheckPointIndices(iTheta, iPhi);
    return _v[_alpha, iTheta, iPhi];
  }

  // Constructors.
  VectorFieldComponent() = default;

  // Construct from grid initialising values to zero.
  VectorFieldComponent(VectorFieldBase<Derived>& v, Int alpha)
      : _v{v}, _alpha{alpha} {}

  // Default copy and move constructors.
  VectorFieldComponent(const VectorFieldComponent&) = default;
  VectorFieldComponent(VectorFieldComponent&&) = default;

  // Default copy and move assigment.
  VectorFieldComponent& operator=(const VectorFieldComponent&) = default;
  VectorFieldComponent& operator=(VectorFieldComponent&&) = default;

  // Use assignment defined in base class.
  using ScalarFieldBase<VectorFieldComponent<Derived>>::operator=;

 private:
  Int _alpha;
  VectorFieldBase<Derived>& _v;
};

}  // namespace GSHTrans

#endif  // #ifndef GSH_TRANS_SCALAR_FIELD_GUARD_H
