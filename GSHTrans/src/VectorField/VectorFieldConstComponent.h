#ifndef GSH_TRANS_VECTOR_FIELD_CONST_COMPONENT_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_CONST_COMPONENT_GUARD_H

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
class VectorFieldConstComponent;

// Set traits for the component class.
namespace Internal {

template <typename Derived>
struct Traits<VectorFieldConstComponent<Derived>> {
  using Int = typename Derived::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

// Define the component class.
template <typename Derived>
class VectorFieldConstComponent
    : public ScalarFieldBase<VectorFieldConstComponent<Derived>> {
 public:
  using Grid =
      typename Internal::Traits<VectorFieldConstComponent<Derived>>::Grid;
  using Value =
      typename Internal::Traits<VectorFieldConstComponent<Derived>>::Value;
  using Int =
      typename Internal::Traits<VectorFieldConstComponent<Derived>>::Int;
  using Real =
      typename Internal::Traits<VectorFieldConstComponent<Derived>>::Real;
  using Complex =
      typename Internal::Traits<VectorFieldConstComponent<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<VectorFieldConstComponent<Derived>>::Scalar;
  using Writeable =
      typename Internal::Traits<VectorFieldConstComponent<Derived>>::Writeable;

  // Return the grid.
  auto GetGrid() const { return _v.GetGrid(); }

  // Read access to the data.
  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _v[_alpha, iTheta, iPhi];
  }

  // Constructors.
  VectorFieldConstComponent() = default;

  // Construct from grid initialising values to zero.
  VectorFieldConstComponent(const VectorFieldBase<Derived>& v, Int alpha)
      : _v{v}, _alpha{alpha} {}

  // Default copy and move constructors.
  VectorFieldConstComponent(const VectorFieldConstComponent&) = default;
  VectorFieldConstComponent(VectorFieldConstComponent&&) = default;

  // Default copy and move assigment.
  VectorFieldConstComponent& operator=(const VectorFieldConstComponent&) =
      default;
  VectorFieldConstComponent& operator=(VectorFieldConstComponent&&) = default;

  // Use assignment defined in base class.
  using ScalarFieldBase<VectorFieldConstComponent<Derived>>::operator=;

 private:
  Int _alpha;
  const VectorFieldBase<Derived>& _v;
};

}  // namespace GSHTrans

#endif
