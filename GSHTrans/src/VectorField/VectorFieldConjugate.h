#ifndef GSH_TRANS_VECTOR_FIELD_CONJUGATE_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_CONJUGATE_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "VectorFieldBase.h"
#include "VectorFieldConstComponent.h"

namespace GSHTrans {

// Forward declare class.
template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, ComplexValued>;
}
class VectorFieldConjugate;

// Set traits.
namespace Internal {

template <typename Derived>
struct Traits<VectorFieldConjugate<Derived>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived::Grid;
  using Value = ComplexValued;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = Complex;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, ComplexValued>;
}
class VectorFieldConjugate
    : public VectorFieldBase<VectorFieldConjugate<Derived>> {
 public:
  using Int = typename Internal::Traits<VectorFieldConjugate<Derived>>::Int;
  using Grid = typename Internal::Traits<VectorFieldConjugate<Derived>>::Grid;
  using Value = typename Internal::Traits<VectorFieldConjugate<Derived>>::Value;
  using Real = typename Internal::Traits<VectorFieldConjugate<Derived>>::Real;
  using Complex =
      typename Internal::Traits<VectorFieldConjugate<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<VectorFieldConjugate<Derived>>::Scalar;
  using Writeable =
      typename Internal::Traits<VectorFieldConjugate<Derived>>::Writeable;

  // Return grid.
  auto GetGrid() const { return _u.GetGrid(); }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return MinusOneToPower<Real>(alpha) * std::conj(_u[-alpha, iTheta, iPhi]);
  }

  // Read access to component.
  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldConstComponent(*this, alpha);
  }

  // Constructors.
  VectorFieldConjugate() = delete;
  VectorFieldConjugate(const VectorFieldBase<Derived>& u) : _u{u} {}

  VectorFieldConjugate(const VectorFieldConjugate&) = default;
  VectorFieldConjugate(VectorFieldConjugate&&) = default;

  // Assignment.
  VectorFieldConjugate& operator=(const VectorFieldConjugate&) = default;
  VectorFieldConjugate& operator=(VectorFieldConjugate&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
};

}  // namespace GSHTrans

#endif