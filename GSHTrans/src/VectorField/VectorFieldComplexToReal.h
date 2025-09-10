#ifndef GSH_TRANS_VECTOR_FIELD_COMPLEX_TO_REAL_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_COMPLEX_TO_REAL_GUARD_H

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
class VectorFieldComplexToReal;

// Set traits.
namespace Internal {

template <typename Derived>
struct Traits<VectorFieldComplexToReal<Derived>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived::Grid;
  using Value = RealValued;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = Real;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, ComplexValued>;
}
class VectorFieldComplexToReal
    : public VectorFieldBase<VectorFieldComplexToReal<Derived>> {
 public:
  using Int = typename Internal::Traits<VectorFieldComplexToReal<Derived>>::Int;
  using Grid =
      typename Internal::Traits<VectorFieldComplexToReal<Derived>>::Grid;
  using Value =
      typename Internal::Traits<VectorFieldComplexToReal<Derived>>::Value;
  using Real =
      typename Internal::Traits<VectorFieldComplexToReal<Derived>>::Real;
  using Complex =
      typename Internal::Traits<VectorFieldComplexToReal<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<VectorFieldComplexToReal<Derived>>::Scalar;
  using Writeable =
      typename Internal::Traits<VectorFieldComplexToReal<Derived>>::Writeable;

  // Return grid.
  auto GetGrid() const { return _u.GetGrid(); }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
    auto u = half * (_u[alpha, iTheta, iPhi] +
                     MinusOneToPower<Real>(alpha) *
                         std::conj(_u[-alpha, iTheta, iPhi]));
    return alpha < 0 ? std::imag(u) : std::real(u);
  }

  // Read access to component.
  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldConstComponent(*this, alpha);
  }

  // Constructors.
  VectorFieldComplexToReal() = delete;
  VectorFieldComplexToReal(const VectorFieldBase<Derived>& u) : _u{u} {}

  VectorFieldComplexToReal(const VectorFieldComplexToReal&) = default;
  VectorFieldComplexToReal(VectorFieldComplexToReal&&) = default;

  // Assignment.
  VectorFieldComplexToReal& operator=(const VectorFieldComplexToReal&) =
      default;
  VectorFieldComplexToReal& operator=(VectorFieldComplexToReal&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
};

}  // namespace GSHTrans

#endif