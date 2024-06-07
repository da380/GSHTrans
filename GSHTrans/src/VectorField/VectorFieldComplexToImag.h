#ifndef GSH_TRANS_VECTOR_FIELD_COMPLEX_TO_IMAG_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_COMPLEX_TO_IMAG_GUARD_H

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
class VectorFieldComplexToImag;

// Set traits.
namespace Internal {

template <typename Derived>
struct Traits<VectorFieldComplexToImag<Derived>> {
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
class VectorFieldComplexToImag
    : public VectorFieldBase<VectorFieldComplexToImag<Derived>> {
 public:
  using Int = typename Internal::Traits<VectorFieldComplexToImag<Derived>>::Int;
  using Grid =
      typename Internal::Traits<VectorFieldComplexToImag<Derived>>::Grid;
  using Value =
      typename Internal::Traits<VectorFieldComplexToImag<Derived>>::Value;
  using Real =
      typename Internal::Traits<VectorFieldComplexToImag<Derived>>::Real;
  using Complex =
      typename Internal::Traits<VectorFieldComplexToImag<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<VectorFieldComplexToImag<Derived>>::Scalar;
  using Writeable =
      typename Internal::Traits<VectorFieldComplexToImag<Derived>>::Writeable;

  // Return grid.
  auto GetGrid() const { return _u.GetGrid(); }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
    auto u = half * (_u[alpha, iTheta, iPhi] -
                     MinusOneToPower<Real>(alpha) *
                         std::conj(_u[-alpha, iTheta, iPhi]));
    return alpha < 0 ? -std::real(u) : std::imag(u);
  }

  // Read access to component.
  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldConstComponent(*this, alpha);
  }

  // Constructors.
  VectorFieldComplexToImag() = delete;
  VectorFieldComplexToImag(const VectorFieldBase<Derived>& u) : _u{u} {}

  VectorFieldComplexToImag(const VectorFieldComplexToImag&) = default;
  VectorFieldComplexToImag(VectorFieldComplexToImag&&) = default;

  // Assignment.
  VectorFieldComplexToImag& operator=(const VectorFieldComplexToImag&) =
      default;
  VectorFieldComplexToImag& operator=(VectorFieldComplexToImag&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
};

}  // namespace GSHTrans

#endif