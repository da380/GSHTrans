#ifndef GSH_TRANS_VECTOR_FIELD_UNARY_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_UNARY_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "VectorFieldBase.h"
#include "VectorFieldConstComponent.h"

namespace GSHTrans {

// Forward declare class.
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class VectorFieldPointwiseUnary;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<VectorFieldPointwiseUnary<Derived, Function>> {
  using Int = typename Derived::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class VectorFieldPointwiseUnary
    : public VectorFieldBase<VectorFieldPointwiseUnary<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      VectorFieldPointwiseUnary<Derived, Function>>::Int;
  using Grid = typename Internal::Traits<
      VectorFieldPointwiseUnary<Derived, Function>>::Grid;
  using Value = typename Internal::Traits<
      VectorFieldPointwiseUnary<Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      VectorFieldPointwiseUnary<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      VectorFieldPointwiseUnary<Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      VectorFieldPointwiseUnary<Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      VectorFieldPointwiseUnary<Derived, Function>>::Writeable;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _u.GetGrid(); }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const -> Scalar {
    this->CheckCanonicalIndices(alpha);
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[alpha, iTheta, iPhi]);
  }

  // Read access component
  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldConstComponent(*this, alpha);
  }

  // Constructors.
  VectorFieldPointwiseUnary() = delete;
  VectorFieldPointwiseUnary(const VectorFieldBase<Derived>& u, Function f)
      : _u{u}, _f{f} {}

  VectorFieldPointwiseUnary(const VectorFieldPointwiseUnary&) = default;
  VectorFieldPointwiseUnary(VectorFieldPointwiseUnary&&) = default;

  // Assignment.
  VectorFieldPointwiseUnary& operator=(VectorFieldPointwiseUnary&) = default;
  VectorFieldPointwiseUnary& operator=(VectorFieldPointwiseUnary&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
  Function _f;
};

}  // namespace GSHTrans

#endif