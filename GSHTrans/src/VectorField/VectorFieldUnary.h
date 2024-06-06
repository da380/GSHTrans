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
class VectorFieldUnary;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<VectorFieldUnary<Derived, Function>> {
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
class VectorFieldUnary
    : public VectorFieldBase<VectorFieldUnary<Derived, Function>> {
 public:
  using Int =
      typename Internal::Traits<VectorFieldUnary<Derived, Function>>::Int;
  using Grid =
      typename Internal::Traits<VectorFieldUnary<Derived, Function>>::Grid;
  using Value =
      typename Internal::Traits<VectorFieldUnary<Derived, Function>>::Value;
  using Real =
      typename Internal::Traits<VectorFieldUnary<Derived, Function>>::Real;
  using Complex =
      typename Internal::Traits<VectorFieldUnary<Derived, Function>>::Complex;
  using Scalar =
      typename Internal::Traits<VectorFieldUnary<Derived, Function>>::Scalar;
  using Writeable =
      typename Internal::Traits<VectorFieldUnary<Derived, Function>>::Writeable;

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
  VectorFieldUnary() = delete;
  VectorFieldUnary(const VectorFieldBase<Derived>& u, Function&& f)
      : _u{u}, _f{f} {}

  VectorFieldUnary(const VectorFieldUnary&) = default;
  VectorFieldUnary(VectorFieldUnary&&) = default;

  // Assignment.
  VectorFieldUnary& operator=(VectorFieldUnary&) = default;
  VectorFieldUnary& operator=(VectorFieldUnary&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
  Function& _f;
};

}  // namespace GSHTrans

#endif