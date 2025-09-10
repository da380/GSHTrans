#ifndef GSH_TRANS_SCALAR_FIELD_EXPANSION_BINARY_WITH_SCALAR_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPANSION_BINARY_WITH_SCALAR_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldExpansionBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Derived, typename Function>
requires requires() {
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Complex,
                           typename _Derived::Complex>,
      typename _Derived::Complex>;
}
class ScalarFieldExpansionBinaryWithScalar;

// Set traits.
namespace Internal {

template <typename _Derived, typename Function>
struct Traits<ScalarFieldExpansionBinaryWithScalar<_Derived, Function>> {
  using Int = std::ptrdiff_t;
  using Real = typename _Derived::Real;
  using Complex = typename _Derived::Complex;
  using Value = typename _Derived::Value;
  using Scalar = typename _Derived::Scalar;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename _Derived, typename Function>
requires requires() {
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Complex,
                           typename _Derived::Complex>,
      typename _Derived::Complex>;
}
class ScalarFieldExpansionBinaryWithScalar
    : public ScalarFieldExpansionBase<
          ScalarFieldExpansionBinaryWithScalar<_Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<_Derived, Function>>::Int;
  using Value = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<_Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<_Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<_Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<_Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<_Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarFieldExpansion Base.
  auto& Grid() const { return _u.Grid(); }
  auto operator[](Int l, Int m) const { return _f(_u[l, m], _s); }

  // Constructors.
  ScalarFieldExpansionBinaryWithScalar() = delete;
  ScalarFieldExpansionBinaryWithScalar(
      const ScalarFieldExpansionBase<_Derived>& u, Scalar s, Function&& f)
      : _u{u}, _s{s}, _f{f} {}

  ScalarFieldExpansionBinaryWithScalar(
      const ScalarFieldExpansionBinaryWithScalar&) = default;
  ScalarFieldExpansionBinaryWithScalar(ScalarFieldExpansionBinaryWithScalar&&) =
      default;

  // Assignment.
  ScalarFieldExpansionBinaryWithScalar& operator=(
      ScalarFieldExpansionBinaryWithScalar&) = default;
  ScalarFieldExpansionBinaryWithScalar& operator=(
      ScalarFieldExpansionBinaryWithScalar&&) = default;

 private:
  const ScalarFieldExpansionBase<_Derived>& _u;
  Scalar _s;
  Function& _f;
};

}  // namespace GSHTrans

#endif