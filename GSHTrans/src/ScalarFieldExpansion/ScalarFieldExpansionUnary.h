#ifndef GSH_TRANS_SCALAR_FIELD_EXPANSION_UNARY_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPANSION_UNARY_GUARD_H

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
  requires std::invocable<Function, typename _Derived::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Complex>,
      typename _Derived::Complex>;
}
class ScalarFieldExpansionUnary;

// Set traits.
namespace Internal {

template <typename _Derived, typename Function>
struct Traits<ScalarFieldExpansionUnary<_Derived, Function>> {
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
  requires std::invocable<Function, typename _Derived::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Complex>,
      typename _Derived::Complex>;
}
class ScalarFieldExpansionUnary
    : public ScalarFieldExpansionBase<
          ScalarFieldExpansionUnary<_Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldExpansionUnary<_Derived, Function>>::Int;
  using Value = typename Internal::Traits<
      ScalarFieldExpansionUnary<_Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldExpansionUnary<_Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldExpansionUnary<_Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldExpansionUnary<_Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldExpansionUnary<_Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarFieldExpansion Base.
  auto& Grid() const { return _u.Grid(); }
  auto operator[](Int l, Int m) const { return _f(_u[l, m]); }

  // Constructors.
  ScalarFieldExpansionUnary() = delete;
  ScalarFieldExpansionUnary(const ScalarFieldExpansionBase<_Derived>& u,
                            Function&& f)
      : _u{u}, _f{f} {}

  ScalarFieldExpansionUnary(const ScalarFieldExpansionUnary&) = default;
  ScalarFieldExpansionUnary(ScalarFieldExpansionUnary&&) = default;

  // Assignment.
  ScalarFieldExpansionUnary& operator=(ScalarFieldExpansionUnary&) = default;
  ScalarFieldExpansionUnary& operator=(ScalarFieldExpansionUnary&&) = default;

 private:
  const ScalarFieldExpansionBase<_Derived>& _u;
  Function& _f;
};

}  // namespace GSHTrans

#endif