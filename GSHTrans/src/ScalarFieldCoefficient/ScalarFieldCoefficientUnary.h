#ifndef GSH_TRANS_SCALAR_FIELD_COEFFICIENT_UNARY_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_COEFFICIENT_UNARY_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldCoefficientBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Complex>,
      typename Derived::Complex>;
}
class ScalarFieldCoefficientUnary;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<ScalarFieldCoefficientUnary<Derived, Function>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Complex>,
      typename Derived::Complex>;
}
class ScalarFieldCoefficientUnary
    : public ScalarFieldCoefficientBase<
          ScalarFieldCoefficientUnary<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldCoefficientUnary<Derived, Function>>::Int;
  using Grid = typename Internal::Traits<
      ScalarFieldCoefficientUnary<Derived, Function>>::Grid;
  using Value = typename Internal::Traits<
      ScalarFieldCoefficientUnary<Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldCoefficientUnary<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldCoefficientUnary<Derived, Function>>::Complex;
  using Writeable = typename Internal::Traits<
      ScalarFieldCoefficientUnary<Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarFieldCoefficient Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator[](Int l, Int m) const { return _f(_u[l, m]); }

  // Constructors.
  ScalarFieldCoefficientUnary() = delete;
  ScalarFieldCoefficientUnary(const ScalarFieldCoefficientBase<Derived>& u,
                              Function&& f)
      : _u{u}, _f{f} {}

  ScalarFieldCoefficientUnary(const ScalarFieldCoefficientUnary&) = default;
  ScalarFieldCoefficientUnary(ScalarFieldCoefficientUnary&&) = default;

  // Assignment.
  ScalarFieldCoefficientUnary& operator=(ScalarFieldCoefficientUnary&) =
      default;
  ScalarFieldCoefficientUnary& operator=(ScalarFieldCoefficientUnary&&) =
      default;

 private:
  const ScalarFieldCoefficientBase<Derived>& _u;
  Function& _f;
};

}  // namespace GSHTrans

#endif