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
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Complex>,
      typename Derived::Complex>;
}
class ScalarFieldExpansionUnary;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<ScalarFieldExpansionUnary<Derived, Function>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = typename Derived::Value;
  using Scalar = typename Derived::Scalar;
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
class ScalarFieldExpansionUnary
    : public ScalarFieldExpansionBase<
          ScalarFieldExpansionUnary<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldExpansionUnary<Derived, Function>>::Int;
  using Grid = typename Internal::Traits<
      ScalarFieldExpansionUnary<Derived, Function>>::Grid;
  using Value = typename Internal::Traits<
      ScalarFieldExpansionUnary<Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldExpansionUnary<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldExpansionUnary<Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldExpansionUnary<Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldExpansionUnary<Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarFieldExpansion Base.
  auto& GetGrid() const { return _u.GetGrid(); }
  auto operator[](Int l, Int m) const { return _f(_u[l, m]); }

  // Constructors.
  ScalarFieldExpansionUnary() = delete;
  ScalarFieldExpansionUnary(const ScalarFieldExpansionBase<Derived>& u,
                            Function&& f)
      : _u{u}, _f{f} {}

  ScalarFieldExpansionUnary(const ScalarFieldExpansionUnary&) = default;
  ScalarFieldExpansionUnary(ScalarFieldExpansionUnary&&) = default;

  // Assignment.
  ScalarFieldExpansionUnary& operator=(ScalarFieldExpansionUnary&) = default;
  ScalarFieldExpansionUnary& operator=(ScalarFieldExpansionUnary&&) = default;

 private:
  const ScalarFieldExpansionBase<Derived>& _u;
  Function& _f;
};

}  // namespace GSHTrans

#endif