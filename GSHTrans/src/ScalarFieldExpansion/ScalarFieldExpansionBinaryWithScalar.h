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
template <typename Derived, typename Function>
requires requires() {
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Complex,
                           typename Derived::Complex>,
      typename Derived::Complex>;
}
class ScalarFieldExpansionBinaryWithScalar;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<ScalarFieldExpansionBinaryWithScalar<Derived, Function>> {
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
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Complex,
                           typename Derived::Complex>,
      typename Derived::Complex>;
}
class ScalarFieldExpansionBinaryWithScalar
    : public ScalarFieldExpansionBase<
          ScalarFieldExpansionBinaryWithScalar<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<Derived, Function>>::Int;
  using Grid = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<Derived, Function>>::Grid;
  using Value = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldExpansionBinaryWithScalar<Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarFieldExpansion Base.
  auto& GetGrid() const { return _u.GetGrid(); }
  auto operator[](Int l, Int m) const { return _f(_u[l, m], _s); }

  // Constructors.
  ScalarFieldExpansionBinaryWithScalar() = delete;
  ScalarFieldExpansionBinaryWithScalar(
      const ScalarFieldExpansionBase<Derived>& u, Scalar s, Function&& f)
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
  const ScalarFieldExpansionBase<Derived>& _u;
  Scalar _s;
  Function& _f;
};

}  // namespace GSHTrans

#endif